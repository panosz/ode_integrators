//
// Created by Panagiotis Zestanakis on 05/05/18.
//

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

#include <boost/numeric/odeint.hpp>
#include "gnuplot-iostream.h"
#include "utilities.hpp"

using namespace boost::numeric::odeint;

using ActionAngle = std::complex<double>;
using State = std::array<ActionAngle, 2>;

class KarneySystem {
 private:
  const double a_;  //coupling constant
  const double nu_; //relative frequency
 public:
  KarneySystem (double coupling_constant, double relative_frequency)
      :
      a_{coupling_constant}, nu_{relative_frequency}
  { };

  virtual void operator() (const State& qJ, State& dqJdt, const double /* t */) const
  {
    const auto &[qJ1, qJ2]= qJ;
    const auto q1 = qJ1.real();
    const auto J1 = qJ1.imag();

    const auto q2 = qJ2.real();
    //const auto J2 =qJ2.imag();

    const double r1 = std::sqrt(2 * J1);

    const double sin_q1 = std::sin(q1);
    const double cos_q1 = std::cos(q1);

    const double rel_phase = r1 * sin_q1 - q2;
    //const double sin_rel_phase = std::sin(rel_phase);
    const double cos_rel_phase = std::cos(rel_phase);

    const double dq1dt = 1 - a_ / r1 * sin_q1 * cos_rel_phase;
    const double dJ1dt = a_ * r1 * cos_q1 * cos_rel_phase;
    const double dq2dt = nu_;
    const double dJ2dt = -a_ * cos_rel_phase;

    dqJdt[0] = ActionAngle{dq1dt, dJ1dt};
    dqJdt[1] = ActionAngle{dq2dt, dJ2dt};
  }

  double hamiltonian (const State& qJ)
  {
    const auto &[qJ1, qJ2]= qJ;
    const auto q1 = qJ1.real();
    const auto J1 = qJ1.imag();

    const auto q2 = qJ2.real();
    const auto J2 = qJ2.imag();

    const double r1 = std::sqrt(2 * J1);

    const double sin_q1 = std::sin(q1);

    const double rel_phase = r1 * sin_q1 - q2;

    return J1 + nu_ * J2 - a_ * std::sin(rel_phase);
  }

  virtual ~KarneySystem () = default;
};

class KarneySystemWithq1AsIndependentVariable : public KarneySystem {

 public:
  KarneySystemWithq1AsIndependentVariable (KarneySystem ks)
      : KarneySystem(ks)
  { };

  void operator() (const State& qJ, State& dqJdt, const double t) const final
  {
    State parentSystemDerivatives;
    KarneySystem::operator()(qJ, parentSystemDerivatives, t);
    const double dq1dt = parentSystemDerivatives[0].real();

    dqJdt[0] = ActionAngle(1, parentSystemDerivatives[0].imag() / dq1dt);
    dqJdt[1] = parentSystemDerivatives[1] / dq1dt;
  }

};

using ErrorStepperType = runge_kutta_cash_karp54<State>;

class SurfaceCrossEventObserver : public PanosOde::Utils::PoincareSurfaceCrossEventObserver<State>
{
 protected:
  double distance_from_surface (const State& x) const override
  {
    return std::fmod(x[0].real() + M_PI, 2 * M_PI);
  }

  bool event_predicate (double distance) const override // distance decreases, only when we cross the surface
  {
    const bool crossed_surface = (distance - previous_distance_from_surface_) < 0;
    return crossed_surface;
  }

};


int main ()
{

  const double a = 2.1;
  const double nu = 30.23;
  const double Ham = 1200;
  KarneySystem ks1(a, nu);


  //integrate with adaptive step using make_controlled


  std::vector<double> init_r1{44.5, 44.8 ,45.1, 45.3,  45.5, 46.1, 46.5, 47.1,47.5,48.1,48.5,49.1,49.5};
  std::vector<State> init_states;
  for (const auto& r1:init_r1)
    {
      State x_in;
      const double myJ1 = r1 * r1 / 2;
      const double myJ2 = (Ham-myJ1)/nu;
      x_in[0] = ActionAngle(0, myJ1);
      x_in[1] = ActionAngle(0, myJ2);
      init_states.push_back(x_in);
    }
  double abs_err = 1.0e-12;
  double rel_err = 1.0e-10;
  auto controlled_stepper = make_controlled(abs_err, rel_err, ErrorStepperType());
  std::vector<double> t_out;
  std::vector<double> energy_out;
//    auto energy_observer = [&t_out,&energy_out,&ks1]
//        (const State & x, double t) {t_out.push_back(t);
//        energy_out.push_back(ks1.hamiltonian(x));};

  std::vector<double> w2_mod_2_pi;
  std::vector<double> r1;

  auto karney_observer = [&w2_mod_2_pi, &r1, &ks1]
      (const State& x, double /*t*/)
  {
      w2_mod_2_pi.push_back(fmod(x[1].real(), 2 * M_PI));
      r1.push_back(std::sqrt(2 * x[0].imag()));
  };


  auto karney_event_observer = PanosOde::Utils::makeEventFunctor(karney_observer, SurfaceCrossEventObserver(),
                                                ErrorStepperType(), KarneySystemWithq1AsIndependentVariable(ks1));


  Gnuplot g1;
  g1 << "plot '-' with dots lc rgb \"black\" \n";
  for (auto x_init:init_states)
    {
      std::cout << "adaptive integration output(1): " << ks1.hamiltonian(x_init) << std::endl;

      integrate_adaptive(controlled_stepper, ks1, x_init, 0.0, 1000.0, 0.01, karney_event_observer);
      std::cout << "adaptive integration output(2): " << ks1.hamiltonian(x_init) << std::endl;
    }
  g1.send1d(std::make_tuple(w2_mod_2_pi,r1));


  int keyboard;
  std::cin >> keyboard;
  return 0;
}