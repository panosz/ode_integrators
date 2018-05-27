//
// Created by Panagiotis Zestanakis on 25/05/18.
//
#include "henon_heiles_impl.hpp"
#include <cmath>
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/roots.hpp>
#include "myUtilities.hpp"
#include <boost/range/algorithm_ext/push_back.hpp>
#include <boost/numeric/odeint.hpp>

double Henon_Heiles::h0(const vector_type& q, const vector_type& p)
{
  const double half = boost::math::double_constants::half;
  const auto &[x, y] = q;
  const auto &[px, py] = p;

  return half * (std::pow(x, 2) + std::pow(y, 2)
                 + std::pow(px, 2) + std::pow(py, 2));
}

double Henon_Heiles::h1 (const vector_type& q)
{
  using boost::math::double_constants::third;
  using boost::math::double_constants::half;

  const auto &[x, y] = q;

  return std::pow(x, 2) * y - third * std::pow(y, 3);
}
vector_type Henon_Heiles::dh0dq (const vector_type& q)
{
  return q;
}
vector_type Henon_Heiles::dh0dp (const vector_type& p)
{
  return p;
}
vector_type Henon_Heiles::dh1dq (const vector_type& q)
{
  const auto &[x, y] = q;

  return vector_type{2 * x * y, std::pow(x, 2) - std::pow(y, 2)};
}
double Henon_Heiles::Hamiltonian (const vector_type& q, const vector_type& p)
{
  return h0(q, p) + h1(q);
}
vector_type Henon_Heiles::dHdp (const vector_type& q, const vector_type& p)
{
  return dh0dp(p);  //dh1dp==0
}
void Henon_Heiles::dpdt (const vector_type& q, vector_type& p_time_derivative)
{
  const auto v0 = dh0dq(q);
  const auto v1 = dh1dq(q);
  p_time_derivative[0] = -v0[0] - v1[0];
  p_time_derivative[1] = -v0[1] - v1[1];
}
void Henon_Heiles::dqdt (const vector_type& p, vector_type& q_time_derivative)
{
  q_time_derivative = p;
}
double kinetic_energy_at_y (double total_energy, double y)
//x==0
{
  using namespace boost::math::double_constants;
  return total_energy - 1. / 2 * std::pow(y, 2) + 1. / 3. * std::pow(y, 3);
}
double derivative_of_kinetic_energy_at_y (double y)
{
  return std::pow(y, 2) - y;
}
std::pair<double, double> kinetic_energy_and_derivative (double total_energy, double y)
{
  return std::make_pair(
      kinetic_energy_at_y(total_energy, y),
      derivative_of_kinetic_energy_at_y(y)
  );
}
double kinetic_energy_root (double total_energy, double yguess)
{
  double ymin = -1;
  double ymax = 10;
  size_t max_iter = 100;
  const int digits = 14;

  auto functor = [total_energy] (double y)
  { return kinetic_energy_and_derivative(total_energy, y); };

  auto root = boost::math::tools::
  newton_raphson_iterate(functor, yguess, ymin, ymax, digits, max_iter);
  return root;
}
vector_type get_q (const State& state)
{
  return vector_type{state[0],state[1]};
}
vector_type get_p (const State& state)
{
  return vector_type{state[2],state[3]};
}
std::vector<State> make_initial_y_px_at_E (int numberOfPositions, double energy_level)
{

  std::vector<State> initial_positions;
  using namespace boost::math::double_constants;

  auto ymin = kinetic_energy_root(energy_level, -1);
  auto ymax = energy_level > sixth ? 1 : kinetic_energy_root(energy_level, 0.5);

  auto distance = ymax - ymin;

  ymin += 0.01 * distance;
  ymax -= 0.01 * distance;



  for (const auto& y: PanosUtilities::linspace(ymin, ymax, numberOfPositions))
    {
      const double px = std::sqrt(2 * kinetic_energy_at_y(energy_level, y));
      initial_positions.push_back(State{0.0, y, px, 0.0});
    }
  return initial_positions;
}
std::ostream& operator<< (std::ostream& os, const vector_type& v)
{
  return os << v[0] << ' ' << v[1];
}
std::ostream& operator<< (std::ostream& os, const State& s)
{
  return os << get_q(s) << " , " << get_p(s) << '\n';
}



std::pair<std::vector<double>, std::vector<double> >
henon_heiles_poincare_surface (double total_energy, double time_integration, size_t noOfInitialPoints)
{
  auto initial_states = make_initial_y_px_at_E(noOfInitialPoints, total_energy);

  for (const auto&s : initial_states)
    {
      const auto q=get_q(s);
      const auto p=get_p(s);
      std::cout << Henon_Heiles::Hamiltonian(q, p) << '\n';
    }

  using ErrorStepperType = boost::numeric::odeint::runge_kutta_cash_karp54<State>;

  double abs_err = 1.0e-16;
  double rel_err = 1.0e-13;
  auto controlled_stepper = boost::numeric::odeint::make_controlled(abs_err, rel_err, ErrorStepperType());


  std::vector<double> y;
  std::vector<double> py;

  auto push_back_observer = [&y, &py]
      (const State& s, double /*t*/)
  {
      y.push_back(s[1]);
      py.push_back(s[3]);
  };


  auto henon_heiles_event_observer = PanosOde::Utils::makeEventFunctor(push_back_observer, SurfaceCrossEventObserver(),
                                                                       ErrorStepperType(), CrossingDerivatives());



  for (auto state :initial_states)
    {
      std::cout << "statring point " << state;

      boost::numeric::odeint::integrate_adaptive
          (controlled_stepper, TrajectoryDerivatives(), state,
           0.0, time_integration, 0.0001, henon_heiles_event_observer);
      std::cout << "point after step " << state;
      const auto q=get_q(state);
      const auto p=get_p(state);
      std::cout << "energy after step= " << Henon_Heiles::Hamiltonian(q, p) << "\n\n";
    }

  return std::make_pair(y,py);

}
void TrajectoryDerivatives::operator() (const State& s, State& dsdt, const double)
{
  const auto q= get_q(s);
  const auto p= get_p(s);
  auto &[dxdt,dydt,dpxdt,dpydt] = dsdt;
  vector_type q_deriv;
  vector_type p_deriv;
  dqdt(p, q_deriv);
  dpdt(q, p_deriv);

  dxdt= q_deriv[0];
  dydt = q_deriv[1];
  dpxdt= p_deriv[0];
  dpydt=p_deriv[1];


}
void CrossingDerivatives::operator() (const State& s, State& dsdt, const double t)
{

  TrajectoryDerivatives::operator()(s,dsdt,t);
  dsdt[0]=1;
  dsdt[1] /=dsdt[0];
  dsdt[2] /=dsdt[0];
  dsdt[3] /=dsdt[0];

}
double SurfaceCrossEventObserver::distance_from_surface (const State& state) const
{
  const auto q = get_q(state);
  const auto &[x, y] = q;
  return x;
}
bool SurfaceCrossEventObserver::event_predicate (double distance) const
{
  const bool crossed_surface = (distance > 0 && previous_distance_from_surface_ < 0);
  return crossed_surface;
}
