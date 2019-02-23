//
// Created by Panagiotis Zestanakis on 18/09/18.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <boost/numeric/odeint.hpp>
#include "../utilities.hpp"
#include <boost/range/iterator_range.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm_ext/push_back.hpp>
#include <myUtilities/zero_crossing.hpp>

using namespace boost::numeric::odeint;

/* The type of container used to hold the state vector */
typedef std::vector<double> State;

class System {
 public:
  System () = default;

  void operator() (const State& x, State& dxdt, const double /*t*/)
  {
    dxdt[0] = -sin(x[1]);
    dxdt[1] = x[0];

  }
};

class SystemCrossSurface {
 private:
  System sys_;
  int index;

 public:
  SystemCrossSurface (System sys, int surface_index)
      : sys_{sys}, index{surface_index}
  { }

  void operator() (const State& x, State& dxdt, const double t)
  {
    sys_(x, dxdt, t);

    const auto normalization = dxdt[index];

    for (auto& p : dxdt)
      p = p / normalization;

  }

};

double event_surface (const State& s)
{
  return s[0];
}

using ErrorStepperType = runge_kutta_cash_karp54<State>;
using ControlledStepperType = controlled_runge_kutta<ErrorStepperType>;

State step_on_surface (SystemCrossSurface scf, State s, ErrorStepperType stepper)
{
  State state_on_surface{0,0};
  double distance_from_surface = -event_surface(s);
  stepper.do_step(scf,s,0,state_on_surface,distance_from_surface);
  std::cout << "step_on_surface. from "<<s[0]<<','<<s[1]<<'\n';
  return state_on_surface;
}

int main ()
{

  const double abs_err = 1.0e-16;
  const double rel_err = 1.0e-16;
  const double max_dt = 0.001;
  const auto controlled_stepper = make_controlled(abs_err, rel_err, max_dt, ErrorStepperType());

  State init_state{1.1, 0};
  auto my_sys = System();
  auto my_cross_sys = SystemCrossSurface(my_sys,0);

  auto orbit_iterators = make_adaptive_range(controlled_stepper, System(), init_state, 0.0, 10.0, 1e-5);
  auto orbit_points = boost::make_iterator_range(orbit_iterators.first, orbit_iterators.second);

  std::vector<State> zero_cross_points{};
  PanosUtilities::zero_cross_transformed(orbit_points, std::back_inserter(zero_cross_points), event_surface);
  std::cout << "we have as many as " << zero_cross_points.size() << " zero cross points\n";

//  boost::push_back(zero_cross_points,orbit_points);

//  std::ofstream outfile("dummy.txt");
//
//  for (auto out:orbit_points)
//    outfile << out.second <<' '<< out.first[0]<<' '<< out.first[1] << '\n';
////
//  outfile.close();


  std::cout << "zero_cross:\n";
  for (auto out:zero_cross_points)
    std::cout << out[0] << ' ' << out[1] << '\n';

  std::vector<State> on_surface_points{};

  const auto my_projection = [my_cross_sys](State s){return step_on_surface(my_cross_sys,s,ErrorStepperType());};

  boost::push_back(on_surface_points, zero_cross_points| boost::adaptors::transformed(my_projection));

  std::cout << "zero_cross projection:\n";
  for (auto out:on_surface_points)
    std::cout << out[0] << ' ' << out[1] << '\n';

}
