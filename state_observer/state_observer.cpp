//
// Created by Panagiotis Zestanakis on 18/09/18.
//

#include <iostream>
#include <fstream>
#include <string_view>
#include <vector>
#include <stdexcept>

#include <boost/numeric/odeint.hpp>
#include "../utilities.hpp"
#include <boost/range/iterator_range.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm_ext/push_back.hpp>
#include <myUtilities/zero_crossing.hpp>
#include <myUtilities/wrap.hpp>

using namespace boost::numeric::odeint;

/* The type of container used to hold the state vector */
typedef std::array<double, 3> State;

std::ostream& operator<< (std::ostream& os, const State& s)
{
  for (const auto& item :s)
    os << item << ' ';
  return os;
}

class System {
 public:
  System () = default;

  void operator() (const State& s, State& dsdt, const double /*t*/)
  {
    dsdt[0] = -sin(s[1]);
    dsdt[1] = s[0];
    dsdt[2] = -cos(s[1]);

  }
};

class SystemNormalToSurface {
 private:
  System sys_;
  unsigned index;

 public:
  SystemNormalToSurface (System sys, unsigned surface_index)
      : sys_{sys}, index{surface_index}
  { }

  void operator() (const State& s, State& dsdt, const double t)
  {
    sys_(s, dsdt, t);

    const auto normalization = dsdt[index];

    for (auto& item : dsdt)
      item = item / normalization;

  }

};

double event_surface (const State& s)
{
  return s[0];
}

using ErrorStepperType = runge_kutta_cash_karp54<State>;
using ControlledStepperType = controlled_runge_kutta<ErrorStepperType>;

template<typename Functor>
State step_on_surface (SystemNormalToSurface systemNormalToSurface, Functor surface, State s, ErrorStepperType stepper)
{
  State state_on_surface{};
  double distance_from_surface = -surface(s);
  stepper.do_step(systemNormalToSurface, s, 0, state_on_surface, distance_from_surface);
  return state_on_surface;
}

struct OrbitCrossOutput {
    State initial_point{};
    std::vector<State> cross_points{};
};

std::ostream& operator<< (std::ostream& os, const OrbitCrossOutput& orbitCrossOutput)
{
  os << "# initial point:\n";
  os << orbitCrossOutput.initial_point << '\n';
  os << "# cross points:\n";
  for (const auto& item : orbitCrossOutput.cross_points)
    os << item << '\n';
  return os;
}

struct IntegrationOptions {
    double abs_err;
    double rel_err;
    double dt;

    IntegrationOptions (double abs_error, double rel_error, double init_time_step)
        : abs_err{abs_error}, rel_err{rel_error}, dt{init_time_step}
    { }
};

template<typename Functor>
OrbitCrossOutput pick_orbit_points_that_cross_surface (System sys,
                                                       Functor surface,
                                                       State init_state,
                                                       double integration_time,
                                                       int direction,
                                                       IntegrationOptions options)
{
  OrbitCrossOutput output{};
  output.initial_point = init_state;

  const double integration_start_time = 0;

  const auto controlled_stepper = make_controlled(options.abs_err, options.rel_err, ErrorStepperType());

  auto orbit_iterators = make_adaptive_range(controlled_stepper,
                                             sys, init_state, integration_start_time, integration_time, options.dt);
  auto orbit_points = boost::make_iterator_range(orbit_iterators.first, orbit_iterators.second);

  PanosUtilities::zero_cross_transformed(orbit_points, std::back_inserter(output.cross_points),
                                         surface, direction);

  return output;
}

template<typename Functor>
OrbitCrossOutput trace_cross_points_on_cross_surface (SystemNormalToSurface systemNormalToSurface,
                                                      Functor surface,
                                                      const OrbitCrossOutput& approximate_cross_output
                                                     )
{
  OrbitCrossOutput output{};
  output.initial_point = approximate_cross_output.initial_point;

  const auto my_projection = [systemNormalToSurface, surface] (State s)
  { return step_on_surface(systemNormalToSurface, surface, s, ErrorStepperType()); };

  boost::push_back(output.cross_points,
                   approximate_cross_output.cross_points | boost::adaptors::transformed(my_projection));

  return output;
}



void write_to_file(const char* filename, const OrbitCrossOutput & orbitCrossOutput)
{
  std::ofstream file(filename);
  if (!file)
    throw std::runtime_error("cannot open output file");

  file<<orbitCrossOutput;

  file.close();
}

int main ()
{

  const auto options = IntegrationOptions(1e-12, 1e-12, 1e-5);

  const int zero_cross_positive_direction = 1;
  const double integration_time = 1000;

  State init_state{1.1, 0, 0};
  auto my_sys = System();
  auto my_cross_sys = SystemNormalToSurface(my_sys, 0);


  const auto approximate_cross_output = pick_orbit_points_that_cross_surface(my_sys, event_surface,
                                                                             init_state, integration_time,
                                                                             zero_cross_positive_direction,
                                                                             options);

  const auto exact_on_surface_output = trace_cross_points_on_cross_surface(my_cross_sys,
                                                                           event_surface,
                                                                           approximate_cross_output
                                                                          );


  std::cout << "we have as many as " << approximate_cross_output.cross_points.size() << " zero cross points\n";

  std::cout << "approximate cross output:\n";
  std::cout << approximate_cross_output << '\n';

  write_to_file("approximate_cross_output.txt", approximate_cross_output);
  std::cout << "exact cross output:\n";
  std::cout << exact_on_surface_output << '\n';

  write_to_file("exact_cross_output.txt", exact_on_surface_output);


}
