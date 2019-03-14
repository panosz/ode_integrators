//
// Created by Panagiotis Zestanakis on 26/02/19.
//

#ifndef ODE_INTEGRATORS_INTEGRATION_UTILITIES_HPP
#define ODE_INTEGRATORS_INTEGRATION_UTILITIES_HPP

#include <boost/numeric/odeint.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm_ext/push_back.hpp>
#include <myUtilities/zero_crossing.hpp>
#include <myUtilities/wrap.hpp>

#include "system_and_poincare_surface.hpp"
#include "samplingCollections.hpp"

using namespace boost::numeric::odeint;

template<typename System>
using ErrorStepperType = runge_kutta_cash_karp54<typename System::StateType, double, typename System::StateType, double, vector_space_algebra>;

template<typename System>
using ControlledStepperType = controlled_runge_kutta<ErrorStepperType<System> >;

template<typename System>
typename System::StateType step_on_surface (SystemAndPoincareSurface<System> sys,
                                            typename System::StateType s,
                                            ErrorStepperType<System> stepper)
{
  typename System::StateType state_on_surface{};

  const auto systemPerpendicularToSurface = [sys] (const typename System::StateType& state,
                                                   typename System::StateType& dsdt, const double t)
  {
      sys.eval_perp_to_surface(state, dsdt, t);
  };

  double distance_from_surface = -sys.surface_eval(s);
  stepper.do_step(systemPerpendicularToSurface, s, 0, state_on_surface, distance_from_surface);
  return state_on_surface;
}

struct IntegrationOptions {
    double abs_err;
    double rel_err;
    double dt;

    IntegrationOptions (double abs_error, double rel_error, double init_time_step)
        : abs_err{abs_error}, rel_err{rel_error}, dt{init_time_step}
    { }
};

template<typename Range, typename State>
class ParticleOrbit
{
 public:
  using StateType =State;
 private:
  StateType init_state_;
  Range range_;
  bool is_consumed_;
 public:
  ParticleOrbit (StateType init_state, Range range)
      : init_state_(init_state), range_(range), is_consumed_{false}
  {};
  StateType init_state () const noexcept
  {
    return init_state_;
  }
  bool is_consumed () const noexcept
  {
    return is_consumed_;
  }

  Range range ()
  /// can be consumed only once
  {
    if (is_consumed_)
      throw std::runtime_error("Orbit may be already consumed");

    is_consumed_=true;
    return range_;
  }
};

template<typename System>
auto make_ParticleOrbit(System sys,
                      typename System::StateType init_state,
                      double integration_time,
                      IntegrationOptions options)
/// \brief returnded orbit has range of type addaptive_range
/// \tparam System
/// \param sys
/// \param init_state
/// \param integration_time
/// \param options
/// \return
{
  const double integration_start_time = 0;

  const auto controlled_stepper = make_controlled(options.abs_err, options.rel_err, ErrorStepperType<System>());



  auto orbit_iterators = make_adaptive_range(controlled_stepper,
                                             sys,
                                             init_state, integration_start_time, integration_time, options.dt);

  auto orbit_range =  boost::make_iterator_range(orbit_iterators.first, orbit_iterators.second);

  return ParticleOrbit(init_state,orbit_range);
}


template<typename System>
auto make_TimeParticleOrbit(System sys,
                        typename System::StateType init_state,
                        double integration_time,
                        IntegrationOptions options)
///returned orbit has range of type addaptive_time_range

{
  const double integration_start_time = 0;

  const auto controlled_stepper = make_controlled(options.abs_err, options.rel_err, ErrorStepperType<System>());



  auto orbit_iterators = make_adaptive_time_range(controlled_stepper,
                                             sys,
                                             init_state, integration_start_time, integration_time, options.dt);

  auto orbit_range =  boost::make_iterator_range(orbit_iterators.first, orbit_iterators.second);

  return ParticleOrbit(init_state,orbit_range);
}





template<typename OrbitType>
OrbitCrossOutput<typename OrbitType::StateType>
pick_orbit_points_that_cross_surface (OrbitType orbit, Surface surface )
{
  const double MAX_SURFACE_CROSS_DISTANCE = boost::math::double_constants::half_pi;

  OrbitCrossOutput<typename OrbitType::StateType> output{};
  output.initial_point = orbit.init_state();



  auto surface_fun = [surf=surface] (const typename OrbitType::StateType& s)
  { return surf.eval(s); };

  std::cout << "start following orbit" << std::endl;


  PanosUtilities::zero_cross_transformed(orbit.range(), std::back_inserter(output.cross_points),
                                         surface_fun, MAX_SURFACE_CROSS_DISTANCE, surface.direction);

  std::cout << "calculated " << output.cross_points.size() << " cross points" << std::endl;
  return output;
}

template<typename System>
OrbitCrossOutput<typename System::StateType>
trace_cross_points_on_cross_surface (SystemAndPoincareSurface<System> sys,
                                     const OrbitCrossOutput<typename System::StateType>& approximate_cross_output)
{
  OrbitCrossOutput<typename System::StateType> output{};
  output.initial_point = approximate_cross_output.initial_point;

  const auto my_projection = [sys] (typename System::StateType s)
  { return step_on_surface(sys, s, ErrorStepperType<System>()); };

  boost::push_back(output.cross_points,
                   approximate_cross_output.cross_points | boost::adaptors::transformed(my_projection));

  return output;
}

template<typename System>
OrbitCrossOutput<typename System::StateType>
trace_on_poincare_surface (SystemAndPoincareSurface<System> sys_and_pc,
                           typename System::StateType init_state,
                           double integration_time,
                           IntegrationOptions options)
{
  OrbitCrossOutput<typename System::StateType> output{};

  auto orbit = make_ParticleOrbit(sys_and_pc,init_state,integration_time,options);

  const auto approximate_points = pick_orbit_points_that_cross_surface(orbit,sys_and_pc.poincare_surface());

  return trace_cross_points_on_cross_surface(sys_and_pc, approximate_points);

}

template<typename System>
std::vector<OrbitCrossOutput<typename System::StateType>>
trace_on_poincare_surface (SystemAndPoincareSurface<System> sys_and_pc,
                           std::vector<typename System::StateType> init_states,
                           double integration_time,
                           IntegrationOptions options)
{
  std::vector<OrbitCrossOutput<typename System::StateType>> output;

  for (const auto& state: init_states)
    output.push_back(trace_on_poincare_surface(sys_and_pc, state, integration_time, options));

  return output;

}

#endif //ODE_INTEGRATORS_INTEGRATION_UTILITIES_HPP
