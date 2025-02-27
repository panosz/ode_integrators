//
// Created by Panagiotis Zestanakis on 26/02/19.
//

#ifndef ODE_INTEGRATORS_INTEGRATION_UTILITIES_HPP
#define ODE_INTEGRATORS_INTEGRATION_UTILITIES_HPP
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/iterator/times_iterator.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/range/empty.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm_ext/push_back.hpp>
#include <myUtilities/zero_crossing.hpp>

#include "system_and_poincare_surface.hpp"
#include "samplingCollections.hpp"

#include "integration_options.hpp"
#include "orbit_closing_check.hpp"

using namespace boost::numeric::odeint;

template<typename System>
using ErrorStepperType = runge_kutta_cash_karp54<typename System::StateType,
                                                 double,
                                                 typename System::StateType,
                                                 double,
                                                 vector_space_algebra>;

template<typename System>
using ControlledStepperType = controlled_runge_kutta<ErrorStepperType<System> >;



template<typename System>
typename System::StateType step_on_surface(SystemAndPoincareSurface<System> sys,
                                           typename System::StateType s,
                                           ErrorStepperType<System> stepper)
{
  typename System::StateType state_on_surface{};

  const auto systemPerpendicularToSurface =
            [sys](const typename System::StateType& state,
                  typename System::StateType& dsdt,
                  const double t)
  {
      sys.eval_perp_to_surface(state, dsdt, t);
  };

  const double distance_from_surface = -sys.surface_eval(s);
  stepper.do_step(systemPerpendicularToSurface,
                  s,
                  0,
                  state_on_surface,
                  distance_from_surface);

  return state_on_surface;
}


template<typename System>
auto make_OrbitRange (System sys,
                      typename System::StateType& init_state,
                      double integration_time,
                      IntegrationOptions options)
/// \brief returned orbit has range of type adaptive_range
/// \tparam System
/// \param sys
/// \param init_state
/// \param integration_time
/// \param options
/// \return
{
  const double integration_start_time = 0;

  const auto controlled_stepper = make_controlled(options.abs_err,
                                                  options.rel_err,
                                                  ErrorStepperType<System>());

  auto orbit_iterators = make_adaptive_range(controlled_stepper,
                                             sys,
                                             init_state,
                                             integration_start_time,
                                             integration_time, options.dt);

  return boost::make_iterator_range(orbit_iterators.first,
                                    orbit_iterators.second);
}

template<typename System>
auto make_OrbitTimeRange (System sys,
                          typename System::StateType& init_state,
                          double integration_time,
                          IntegrationOptions options)
///returned orbit has range of type adaptive_time_range

{
  const double integration_start_time = 0;

  const auto controlled_stepper = make_controlled(options.abs_err,
                                                  options.rel_err,
                                                  ErrorStepperType<System>());

  auto orbit_iterators = make_adaptive_time_range(controlled_stepper,
                                                  sys,
                                                  init_state,
                                                  integration_start_time,
                                                  integration_time,
                                                  options.dt);

  return boost::make_iterator_range(orbit_iterators.first,
                                    orbit_iterators.second);
}

template<typename System, typename TimeIterator>
auto make_OrbitTimesRange (System sys,
                           typename System::StateType& init_state,
                           TimeIterator t_start ,
                           TimeIterator t_end ,
                           IntegrationOptions options)
///returned orbit has range of type times_range

{

  const auto controlled_stepper = make_controlled(options.abs_err,
                                                  options.rel_err,
                                                  ErrorStepperType<System>());



  auto orbit_iterators = boost::numeric::odeint::make_times_range(controlled_stepper,
                                          sys,
                                          init_state,
                                          t_start,
                                          t_end,
                                          options.dt);

  return boost::make_iterator_range(orbit_iterators.first,
                                    orbit_iterators.second);
}

template<typename System>
class ParticleOrbit {
 public:
  using StateType = typename System::StateType;
 private:
  System sys_;
  StateType init_state_;
  double integration_time_;
  IntegrationOptions options_;
  mutable StateType state_for_range_;
  mutable StateType state_for_time_range_;
  mutable StateType state_for_times_range_;
 public:

  ParticleOrbit (System sys,
                 StateType init_state,
                 double integration_time,
                 IntegrationOptions options)
               : sys_{sys},
                 init_state_{init_state},
                 integration_time_{integration_time},
                 options_{options},
                 state_for_range_{init_state},
                 state_for_time_range_{init_state},
                 state_for_times_range_{init_state}
  { }

  StateType init_state () const
  {
    return init_state_;
  }
  System sys () const
  {
    return sys_;
  }

  auto range ()
  /// should not follow more than two instances of the same orbit range
  /// simultaneously
  {
    state_for_range_ = init_state_;

    return make_OrbitRange(sys_,
                           state_for_range_,
                           integration_time_,
                           options_);
  }

  auto time_range ()
  /// should not follow more than two instances of the same orbit range
  /// simultaneously
  {
    state_for_time_range_ = init_state_;
    return make_OrbitTimeRange(sys_, state_for_time_range_, integration_time_, options_);
  }

template<typename TimeIterator>
  auto times_range(TimeIterator t_start,
                   TimeIterator t_end)
  {
    state_for_times_range_ = init_state_;
    return make_OrbitTimesRange(sys_, state_for_times_range_, t_start, t_end, options_);
  }


};


template<typename System>
auto make_ParticleOrbit (System sys,
                         typename System::StateType init_state,
                         double integration_time,
                         IntegrationOptions options)
{
  return ParticleOrbit<System>(sys, init_state, integration_time, options);
}


template<typename OrbitType>
std::vector<typename OrbitType::StateType>
rough_cross_points (OrbitType& orbit, Surface surface)
{
  const double MAX_SURFACE_CROSS_DISTANCE = boost::math::double_constants::half_pi;

  std::vector<typename OrbitType::StateType> output{};

  auto surface_fun = [surf = surface] (const typename OrbitType::StateType& s)
  { return surf.eval(s); };

  std::cout << "start following orbit" << std::endl;

  PanosUtilities::zero_cross_transformed(orbit.range(),
                                         std::back_inserter(output),
                                         surface_fun,
                                         MAX_SURFACE_CROSS_DISTANCE,
                                         surface.direction);

  std::cout << "calculated " << output.size() << " cross points" << std::endl;
  return output;
}


template<typename System>
typename System::StateType
accurate_from_rough_cross_point (SystemAndPoincareSurface<System> sys,
                                 const typename System::StateType&
                                                       rough_cross_point)
{
  return step_on_surface(sys, rough_cross_point, ErrorStepperType<System>());
}


template<typename System>
std::vector<typename System::StateType>
accurate_from_rough_cross_points (SystemAndPoincareSurface<System> sys,
                                  const std::vector<typename System::StateType>&
                                                        rough_cross_input)
{
  std::vector<typename System::StateType> fine_cross_output{};

  const auto refine_point = [sys] (typename System::StateType s)
  { return accurate_from_rough_cross_point(sys, s); };

  boost::push_back(fine_cross_output,
                   rough_cross_input | boost::adaptors::transformed(refine_point));

  return fine_cross_output;
}

template<typename System>
OrbitCrossOutput
trace_on_poincare_surface (SystemAndPoincareSurface<System> sys_and_pc,
                           typename System::StateType init_state,
                           double integration_time,
                           IntegrationOptions options)
{
  auto orbit = make_ParticleOrbit(sys_and_pc,
                                  init_state,
                                  integration_time,
                                  options);

  const auto approximate_points =
                rough_cross_points(orbit,
                                   sys_and_pc.poincare_surface());

  auto cross_points = accurate_from_rough_cross_points(sys_and_pc,
                                                       approximate_points);

  return make_OrbitCrossOutput(init_state, cross_points);
}

template<typename System>
std::vector<OrbitCrossOutput>
trace_on_poincare_surface (SystemAndPoincareSurface<System> sys_and_pc,
                           std::vector<typename System::StateType> init_states,
                           double integration_time,
                           IntegrationOptions options)
{
  std::vector<OrbitCrossOutput> output;

  for (const auto& state: init_states)
  {
    output.push_back(trace_on_poincare_surface(sys_and_pc,
                                               state,
                                               integration_time,
                                               options));
  }

  return output;

}

namespace
{

    template<typename T>
    class BlackHoleContainer {
     public:
      using value_type = T;
      BlackHoleContainer () = default;
      template<typename S>
      void push_back (const S&) const noexcept
      { }
      void pop_back () const noexcept
      { }
    };


    template<typename System>
    auto system_and_poincare_surface_for_closed_orbit_integration
                            (System sys,
                             const typename System::StateType& init_state)
    {
      const auto zero_cross_position =
         init_state[static_cast<unsigned >(CoordinateTag::q)];

      typename System::StateType init_derivatives{};

      sys(init_state, init_derivatives, 0);

      const auto init_dqdt =
        init_derivatives[static_cast<unsigned >(CoordinateTag::q)];

      //direction is the sign of dqdt
      const int zero_cross_direction = (init_dqdt > 0) - (init_dqdt < 0);


      auto my_poincare_surface = Surface{CoordinateTag::q,
                                         zero_cross_position,
                                         zero_cross_direction};

      return make_system_and_poincare_surface(sys, my_poincare_surface);
    }

        template<typename System, typename OutputContainer>
    auto
    integrate_along_closed_orbit_impl (System sys,
                                       const typename System::StateType&
                                                       first_point,
                                       OutputContainer& outputContainer,
                                       double max_time,
                                       IntegrationOptions integrationOptions)
    {
      auto my_system_and_pc =
        system_and_poincare_surface_for_closed_orbit_integration(sys,
                                                                 first_point);

      auto orbit1 = make_ParticleOrbit(my_system_and_pc,
                                       first_point,
                                       max_time,
                                       integrationOptions);

      auto orbit_range = orbit1.range();

      auto close_enough_to_initial_point =
        OrbitClosing::make_P_Checker(first_point, integrationOptions);

      const auto begin_range = orbit_range.begin();
      const auto end_range = orbit_range.end();
      auto range_it = begin_range;

      while (true)
        {
          const auto poinc_surf_direction =
            my_system_and_pc.poincare_surface().direction;

          range_it =
            PanosUtilities::copy_until_zero_cross_transformed(
                                           range_it,
                                           end_range,
                                           std::back_inserter(outputContainer),
                                           my_system_and_pc.poincare_surface(),
                                           1.0,
                                           poinc_surf_direction);

          if (range_it == orbit_range.end())
            {

              throw std::runtime_error("orbit did not close yet");
            }

          else
            {
              const auto refined_point =
                accurate_from_rough_cross_point(my_system_and_pc, *range_it);

              if (close_enough_to_initial_point(refined_point))
                {
                  outputContainer.pop_back();
                  outputContainer.push_back(refined_point);
                  return refined_point;
                }
            }

        }

    }
}

template<typename System, typename Range>
std::vector<typename System::StateType>
orbit_points_at_times(System sys,
                      const typename System::StateType& first_point,
                      const Range& times,
                      IntegrationOptions integrationOptions)
{
  std::vector<typename System::StateType> output{};
  if (boost::empty(times))
    return output;

  const double max_time = 0; // a dummy value
  auto orbit = make_ParticleOrbit(sys,
                                  first_point,
                                  max_time,
                                  integrationOptions);

  const auto orbit_points = orbit.times_range(std::cbegin(times),std::cend(times));



  boost::range::push_back(output,orbit_points);

  return output;

}

template<typename System>
std::vector<typename System::StateType>
integrate_along_closed_orbit (System sys,
                              const typename System::StateType& first_point,
                              double max_time,
                              IntegrationOptions integrationOptions)
{

  std::vector<typename System::StateType> output{};
  integrate_along_closed_orbit_impl(sys,
                                    first_point,
                                    output,
                                    max_time,
                                    integrationOptions);

  return output;
}


template<typename System>
typename System::StateType
last_point_on_closed_orbit (System sys,
                            const typename System::StateType& first_point,
                            double max_time,
                            IntegrationOptions integrationOptions)
{

  BlackHoleContainer<typename System::StateType> no_output{};
  return integrate_along_closed_orbit_impl(sys,
                                           first_point,
                                           no_output,
                                           max_time,
                                           integrationOptions);
}

#endif //ODE_INTEGRATORS_INTEGRATION_UTILITIES_HPP
