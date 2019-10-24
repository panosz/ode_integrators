//
// Created by Panagiotis Zestanakis on 2019-09-11.
//

#ifndef ODE_INTEGRATORS_ACTION_INTEGRATION_HPP
#define ODE_INTEGRATORS_ACTION_INTEGRATION_HPP

#include "action_integration_result.hpp"
#include "integration_utilities.hpp"
#include "hamiltonian_dynamic_system.hpp"

template<typename System>
ActionIntegrationResult
action_integration (System sys,
                    const typename DS::ExtendedSpaceState& first_point,
                    double max_time,
                    IntegrationOptions integrationOptions)
{
  using boost::math::double_constants::one_div_two_pi;

  const auto last_point = last_point_on_closed_orbit(sys, first_point, max_time, integrationOptions);
  const auto delta_point = last_point - first_point;

  SpecialIntegrals integrals{};

  integrals.beta = delta_point[static_cast<unsigned>(CoordinateTag::beta)];
  integrals.beta1 = delta_point[static_cast<unsigned>(CoordinateTag::beta1)];
  integrals.beta2 = delta_point[static_cast<unsigned>(CoordinateTag::beta2)];

  integrals.gamma = delta_point[static_cast<unsigned>(CoordinateTag::gamma)];
  integrals.gamma1 = delta_point[static_cast<unsigned>(CoordinateTag::gamma1)];
  integrals.gamma2 = delta_point[static_cast<unsigned>(CoordinateTag::gamma2)];

  const auto delta_J = delta_point[static_cast<unsigned>(CoordinateTag::J)];
  const auto theta_period = delta_point[static_cast<unsigned>(CoordinateTag::t)];
  const auto delta_phi = delta_point[static_cast<unsigned>(CoordinateTag::phi)];

  const auto normalized_delta_Action = delta_J * one_div_two_pi;

  return ActionIntegrationResult{normalized_delta_Action, theta_period, delta_phi, integrals};
}

template<typename System>
ActionIntegrationResult
action_integration (System sys,
                    const typename DS::PhaseSpaceState & first_point,
                    double max_time,
                    IntegrationOptions integrationOptions)
{
  return action_integration(sys, DS::phase_to_extended_space_state(first_point), max_time, integrationOptions);
}
#endif //ODE_INTEGRATORS_ACTION_INTEGRATION_HPP
