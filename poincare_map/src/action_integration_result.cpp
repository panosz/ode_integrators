//
// Created by Panagiotis Zestanakis on 10/09/19.
//

#include "action_integration_result.hpp"
double ActionIntegrationResult::d2K_dF2 () const
{
  using boost::math::pow;
  const auto G = g_factor();
  const auto T = theta_period();
  const auto beta2 = integrals_.beta2;
  const auto gamma2 = integrals_.gamma2;

  return G * d2K_dJdF() + gamma2 / T - G / pow<2>(T) * beta2;
}
double ActionIntegrationResult::d2K_dJdF () const
{
  using boost::math::pow;
  using boost::math::double_constants::one_div_two_pi;

  const auto beta1 = integrals_.beta1;
  const auto gamma1 = integrals_.gamma1;
  const auto G = g_factor();
  const auto T = theta_period();
  const auto omega_sq = pow<2>(omega_theta());

  return omega_sq * one_div_two_pi * (gamma1 - G / T * beta1);
}
double ActionIntegrationResult::dT_dJ () const
{
  //\[Omega]/(2 \[Pi]) beta1Integral;

  return integrals_.beta1 / theta_period();
}
double ActionIntegrationResult::dT_dF () const
{
  using boost::math::double_constants::one_div_two_pi;
  //dTdF = 1/(2 \[Pi]) beta2Integral + myG*dTdJ;
  const auto G = g_factor();
  const auto beta2 = integrals_.beta2;

  return one_div_two_pi * beta2 + G * dT_dJ();
}
ActionIntegrationResult::ActionIntegrationResult (double Action, double theta_period, double delta_phi, SpecialIntegrals integrals)
    : Action_{Action}, theta_period_{theta_period}, delta_phi_{delta_phi}, integrals_{integrals}
{ }
double ActionIntegrationResult::Action () const noexcept
{
  return Action_;
}
double ActionIntegrationResult::theta_period () const noexcept
{
  return theta_period_;
}
double ActionIntegrationResult::delta_phi () const noexcept
{
  return delta_phi_;
}
double ActionIntegrationResult::omega_theta () const
{
  using boost::math::double_constants::two_pi;
  return two_pi / theta_period();
}
double ActionIntegrationResult::g_factor () const noexcept
{
  using boost::math::double_constants::one_div_two_pi;
  return delta_phi() * one_div_two_pi;
}
double ActionIntegrationResult::omega_phi () const
{
  return delta_phi() / theta_period();
}
double ActionIntegrationResult::domega_dJ () const
{
  using boost::math::double_constants::two_pi;
  using boost::math::pow;
  //d\[Omega]dJ = -((2 \[Pi])/T^2) dTdJ;
  return -two_pi * dT_dJ() / pow<2>(theta_period());
}
double ActionIntegrationResult::d2K_dJ2 () const
{
  using boost::math::double_constants::two_pi;
  using boost::math::pow;

  return -pow<3>(omega_theta()) * integrals_.beta1 / pow<2>(two_pi);
}
double ActionIntegrationResult::one_div_two_pi_gamma () const
{
  using boost::math::double_constants::one_div_two_pi;

  return one_div_two_pi * integrals_.gamma;
}
double ActionIntegrationResult::domega_dF () const
{
  using boost::math::double_constants::two_pi;
  using boost::math::pow;
  //d\[Omega]dF = -((2 \[Pi])/T^2) dTdF;

  return -two_pi / pow<2>(theta_period()) * dT_dF();
}
