//
// Created by Panagiotis Zestanakis on 10/09/19.
//

#ifndef ODE_INTEGRATORS_ACTIONINTEGRATIONRESULT_HPP
#define ODE_INTEGRATORS_ACTIONINTEGRATIONRESULT_HPP

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/pow.hpp>


struct SpecialIntegrals {
    double beta = 0;
    double gamma = 0;
    double beta1 = 0;
    double beta2 = 0;
    double gamma1 = 0;
    double gamma2 = 0;
};


class ActionIntegrationResult {
 private:
  double Action_;
  double theta_period_;
  double delta_phi_;
  SpecialIntegrals integrals_;

  double dT_dJ () const;

  double dT_dF () const;

 public:
  ActionIntegrationResult (double Action, double theta_period, double delta_phi, SpecialIntegrals integrals);
  double Action () const noexcept;
  double theta_period () const noexcept;
  double delta_phi () const noexcept;
  double omega_theta () const;
  double g_factor () const noexcept;
  double omega_phi () const;
  double domega_dJ () const;
  double d2K_dJ2 () const;
  double one_div_two_pi_gamma () const;
  double domega_dF () const;
  double d2K_dJdF () const;
  double d2K_dF2 () const;
};
#endif //ODE_INTEGRATORS_ACTIONINTEGRATIONRESULT_HPP
