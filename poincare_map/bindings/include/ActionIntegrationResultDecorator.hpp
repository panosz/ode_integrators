#ifndef ACTIONINTEGRATIONRESULTDECORATOR_HPP_ETQNS6TA
#define ACTIONINTEGRATIONRESULTDECORATOR_HPP_ETQNS6TA
#include "action_integration_result.hpp"
#include "state_bindings.hpp"

namespace ActionIntegrationBindings{
 class ActionIntegrationResultDecorator
  {
    private:
      ActionIntegrationResult air_;
    public:
      ActionIntegrationResultDecorator(const ActionIntegrationResult& air);
      ActionIntegrationResultDecorator(double x, double y, double z, SpecialIntegrals si);
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
      boost::python::numpy::ndarray hessian() const;


  };

  void export_ActionIntegrationResultDecorator();
}
#endif /* end of include guard: ACTIONINTEGRATIONRESULTDECORATOR_HPP_ETQNS6TA */
