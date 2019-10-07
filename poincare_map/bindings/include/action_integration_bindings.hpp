#ifndef PANOS_ACTION_INTEGRATION_BINDINGS_HPP
#define PANOS_ACTION_INTEGRATION_BINDINGS_HPP
#include "action_integration.hpp"
#include "action_integration_result.hpp"
#include "state_bindings.hpp"


namespace ActionIntegrationBindings{

  namespace p = boost::python;
  namespace np = boost::python::numpy;

  ///ToDo: Provide interface for IntegrationOptions
  const auto DefaultIntegrationOptions = IntegrationOptions(1e-12, 1e-12, 1e-5);
  const double DefaultIntegrationTime = 1000;


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
      np::ndarray hessian() const;


  };

  ActionIntegrationResultDecorator integrate_E_H_O(const np::ndarray& ndar,
                                                   double mass,
                                                   double integration_time,
                                                   IntegrationOptions options);

  ActionIntegrationResultDecorator integrate_E_Pendulum(const np::ndarray& ndar,
                                                        double mass,
                                                        double integration_time,
                                                        IntegrationOptions options);
  void export_ActionIntegrationResultDecorator();
  void export_integrate_E_H_O();
  void export_integrate_E_Pendulum();
  void export_IntegrationOptions();


}


#endif /* ifndef PANOS_ACTION_INTEGRATION_BINDINGS_HPP */
