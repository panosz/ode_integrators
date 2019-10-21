#ifndef PANOS_ACTION_INTEGRATION_BINDINGS_HPP
#define PANOS_ACTION_INTEGRATION_BINDINGS_HPP
#include "action_integration.hpp"
#include "ActionIntegrationResultDecorator.hpp"


namespace ActionIntegrationBindings{

  namespace p = boost::python;
  namespace np = boost::python::numpy;




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


}


#endif /* ifndef PANOS_ACTION_INTEGRATION_BINDINGS_HPP */
