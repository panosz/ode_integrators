#include <exception>
#include "action_integration_bindings.hpp"
#include "hamiltonian_dynamic_system.hpp"
#include "action_integration_result.hpp"

namespace ActionIntegrationBindings{

  template<typename Hamiltonian>
  ActionIntegrationResultDecorator
  action_integration_binding_impl(const np::ndarray& ndar,
                                  double mass,
                                  double integration_time,
                                  IntegrationOptions options)
  {
    const auto myHam = Hamiltonian(mass);
    auto my_sys = DS::makeActionDynamicSystem(myHam);
    const auto s = StateBindings::ndarray_to_phase_space_state(ndar);


    try
    {
      const auto res = action_integration(my_sys, s, integration_time, options);
    }
    catch(...)
    {
      throw std::logic_error("just integrated");
    }

    return action_integration(my_sys, s, integration_time, options);


  }




  ActionIntegrationResultDecorator integrate_E_H_O(const np::ndarray& ndar,
                                                   double mass,
                                                   double integration_time,
                                                   IntegrationOptions options)
  {
    return action_integration_binding_impl<DS::UnperturbedExtendedOscillatorHamiltonian>
      (ndar, mass, integration_time, options);
  }



  ActionIntegrationResultDecorator integrate_E_Pendulum(const np::ndarray& ndar,
                                                        double mass,
                                                        double integration_time,
                                                        IntegrationOptions options)
  {
    return action_integration_binding_impl<DS::UnperturbedExtendedPendulumHamiltonian>
      (ndar, mass, integration_time, options);
  }


/**************************************************
*  ActionIntegrationResultDecorator definitions  *
**************************************************/

/*********************************************************
*  end of ActionIntegrationResultDecorator definitions  *
*********************************************************/


  void export_integrate_E_H_O()
  {
    p::def("integrate_E_H_O",
            integrate_E_H_O,
            (p::arg("s"),
             p::arg("mass"),
             p::arg("integration_time"),
             p::arg("integration_options")));
  }


  void export_integrate_E_Pendulum()
  {
    p::def("integrate_E_Pendulum",
            integrate_E_Pendulum,
            (p::arg("s"),
             p::arg("mass"),
             p::arg("integration_time"),
             p::arg("integration_options")));
  }
}

// export bindings

namespace p = boost::python;
namespace np = boost::python::numpy;
