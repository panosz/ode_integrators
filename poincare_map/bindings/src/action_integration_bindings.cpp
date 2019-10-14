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
  ActionIntegrationResultDecorator::
    ActionIntegrationResultDecorator(const ActionIntegrationResult& air):
      air_{air}{}

  ActionIntegrationResultDecorator::
    ActionIntegrationResultDecorator(double x, double y, double z, SpecialIntegrals si):
      air_{x,y,z,si}{}

  double ActionIntegrationResultDecorator::Action () const noexcept
  {return air_.Action();}

  double ActionIntegrationResultDecorator::theta_period () const noexcept
  {return air_.theta_period();}

  double ActionIntegrationResultDecorator::delta_phi () const noexcept
  {return air_.delta_phi();}

  double ActionIntegrationResultDecorator::omega_theta () const
  {return air_.omega_theta();}

  double ActionIntegrationResultDecorator::g_factor () const noexcept
  {return air_.g_factor();}

  double ActionIntegrationResultDecorator::omega_phi () const
  {return air_.omega_phi();}

  double ActionIntegrationResultDecorator::domega_dJ () const
  {return air_.domega_dJ();}

  double ActionIntegrationResultDecorator::d2K_dJ2 () const
  {return air_.d2K_dJ2();}

  double ActionIntegrationResultDecorator::one_div_two_pi_gamma () const
  {return air_.one_div_two_pi_gamma();}

  double ActionIntegrationResultDecorator::domega_dF () const
  {return air_.domega_dF();}

  double ActionIntegrationResultDecorator::d2K_dJdF () const
  {return air_.d2K_dJdF();}

  double ActionIntegrationResultDecorator::d2K_dF2 () const
  {return air_.d2K_dF2();}

  np::ndarray ActionIntegrationResultDecorator::hessian() const
  {
    return StateBindings::args_to_Hessian(air_.d2K_dJ2(),
                                          air_.d2K_dJdF(),
                                          air_.d2K_dF2());
  }
/*********************************************************
*  end of ActionIntegrationResultDecorator definitions  *
*********************************************************/

  void export_ActionIntegrationResultDecorator()
  {
    p::class_<ActionIntegrationResultDecorator>("ActionIntegrationResult",
        p::init<double,double,double,SpecialIntegrals>())
      .def("Action",&ActionIntegrationResultDecorator::Action)
      .def("theta_period",&ActionIntegrationResultDecorator::theta_period)
      .def("delta_phi",&ActionIntegrationResultDecorator::delta_phi)
      .def("omega_theta",&ActionIntegrationResultDecorator::omega_theta)
      .def("g_factor",&ActionIntegrationResultDecorator::g_factor)
      .def("omega_phi",&ActionIntegrationResultDecorator::omega_phi)
      .def("domega_dJ",&ActionIntegrationResultDecorator::domega_dJ)
      .def("d2K_dJ2",&ActionIntegrationResultDecorator::d2K_dJ2)
      .def("one_div_two_pi_gamma",&ActionIntegrationResultDecorator::one_div_two_pi_gamma)
      .def("domega_dF",&ActionIntegrationResultDecorator::domega_dF)
      .def("d2K_dJdF",&ActionIntegrationResultDecorator::d2K_dJdF)
      .def("d2K_dF2",&ActionIntegrationResultDecorator::d2K_dF2)
      .def("hessian",&ActionIntegrationResultDecorator::hessian);
  }


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
