#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <string>
#include "action_integration_result.hpp"
#include "action_integration.hpp"
#include "hamiltonian_dynamic_system.hpp"
#include "state_bindings.hpp"

namespace p = boost::python;
namespace np = boost::python::numpy;

namespace { // Avoid cluttering the global namespace.


  ActionIntegrationResult integrate_E_H_O_impl(const np::ndarray& ndar, double mass, double integration_time)
  {

  const auto options = IntegrationOptions(1e-12, 1e-12, 1e-5);
  const auto myHam = DS::UnperturbedExtendedOscillatorHamiltonian(mass);
  auto my_sys = DS::makeUnperturbedDynamicSystem(myHam);
  const auto s = StateBindings::ndarray_to_phase_space_state(ndar);

  return action_integration(my_sys, s, integration_time, options);

  }


  ActionIntegrationResult integrate_E_Pendulum_impl(const np::ndarray& ndar, double mass, double integration_time)
  {

  const auto options = IntegrationOptions(1e-12, 1e-12, 1e-5);
  const auto myHam = DS::UnperturbedExtendedPendulumHamiltonian(mass);
  auto my_sys = DS::makeUnperturbedDynamicSystem(myHam);
  const auto s = StateBindings::ndarray_to_phase_space_state(ndar);

  return action_integration(my_sys, s, integration_time, options);
  }



}



class ActionIntegrationResultDecorator
{
  private:
    ActionIntegrationResult air_;
  public:
    ActionIntegrationResultDecorator(const ActionIntegrationResult& air):
      air_{air}{};
    ActionIntegrationResultDecorator(double x, double y, double z, SpecialIntegrals si):
      air_{x,y,z,si}{};
    double Action () const noexcept
    {return air_.Action();}
    double theta_period () const noexcept
    {return air_.theta_period();}
    double delta_phi () const noexcept
    {return air_.delta_phi();}
    double omega_theta () const
    {return air_.omega_theta();}
    double g_factor () const noexcept
    {return air_.g_factor();}
    double omega_phi () const
    {return air_.omega_phi();}
    double domega_dJ () const
    {return air_.domega_dJ();}
    double d2K_dJ2 () const
    {return air_.d2K_dJ2();}
    double one_div_two_pi_gamma () const
    {return air_.one_div_two_pi_gamma();}
    double domega_dF () const
    {return air_.domega_dF();}
    double d2K_dJdF () const
    {return air_.d2K_dJdF();}
    double d2K_dF2 () const
    {return air_.d2K_dF2();}
    np::ndarray hessian() const
    {return StateBindings::args_to_Hessian(air_.d2K_dJ2(), air_.d2K_dJdF(), air_.d2K_dF2());}


};


ActionIntegrationResultDecorator integrate_E_H_O(const np::ndarray& ndar, double mass, double integration_time)
{
  return integrate_E_H_O_impl(ndar, mass, integration_time);
}

ActionIntegrationResultDecorator integrate_E_Pendulum(const np::ndarray& ndar, double mass, double integration_time)
{
  return integrate_E_Pendulum_impl(ndar, mass, integration_time);
}

BOOST_PYTHON_MODULE(action_integration_ext)
{
  np::initialize();  // have to put this in any module that uses Boost.NumPy
  export_args_to_nd_array();
  export_args_to_Hessian();


  p::class_<ActionIntegrationResultDecorator>("ActionIntegrationResult", p::init<double,double,double,SpecialIntegrals>())
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

  p::def("integrate_E_H_O",integrate_E_H_O,(p::arg("s"),p::arg("mass"),p::arg("integration_time")=1000));
  p::def("integrate_E_Pendulum",integrate_E_Pendulum,(p::arg("s"),p::arg("mass"),p::arg("integration_time")=1000));
}

