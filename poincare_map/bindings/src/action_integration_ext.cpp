#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <string>
#include "hamiltonian_dynamic_system.hpp"
#include "state_bindings.hpp"
#include "action_integration_bindings.hpp"
#include "integration_options_bindings.hpp"
#include "orbit_integration_bindings.hpp"
#include "dynamic_system_bindings.hpp"

namespace p = boost::python;
namespace np = boost::python::numpy;



BOOST_PYTHON_MODULE(action_integration_ext)
{
  np::initialize();  // have to put this in any module that uses Boost.NumPy
  export_args_to_nd_array();
  export_args_to_Hessian();
  IntegrationOptionsBindings::export_IntegrationOptions();
  ActionIntegrationBindings::export_ActionIntegrationResultDecorator();
  ActionIntegrationBindings::export_integrate_E_H_O();
  ActionIntegrationBindings::export_integrate_E_Pendulum();
  StateBindings::ArmaSB::export_iterable_to_ndarray_for_testing();

  OrbitIntegrationBindings::export_closed_pendulum_orbit();
  OrbitIntegrationBindings::export_closed_harmonic_osc_orbit();
  DynamicSystemBindings::export_pendulum_dynamic_system();
  DynamicSystemBindings::export_harmonic_osc_dynamic_system();



}

