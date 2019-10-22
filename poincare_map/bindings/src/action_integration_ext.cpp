#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <string>
#include "hamiltonian_dynamic_system.hpp"
#include "state_bindings.hpp"
#include "integration_options_bindings.hpp"
#include "dynamic_system_bindings.hpp"
#include "ActionIntegrationResultDecorator.hpp"
#include "hamiltonians_bindings.hpp"

namespace p = boost::python;
namespace np = boost::python::numpy;



BOOST_PYTHON_MODULE(action_integration_ext)
{
  np::initialize();  // have to put this in any module that uses Boost.NumPy
  StateBindings::export_args_to_nd_array();
  StateBindings::export_args_to_Hessian();

  IntegrationOptionsBindings::export_IntegrationOptions();
  ActionIntegrationBindings::export_ActionIntegrationResultDecorator();
  StateBindings::ArmaSB::export_iterable_to_ndarray_for_testing();

  HamiltoniansBindings::export_hamiltonians();

  DynamicSystemBindings::export_dynamic_systems();



}

