#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <string>
#include "hamiltonian_dynamic_system.hpp"
#include "state_bindings.hpp"
#include "action_integration_bindings.hpp"

namespace p = boost::python;
namespace np = boost::python::numpy;

namespace { // Avoid cluttering the global namespace.




}






BOOST_PYTHON_MODULE(action_integration_ext)
{
  np::initialize();  // have to put this in any module that uses Boost.NumPy
  export_args_to_nd_array();
  export_args_to_Hessian();
  ActionIntegrationBindings::export_ActionIntegrationResultDecorator();
  ActionIntegrationBindings::export_integrate_E_H_O();
  ActionIntegrationBindings::export_integrate_E_Pendulum();
  ActionIntegrationBindings::export_IntegrationOptions();
  StateBindings::ArmaSB::export_iterable_to_ndarray_for_testing();


}

