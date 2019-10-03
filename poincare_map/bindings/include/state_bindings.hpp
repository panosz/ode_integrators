#ifndef PANOS_STATE_BINDINGS_HPP
#define PANOS_STATE_BINDINGS_HPP
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include "hamiltonian_dynamic_system.hpp"


void export_args_to_nd_array();
void export_args_to_Hessian();

namespace StateBindings {

  namespace p = boost::python;
  namespace np = boost::python::numpy;

  DS::PhaseSpaceState ndarray_to_phase_space_state(const np::ndarray& ndar);

  np::ndarray args_to_nd_array(double, double, double);
  np::ndarray args_to_Hessian(double xx, double xy, double yy);

}
#endif /* ifndef PANOS_STATE_BINDINGS_HPP */
