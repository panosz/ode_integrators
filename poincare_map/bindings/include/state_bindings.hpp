#ifndef PANOS_STATE_BINDINGS_HPP
#define PANOS_STATE_BINDINGS_HPP
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include "hamiltonian_dynamic_system.hpp"


namespace StateBindings {

  namespace p = boost::python;
  namespace np = boost::python::numpy;

  DS::PhaseSpaceState ndarray_to_phase_space_state(const np::ndarray& ndar);

  np::ndarray args_to_nd_array(double, double, double);
  np::ndarray args_to_Hessian(double xx, double xy, double yy);
  void export_args_to_nd_array();
  void export_args_to_Hessian();


  namespace ArmaSB{


    np::ndarray copy_to_nd_array(const std::vector<arma::mat>& A);

    np::ndarray copy_to_nd_array(const std::vector<DS::ExtendedSpaceState>& vds);

    np::ndarray copy_to_nd_array(const std::vector<DS::PhaseSpaceState>& vds);

    void export_iterable_to_ndarray_for_testing();


  }
}
#endif /* ifndef PANOS_STATE_BINDINGS_HPP */
