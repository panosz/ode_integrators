#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <armadillo>
#include <iostream>
#include "state_bindings.hpp"

namespace{

  using namespace StateBindings;

  std::vector<double> iterable_to_vector_double(const p::object& python_iterable)
  {

    auto input_begin = p::stl_input_iterator<double>(python_iterable);
    auto input_end = p::stl_input_iterator<double>();

    const auto output_vector = std::vector<double>(input_begin,input_end);

    return output_vector;
  }

  std::vector<arma::mat> iterable_to_vector_of_arma_mat(const p::object& python_iterable)
  {
    std::vector<arma::mat> output{};
    auto input_begin = p::stl_input_iterator<p::object>(python_iterable);
    auto input_end = p::stl_input_iterator<p::object>();

    for (auto i = input_begin; i != input_end; ++i)
    {
      output.push_back(iterable_to_vector_double(*i));
    }

    return output;

  }


  template<typename ArmaState>
  np::ndarray vector_of_Arma_State_to_nd_array_naive(const std::vector<ArmaState>& A)
  {

    // construct a type for C++ double
    const np::dtype dtype = np::dtype::get_builtin<double>();

    const auto n_rows = A.size();

    if (n_rows == 0) {
      const p::tuple shape = p::make_tuple(0);
      return np::empty(shape, dtype);
    }

    const auto n_cols = A[0].size();

    if (n_cols == 0) {
      const p::tuple shape = p::make_tuple(n_rows,0);
      return np::empty(shape, dtype);
    }

    //construct the shape as a tuple.
    const p::tuple shape = p::make_tuple(n_rows, n_cols);

    np::ndarray a = np::empty(shape, dtype);

    for (arma::uword i = 0; i < n_rows; ++i)
    {
      for (arma::uword j = 0; j < n_cols; ++j)
      {
        a[i][j] = A[i].at(j);
      }
    }
    return a;
  }

}


namespace StateBindings {

  namespace ArmaSB {
    namespace p = StateBindings::p;
    namespace np = StateBindings::np;


    np::ndarray vector_of_arma_mat_to_nd_array_naive(const std::vector<arma::mat>& A)
    {
      return vector_of_Arma_State_to_nd_array_naive<arma::mat>(A);
    }

    np::ndarray vector_of_arma_dynamic_state_to_nd_array_naive(const std::vector<DS::ExtendedSpaceState>& vds)
    {
      return vector_of_Arma_State_to_nd_array_naive<DS::ExtendedSpaceState>(vds);
    }

    np::ndarray vector_of_arma_dynamic_state_to_nd_array_naive(const std::vector<DS::PhaseSpaceState>& vds)
    {
      return vector_of_Arma_State_to_nd_array_naive<DS::PhaseSpaceState>(vds);
    }


    np::ndarray iterable_to_ndarray_for_testing(const p::object& python_iterable)
    {
      const auto helper_vector = iterable_to_vector_of_arma_mat(python_iterable);
      return vector_of_arma_mat_to_nd_array_naive(helper_vector);
    }

    /*********************
    *  export bindings  *
    *********************/

    void export_iterable_to_ndarray_for_testing()
    {
      p::def("iterable_to_ndarray_for_testing", iterable_to_ndarray_for_testing);
    }

  }
}
