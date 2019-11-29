#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <armadillo>
#include <iostream>
#include "state_bindings.hpp"

namespace{

  using namespace StateBindings;


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
  np::ndarray copy_Arma_to_nd_array_naive(const std::vector<ArmaState>& A)
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
  std::vector<double> iterable_to_vector_double(const p::object& python_iterable)
  {

    auto input_begin = p::stl_input_iterator<double>(python_iterable);
    auto input_end = p::stl_input_iterator<double>();

    const auto output_vector = std::vector<double>(input_begin,input_end);

    return output_vector;
  }

  np::ndarray copy_to_nd_array(const std::vector<double>& vect)
  {
    // construct a type for C++ double
    const np::dtype dtype = np::dtype::get_builtin<double>();
    const p::tuple shape = p::make_tuple(vect.size());

    //construct uninitialized array
    np::ndarray arr = np::empty(shape, dtype);

    boost::range::copy(vect,
        reinterpret_cast<double*>(arr.get_data()));

    return arr;
  }

  namespace Testing
  {
    namespace p = StateBindings::p;
    namespace np = StateBindings::np;
      np::ndarray iterable_1D_to_ndarray(const p::object& python_iterable)
      {
        const auto helper_vector = iterable_to_vector_double(python_iterable);
        return copy_to_nd_array(helper_vector);
      }

      const char* iterable_1D_to_ndarray_doc =
               """iterable_1D_to_ndarray(iterable)\n"
               "create a 1D ndarray from 'iterable' via intermediate construction\n"
               "of an 'std::vector<double>' instance.\n"
               "To be used to test copying from and to std::vector\n"
               "\n"
               "Parameters\n"
               "----------\n"
               "iterable: object, iterable\n"
               "\t must be an object convertible to a 1D array\n"
               "\t i.e. an iterable of doubles\n"
               "Retruns\n"
               "-------\n"
               "ndarray\n"
               """";
  }


  namespace ArmaSB {
    namespace p = StateBindings::p;
    namespace np = StateBindings::np;


    np::ndarray copy_to_nd_array(const std::vector<arma::mat>& A)
    {
      return copy_Arma_to_nd_array_naive<arma::mat>(A);
    }

    np::ndarray copy_to_nd_array(const std::vector<DS::ExtendedSpaceState>& vds)
    {
      return copy_Arma_to_nd_array_naive<DS::ExtendedSpaceState>(vds);
    }

    np::ndarray copy_to_nd_array(const std::vector<DS::PhaseSpaceState>& vds)
    {
      return copy_Arma_to_nd_array_naive<DS::PhaseSpaceState>(vds);
    }


    namespace Testing
    {
      np::ndarray iterable_2D_to_ndarray(const ArmaSB::p::object& python_iterable)
      {
        const auto helper_vector = iterable_to_vector_of_arma_mat(python_iterable);
        return copy_to_nd_array(helper_vector);
      }

      const char* iterable_2D_to_ndarray_doc =
               """iterable_2D_to_ndarray(iterable)\n"
               "create a 2D ndarray from 'iterable' via intermediate construction\n"
               "of an 'Arma::mat' instance.\n"
               "To be used to test copying from and to Arma::mat\n"
               "\n"
               "Parameters\n"
               "----------\n"
               "iterable: object, iterable\n"
               "\t must be an object convertible to a 2D array\n"
               "\t i.e. an iterable of same length iterables\n"
               "Retruns\n"
               "-------\n"
               "ndarray\n"
               """";
    }


  }

  void export_methods_for_testing()
  {
    p::def("iterable_2D_to_ndarray",
            ArmaSB::Testing::iterable_2D_to_ndarray,
            p::args("iterable"),
            ArmaSB::Testing::iterable_2D_to_ndarray_doc
            );

    p::def("iterable_1D_to_ndarray",
            Testing::iterable_1D_to_ndarray,
            p::args("iterable"),
            Testing::iterable_1D_to_ndarray_doc
            );
  }
}
