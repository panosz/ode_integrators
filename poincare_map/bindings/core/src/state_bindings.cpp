#include "state_bindings.hpp"

namespace StateBindings {

  DS::PhaseSpaceState ndarray_to_phase_space_state(const np::ndarray& ndar)
  {

    if (ndar.get_dtype() != np::dtype::get_builtin<double>()) {
        PyErr_SetString(PyExc_TypeError, "Incorrect array data type");
        p::throw_error_already_set();
    }

    if (ndar.get_nd() != 1) {
        PyErr_SetString(PyExc_TypeError, "Incorrect number of array dimensions");
        p::throw_error_already_set();
    }

    if (ndar.shape(0) != 4) {
        PyErr_SetString(PyExc_TypeError, "Incorrect length input");
        p::throw_error_already_set();
    }

  auto s = DS::PhaseSpaceState{};
  s[0] = p::extract<double>(ndar[0]);
  s[1] = p::extract<double>(ndar[1]);
  s[2] = p::extract<double>(ndar[2]);
  s[3] = p::extract<double>(ndar[3]);
  return s;
  }


  std::array<double,3> args_to_array(double x, double y, double z)
  {

    return std::array<double,3>{x,y,z};
  }

  np::ndarray args_to_nd_array(double x, double y, double z)
  {
    //construct the data
    std::array<double,3> std_array = args_to_array(x,y,z);

    //construct the shape as a tuple
    p::tuple shape = p::make_tuple(3);

    // construct a type for C++ double
    np::dtype dtype = np::dtype::get_builtin<double>();

    // Construct an array with the above shape and type
    np::ndarray a = np::zeros(shape, dtype);
    std::copy(std_array.cbegin(), std_array.cend(),
            reinterpret_cast<double *>(a.get_data()));
    // Construct an empty array with the above shape and dtype as well
    // np::ndarray b = np::empty(shape,dtype);
    //

    return a;
  }

  np::ndarray args_to_2_by_2_array(double xx, double xy, double yx, double yy)
  {

    //construct the shape as a tuple
    p::tuple shape = p::make_tuple(2,2);

    // construct a type for C++ double
    np::dtype dtype = np::dtype::get_builtin<double>();

    // Construct an array with the above shape and type
    np::ndarray hes = np::zeros(shape, dtype);

    // Copy the data into the array
    hes[0][0]=xx;
    hes[0][1]=xy;
    hes[1][0]=yx;
    hes[1][1]=yy;

    return hes;
  }

  np::ndarray args_to_Hessian(double xx, double xy, double yy)
  {
    return args_to_2_by_2_array(xx,xy,xy,yy);
  }

  // export bindings
  void export_args_to_nd_array()
  {
    p::def("args_to_array", StateBindings::args_to_nd_array);
  }
  void export_args_to_Hessian()
  {
    p::def("args_to_Hessian", StateBindings::args_to_Hessian);
  }
}



