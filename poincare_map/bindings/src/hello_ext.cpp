#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <string>

namespace p = boost::python;
namespace np = boost::python::numpy;

namespace { // Avoid cluttering the global namespace.
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

  class array_builder
  {
    public:
      array_builder() = default;
      np::ndarray to_array(double x, double y, double z) const
      {
        return args_to_nd_array(x, y, z);
      }
  };


  // A friendly class.
  class MyPoint
  {
  public:
    double x_;
    double y_;
    MyPoint (double x, double y):
      x_{x},y_{y}
    {};
  };

  //A friendly class factory.

  MyPoint make_myPoint(double x, double y)
  {
    return MyPoint(x,y);
  }
}



BOOST_PYTHON_MODULE(hello_ext)
{
  np::initialize();  // have to put this in any module that uses Boost.NumPy
  p::def("args_to_array", args_to_nd_array);
  p::def("args_to_Hessian", args_to_Hessian);

  p::class_<array_builder>("array_builder")
    // Add a regular member function.
    .def("to_array", &array_builder::to_array);

  p::class_<MyPoint>("MyPythonPoint",p::init<double,double>())
    .def_readwrite("x",&MyPoint::x_)
    .def_readwrite("y",&MyPoint::y_) ;
  //exposed factory functions, are python factory functions for the exposed python class
  p::def("make_myPoint",make_myPoint);
}

