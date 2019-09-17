#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <string>
#include "action_integration_result.hpp"

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
    {return args_to_Hessian(air_.d2K_dJ2(), air_.d2K_dJdF(), air_.d2K_dF2());}

};


ActionIntegrationResultDecorator make_action_integration_result()
{
  return ActionIntegrationResultDecorator(1, 2, 3, SpecialIntegrals());
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

  p::def("make_action_integration_result",make_action_integration_result);
}

