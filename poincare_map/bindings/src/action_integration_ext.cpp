#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <string>
#include "action_integration_result.hpp"
#include "action_integration.hpp"
#include "hamiltonian_dynamic_system.hpp"

namespace p = boost::python;
namespace np = boost::python::numpy;

namespace { // Avoid cluttering the global namespace.

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

  ActionIntegrationResult integrate_E_H_O_impl(const np::ndarray& ndar, double mass, double integration_time)
  {

  const auto options = IntegrationOptions(1e-12, 1e-12, 1e-5);
  const auto myHam = DS::UnperturbedExtendedOscillatorHamiltonian(mass);
  auto my_sys = DS::makeActionDynamicSystem(myHam);
  const auto s = ndarray_to_phase_space_state(ndar);

  return action_integration(my_sys, s, integration_time, options);

  }


  ActionIntegrationResult integrate_E_Pendulum_impl(const np::ndarray& ndar, double mass, double integration_time)
  {

  const auto options = IntegrationOptions(1e-12, 1e-12, 1e-5);
  const auto myHam = DS::UnperturbedExtendedPendulumHamiltonian(mass);
  auto my_sys = DS::makeActionDynamicSystem(myHam);
  const auto s = ndarray_to_phase_space_state(ndar);

  return action_integration(my_sys, s, integration_time, options);
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


ActionIntegrationResultDecorator integrate_E_H_O(const np::ndarray& ndar, double mass, double integration_time)
{
  return integrate_E_H_O_impl(ndar, mass, integration_time);
}

ActionIntegrationResultDecorator integrate_E_Pendulum(const np::ndarray& ndar, double mass, double integration_time)
{
  return integrate_E_Pendulum_impl(ndar, mass, integration_time);
}

BOOST_PYTHON_MODULE(action_integration_ext)
{
  np::initialize();  // have to put this in any module that uses Boost.NumPy
  p::def("args_to_array", args_to_nd_array);
  p::def("args_to_Hessian", args_to_Hessian);


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

  p::def("integrate_E_H_O",integrate_E_H_O,(p::arg("s"),p::arg("mass"),p::arg("integration_time")=1000));
  p::def("integrate_E_Pendulum",integrate_E_Pendulum,(p::arg("s"),p::arg("mass"),p::arg("integration_time")=1000));
}

