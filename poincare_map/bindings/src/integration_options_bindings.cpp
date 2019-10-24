#include <boost/python.hpp>
#include "integration_options_bindings.hpp"

void IntegrationOptionsBindings::export_IntegrationOptions()
{
  namespace p=boost::python;

  const auto integrationOptions_docstring =\
             """Holds options for the integration utilities\n"""
             "\n"
             "IntegrationOptions(abs_err,rel_err,init_time_step)\n"
             "\n"
             "Parameters\n"
             "----------\n"
             "abs_err: float\n"
             "max absolute error for the integration scheme\n"
             "rel_error: float\n"
             "max relative error for the integration scheme\n"
             "init_time_step: float\n"
             "initial time step for the integration scheme\n"
             "\n"
             "Attributes\n"
             "----------\n"
             "abs_err: float\n"
             "max absolute error for the integration scheme\n"
             "rel_error: float\n"
             "max relative error for the integration scheme\n"
             "dt: float\n"
             "initial time step for the integration scheme\n"
             """";

  p::class_<IntegrationOptions>("IntegrationOptions",
      integrationOptions_docstring,
      p::init<double,double,double>(
        (p::arg("abs_err"),
         p::arg("rel_err"),
         p::arg("init_time_step"))))
    .def_readwrite("abs_err",&IntegrationOptions::abs_err)
    .def_readwrite("rel_err",&IntegrationOptions::rel_err)
    .def_readwrite("dt",&IntegrationOptions::dt);
}



