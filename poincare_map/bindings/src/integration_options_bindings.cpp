#include <boost/python.hpp>
#include "integration_options_bindings.hpp"

void IntegrationOptionsBindings::export_IntegrationOptions()
{
  namespace p=boost::python;

  const auto integrationOptions_docstring =\
             """Holds options for the integration utilities\n"""
             "\n"
             "IntegrationOptions(abs_err,rel_err,dt,coefficient)\n"
             "\n"
             "Parameters\n"
             "----------\n"
             "abs_err: float\n"
             "max absolute error for the integration scheme\n"
             "rel_error: float\n"
             "max relative error for the integration scheme\n"
             "dt: float\n"
             "initial time step for the integration scheme\n"
             "coefficient: float\n"
             "tolerance coefficient. Relevant for the closing orbit detection\n"
             "strategies. Multiplies the integration tolerance to determine\n"
             "the orbit closing tolerance\n"
             "\n"
             "Attributes\n"
             "----------\n"
             "abs_err: float\n"
             "max absolute error for the integration scheme\n"
             "rel_error: float\n"
             "max relative error for the integration scheme\n"
             "dt: float\n"
             "initial time step for the integration scheme\n"
             "coefficient: float\n"
             "tolerance coefficient. Relevant for the closing orbit detection\n"
             """";

  p::class_<IntegrationOptions>("IntegrationOptions",
      integrationOptions_docstring,
      p::init<double,double,double,double>(
        (p::arg("abs_err"),
         p::arg("rel_err"),
         p::arg("dt"),
         p::arg("coefficient"))))
    .def(p::init<double,double,double>(
        (p::arg("abs_err"),
         p::arg("rel_err"),
         p::arg("dt"))))
    .def(p::init<double,double>(
        (p::arg("abs_err"),
         p::arg("rel_err"))))
    .def(p::init<>())
    .def_readwrite("abs_err",&IntegrationOptions::abs_err)
    .def_readwrite("rel_err",&IntegrationOptions::rel_err)
    .def_readwrite("dt",&IntegrationOptions::dt)
    .def_readwrite("coefficient",&IntegrationOptions::coefficient);
}



