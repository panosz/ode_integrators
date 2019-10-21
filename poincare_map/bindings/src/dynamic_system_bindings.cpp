#include "dynamic_system_bindings.hpp"
#include "action_integration.hpp"
#include "ActionIntegrationResultDecorator.hpp"
#include "hamiltonians_bindings.hpp"

namespace{
  namespace p = boost::python;
  namespace np = boost::python::numpy;

  template<typename Hamiltonian>
  class DynamicSystem
  {
    private:
        Hamiltonian ham_{};
        DS::ActionDynamicSystem<Hamiltonian> action_system_{};
        DS::UnperturbedDynamicSystem<Hamiltonian> unperturbed_system_{};
        using HDec = HamiltoniansBindings::Decorator<Hamiltonian>;

    public:
        explicit DynamicSystem(double mass):
          ham_{mass},
          action_system_{ham_},
          unperturbed_system_{ham_}
        {};

        ActionIntegrationBindings::ActionIntegrationResultDecorator
        action_integrals(const np::ndarray& ndar,
                         double integration_time,
                         IntegrationOptions options) const
        {
          const auto s = StateBindings::ndarray_to_phase_space_state(ndar);
          return action_integration(action_system_, s, integration_time, options);
        }

        HDec get_ham() const noexcept
        {
          return HDec{ham_};
        }


        np::ndarray
        closed_orbit(const np::ndarray& starting_point,
                     double integration_time,
                     IntegrationOptions options) const
        {
          const auto input_state =
            StateBindings::ndarray_to_phase_space_state(starting_point);

          const auto orbit = integrate_along_closed_orbit(unperturbed_system_,
                                                          input_state,
                                                          integration_time,
                                                          options);

          return StateBindings::ArmaSB::
                vector_of_arma_dynamic_state_to_nd_array_naive(orbit);
        }

  };

  const auto pendulum_docstring =
             """The unpertrurbed extended pendulum class\n"""
             "\n"
             "PengulumDynamicSystem(mass)\n"
             "\n"
             "Parameters\n"
             "----------\n"
             "mass: float\n"
             "\tthe mass of the system\n"
             """";

  const auto harmonic_osc_docstring =
             """The unpertrurbed extended harmonic oscillator class\n"""
             "\n"
             "HarmonicOscDynamicSystem(mass)\n"
             "\n"
             "Parameters\n"
             "----------\n"
             "mass: float\n"
             "\tthe mass of the system\n"
             """";

  const auto action_integrals_docstring =
             """action_integrals(s, time, options)\n"
             "calculate the action integrals along a closed orbit\n"
             "\n"
             "Parameters\n"
             "----------\n"
             "s: ndarray,\n"
             "\tthe starting point of the orbit\n"
             "time: ndarray,\n"
             "\tthe maximum integration time\n"
             "\tif in is exceeded before the orbit closes, an exception is thrown\n"
             "options: IntegrationOptions,\n"
             "\tthe integration options\n"
             """";

  const auto closed_orbit_docstring =
             """closed_orbit(s, time, options)\n"
             "calculate the closed orbit passing from point 's'\n"
             "\n"
             "Parameters\n"
             "----------\n"
             "s: ndarray,\n"
             "\tthe starting point of the orbit\n"
             "time: ndarray,\n"
             "\tthe maximum integration time\n"
             "\tif in is exceeded before the orbit closes, an exception is thrown\n"
             "options: IntegrationOptions,\n"
             "\tthe integration options\n"
             """";

}

namespace DynamicSystemBindings{
  using PendulumDynamicSystem =
    DynamicSystem<DS::UnperturbedExtendedPendulumHamiltonian>;
  using HarmonicOscDynamicSystem =
    DynamicSystem<DS::UnperturbedExtendedOscillatorHamiltonian>;

  void export_pendulum_dynamic_system()
  {
    p::class_<PendulumDynamicSystem>("PendulumDynamicSystem",
                                     pendulum_docstring,
                                     p::init<double>())
      .add_property("hamiltonian",&PendulumDynamicSystem::get_ham)
      .def("action_integrals",
           &PendulumDynamicSystem::action_integrals,
           (p::arg("s"),
            p::arg("time"),
            p::arg("options")),
           action_integrals_docstring
           )
      .def("closed_orbit",
           &PendulumDynamicSystem::closed_orbit,
           (p::arg("s"),
            p::arg("time"),
            p::arg("options")),
           closed_orbit_docstring);
  }

  void export_harmonic_osc_dynamic_system()
  {
    p::class_<HarmonicOscDynamicSystem>("HarmonicOscDynamicSystem",
                                        harmonic_osc_docstring,
                                        p::init<double>())
      .add_property("hamiltonian",&HarmonicOscDynamicSystem::get_ham)
      .def("action_integrals",
           &HarmonicOscDynamicSystem::action_integrals,
           (p::arg("s"),
            p::arg("time"),
            p::arg("options")),
           action_integrals_docstring)
      .def("closed_orbit",
           &HarmonicOscDynamicSystem::closed_orbit,
           (p::arg("s"),
            p::arg("time"),
            p::arg("options")),
           closed_orbit_docstring);
  }
}
