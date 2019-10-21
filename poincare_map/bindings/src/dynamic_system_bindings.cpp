#include "dynamic_system_bindings.hpp"
#include "action_integration.hpp"
#include "ActionIntegrationResultDecorator.hpp"

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

}

namespace DynamicSystemBindings{
  using PendulumDynamicSystem =
    DynamicSystem<DS::UnperturbedExtendedPendulumHamiltonian>;
  using HarmonicOscDynamicSystem =
    DynamicSystem<DS::UnperturbedExtendedOscillatorHamiltonian>;

  void export_pendulum_dynamic_system()
  {
    p::class_<PendulumDynamicSystem>("PendulumDynamicSystem",
        p::init<double>())
      .def("action_integrals",&PendulumDynamicSystem::action_integrals)
      .def("closed_orbit",&PendulumDynamicSystem::closed_orbit);
  }

  void export_harmonic_osc_dynamic_system()
  {
    p::class_<HarmonicOscDynamicSystem>("HarmonicOscDynamicSystem",
        p::init<double>())
      .def("action_integrals",&HarmonicOscDynamicSystem::action_integrals)
      .def("closed_orbit",&HarmonicOscDynamicSystem::closed_orbit);
  }
}
