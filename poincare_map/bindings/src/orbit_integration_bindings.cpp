#include "orbit_integration_bindings.hpp"
#include "hamiltonian_dynamic_system.hpp"

namespace OrbitIntegrationBindings
{


  template<typename Hamiltonian>
  np::ndarray
  integrate_along_closed_orbit_binding_impl(const np::ndarray& starting_point,
                                            double mass,
                                            double integration_time,
                                            IntegrationOptions options)
  {
    const auto myHam = Hamiltonian(mass);
    auto my_sys = DS::makeUnperturbedDynamicSystem(myHam);
    const auto s = StateBindings::ndarray_to_phase_space_state(starting_point);

    const auto orbit = integrate_along_closed_orbit(my_sys,
                                                    s,
                                                    integration_time,
                                                    options);

    return StateBindings::ArmaSB::vector_of_arma_dynamic_state_to_nd_array_naive(orbit);
  }

  np::ndarray
  integrate_along_closed_harmonic_osc_orbit_binding(const np::ndarray& starting_point,
                                                    double mass,
                                                    double integration_time,
                                                    IntegrationOptions options)
  {
    using Ham = DS::UnperturbedExtendedOscillatorHamiltonian;

    return integrate_along_closed_orbit_binding_impl<Ham> (starting_point,
                                                           mass,
                                                           integration_time,
                                                           options);
  }


  np::ndarray
  integrate_along_closed_pendulum_orbit_binding(const np::ndarray& starting_point,
                                                double mass,
                                                double integration_time,
                                                IntegrationOptions options)
  {
    using Ham = DS::UnperturbedExtendedPendulumHamiltonian;

    return integrate_along_closed_orbit_binding_impl<Ham> (starting_point,
                                                           mass,
                                                           integration_time,
                                                           options);
  }

  void export_closed_harmonic_osc_orbit()
  {
    p::def("closed_harmonic_osc_orbit",
            integrate_along_closed_harmonic_osc_orbit_binding,
            (p::arg("s"),
             p::arg("mass"),
             p::arg("integration_time"),
             p::arg("integration_options")));
  }


  void export_closed_pendulum_orbit()
  {
    p::def("closed_pendulum_orbit",
            integrate_along_closed_pendulum_orbit_binding,
            (p::arg("s"),
             p::arg("mass"),
             p::arg("integration_time"),
             p::arg("integration_options")));
  }


}
