#include <boost/python/tuple.hpp>
#include "dynamic_system_bindings.hpp"
#include "integration_utilities.hpp"
#include "action_integration.hpp"
#include "ActionIntegrationResultDecorator.hpp"
#include "hamiltonians_bindings.hpp"
#include "UnperturbedExtendedOscillatorHamiltonian.hpp"
#include "UnperturbedExtendedPendulumHamiltonian.hpp"


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

        auto follow_orbit(const np::ndarray& ndar,
                             double integration_time,
                             IntegrationOptions options) const
        {
          const auto init_s = StateBindings::ndarray_to_phase_space_state(ndar);

          auto orbit = make_ParticleOrbit(unperturbed_system_,
                                          init_s,
                                          integration_time,
                                          options);

          std::vector<double> time{};
          std::vector<DS::PhaseSpaceState> positions{};
          for (const auto & p_t : orbit.time_range())
          {
            auto const &[pos, t] = p_t;
            positions.push_back(pos);
            time.push_back(t);
          }

          return std::make_pair(positions, time);
        }


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
                copy_to_nd_array(orbit);
        }

        p::tuple
          orbit(const np::ndarray& starting_point,
                double integration_time,
                IntegrationOptions options) const
          {
            const auto [vector, time] = follow_orbit(starting_point,
                                                     integration_time,
                                                     options);

            return p::make_tuple(StateBindings::ArmaSB::copy_to_nd_array(vector),
                                 StateBindings::copy_to_nd_array(time));
          }

        auto orbit_at_times(const np::ndarray& ndar,
                           const p::object& times_iterable,
                           IntegrationOptions options) const
        {
          const auto init_s = StateBindings::ndarray_to_phase_space_state(ndar);
          const auto times_vector =
            StateBindings::iterable_to_vector_double(times_iterable);

          const auto points =  orbit_points_at_times(unperturbed_system_,
                                                     init_s,
                                                     times_vector,
                                                     options);
          return StateBindings::ArmaSB:: copy_to_nd_array(points);

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


template<typename Hamiltonian>
void export_dynamic_system(const char* system_name,
                           const char* system_docstring)
{
  using DSH = DynamicSystem<Hamiltonian>;

    p::class_<DSH>(system_name,
                   system_docstring,
                   p::init<double>())
      .add_property("hamiltonian",&DSH::get_ham)
      .def("action_integrals",
           &DSH::action_integrals,
           (p::arg("s"),
            p::arg("time"),
            p::arg("options")),
           action_integrals_docstring
           )
      .def("closed_orbit",
           &DSH::closed_orbit,
           (p::arg("s"),
            p::arg("time"),
            p::arg("options")),
           closed_orbit_docstring)
      .def("orbit",
           &DSH::orbit,
           (p::arg("s"),
            p::arg("time"),
            p::arg("options")))
      .def("orbit_at_times",
           &DSH::orbit_at_times,
           (p::arg("s"),
            p::arg("times"),
            p::arg("options")));

}

namespace DynamicSystemBindings{
  using PendulumHam = DS::UnperturbedExtendedPendulumHamiltonian;
  using HarmonicOscHam= DS::UnperturbedExtendedOscillatorHamiltonian;

  void export_dynamic_systems()
  {
    export_dynamic_system<PendulumHam>("PendulumDynamicSystem",
                                       pendulum_docstring);
    export_dynamic_system<HarmonicOscHam>("HarmonicOscDynamicSystem",
                                          harmonic_osc_docstring);
  }
}
