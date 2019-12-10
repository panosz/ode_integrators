#ifndef HAMILTONIANS_BINDINGS_HPP_SM3ZJ0HX
#define HAMILTONIANS_BINDINGS_HPP_SM3ZJ0HX
#include "hamiltonian_dynamic_system.hpp"
#include "state_bindings.hpp"

namespace HamiltoniansBindings{

  void export_hamiltonians();

  template<typename Hamiltonian>
    class Decorator
    {
      private:
        Hamiltonian ham_;

        inline DS::PhaseSpaceState read_state(const StateBindings::np::ndarray& nd_s) const
        {
          return StateBindings::ndarray_to_phase_space_state(nd_s);
        }

      public:
        explicit Decorator(const Hamiltonian& h):
          ham_{h}{}

        double value(const StateBindings::np::ndarray& nd_s) const
        {
          const auto s = read_state(nd_s);
          return ham_.value(s);
        }
    };

}

#endif /* end of include guard: HAMILTONIANS_BINDINGS_HPP_SM3ZJ0HX */
