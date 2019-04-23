//
// Created by Panagiotis Zestanakis on 23/04/19.
//

#ifndef ODE_INTEGRATORS_HAMILTONIAN_DYNAMIC_SYSTEM_HPP
#define ODE_INTEGRATORS_HAMILTONIAN_DYNAMIC_SYSTEM_HPP
#include <cmath>
#include "armadillo_state.hpp"

namespace DS
{

    using myState = armadillo_state<6>;

    struct OneForm {
        double p;
        double q;
        OneForm (double P, double Q)
            : p(P), q(Q)
        { }
    };

    class UnperturbedExtendedPendulumHamiltonian {

     private:
      double M_;

     public:
      using StateType=myState;

      explicit UnperturbedExtendedPendulumHamiltonian (double M)
          : M_{M}
      { };

      double operator() (const myState& s) const
      {
        const auto& p = s[0];
        const auto& q = s[1];
        const auto& F = s[2];
        return M_ * p * p / 2 - F * cos(q);
      }

      double dp (const myState& s) const noexcept
      {
        const auto& p = s[0];
        return M_ * p;
      }

      double dq (const myState& s) const noexcept
      {
        const auto& q = s[1];
        const auto& F = s[2];
        return F * sin(q);
      }

      double dF (const myState& s) const noexcept
      {
        const auto& q = s[1];
        return -cos(q);
      }

    };

    template<typename UnperturbedHamiltonian>
    class UnperturbedDynamicSystem {

     private:
      UnperturbedHamiltonian h_;
     public:
      using StateType = typename UnperturbedHamiltonian::StateType;

      explicit UnperturbedDynamicSystem (const UnperturbedHamiltonian& h)
          : h_{h}
      { }

      void operator() (const StateType& s, StateType& dsdt, const double /*t*/) const
      {
        const auto& p = s[0];

        dsdt[0] = -h_.dq(s);
        dsdt[1] = h_.dp(s);
        dsdt[2] = 0;
        dsdt[3] = h_.dF(s);
        dsdt[4] = p *dsdt[1];
        dsdt[5] = 1;
      }
    };

    template<typename UnperturbedHamiltonian>
    UnperturbedDynamicSystem<UnperturbedHamiltonian> makeUnperturbedDynamicSystem(const UnperturbedHamiltonian& h)
    {
      return UnperturbedDynamicSystem<UnperturbedHamiltonian>(h);
    }




}
#endif //ODE_INTEGRATORS_HAMILTONIAN_DYNAMIC_SYSTEM_HPP
