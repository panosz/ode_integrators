//
// Created by Panagiotis Zestanakis on 23/04/19.
//

#ifndef ODE_INTEGRATORS_HAMILTONIAN_DYNAMIC_SYSTEM_HPP
#define ODE_INTEGRATORS_HAMILTONIAN_DYNAMIC_SYSTEM_HPP
#include <cmath>
#include "armadillo_state.hpp"
#include "system_and_poincare_surface.hpp"

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
        const auto& p = s[static_cast<unsigned>(CoordinateTag::p)];
        const auto& q = s[static_cast<unsigned>(CoordinateTag::q)];
        const auto& F = s[static_cast<unsigned>(CoordinateTag::F)];
        return M_ * p * p / 2 - F * cos(q);
      }

      double dp (const myState& s) const noexcept
      {
        const auto& p = s[static_cast<unsigned>(CoordinateTag::p)];
        return M_ * p;
      }

      double dq (const myState& s) const noexcept
      {
        const auto& q = s[static_cast<unsigned>(CoordinateTag::q)];
        const auto& F = s[static_cast<unsigned>(CoordinateTag::F)];
        return F * sin(q);
      }

      double dF (const myState& s) const noexcept
      {
        const auto& q = s[static_cast<unsigned>(CoordinateTag::q)];
        return -cos(q);
      }

    };

    template<typename UnperturbedHamiltonian>
    class UnperturbedDynamicSystem {

     public:
      using StateType = typename UnperturbedHamiltonian::StateType;

     private:
      UnperturbedHamiltonian h_;

      double dpdt(const StateType& s) const
      {
       return -h_.dq(s);
      };
      double dqdt(const StateType& s) const
      {
        return h_.dp(s);
      }
      double dphidt(const StateType& s) const
      {
        return h_.dF(s);
      }

      double oneFormTimeDerivative(const OneForm& of, const StateType& s) const
      {
        return of.p*dpdt(s) + of.q*dqdt(s);
      };

      OneForm dJ(const StateType& s) const
      {
        const double p = s[static_cast<unsigned>(CoordinateTag::p)];
        return OneForm{0,p};
      }

     public:

      explicit UnperturbedDynamicSystem (const UnperturbedHamiltonian& h)
          : h_{h}
      { }

      void operator() (const StateType& s, StateType& dsdt, const double /*t*/) const
      {

        dsdt[static_cast<unsigned>(CoordinateTag::p)] = dpdt(s);
        dsdt[static_cast<unsigned>(CoordinateTag::q)] = dqdt(s);
        dsdt[static_cast<unsigned>(CoordinateTag::F)] = 0;
        dsdt[static_cast<unsigned>(CoordinateTag::phi)] = dphidt(s);
        dsdt[static_cast<unsigned>(CoordinateTag::J)] = oneFormTimeDerivative(dJ(s),s);
        dsdt[static_cast<unsigned>(CoordinateTag::t)] = 1;
      }
    };

    template<typename UnperturbedHamiltonian>
    UnperturbedDynamicSystem<UnperturbedHamiltonian> makeUnperturbedDynamicSystem (const UnperturbedHamiltonian& h)
    {
      return UnperturbedDynamicSystem<UnperturbedHamiltonian>(h);
    }

}
#endif //ODE_INTEGRATORS_HAMILTONIAN_DYNAMIC_SYSTEM_HPP
