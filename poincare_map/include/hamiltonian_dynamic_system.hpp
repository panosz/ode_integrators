//
// Created by Panagiotis Zestanakis on 23/04/19.
//

#ifndef ODE_INTEGRATORS_HAMILTONIAN_DYNAMIC_SYSTEM_HPP
#define ODE_INTEGRATORS_HAMILTONIAN_DYNAMIC_SYSTEM_HPP
#include <cmath>
#include "armadillo_state.hpp"
#include "system_and_poincare_surface.hpp"
#include "fields_and_brackets.hpp"
#include <boost/math/constants/constants.hpp>

namespace DS
{

    using myState = armadillo_state<8>;

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

      FirstDerivatives first_derivatives (const myState& s) const noexcept
      {
        FirstDerivatives derivs{};
        const auto& p = s[static_cast<unsigned>(CoordinateTag::p)];
        const auto& q = s[static_cast<unsigned>(CoordinateTag::q)];
        const auto& F = s[static_cast<unsigned>(CoordinateTag::F)];

        derivs.dp = M_ * p;
        derivs.dq = F * sin(q);
        derivs.dF = -cos(q);
        return derivs;
      }

      SecondDerivatives second_derivatives (const myState& s) const noexcept
      {
        SecondDerivatives second_derivs{};

        const auto& p = s[static_cast<unsigned>(CoordinateTag::p)];
        const auto& q = s[static_cast<unsigned>(CoordinateTag::q)];
        const auto& F = s[static_cast<unsigned>(CoordinateTag::F)];

        second_derivs.dp2 = M_;
        second_derivs.dp_dq = 0;
        second_derivs.dp_dF = 0;
        second_derivs.dq2 = F * cos(q);
        second_derivs.dq_dF = sin(q);
        second_derivs.dF2 = 0;

        return second_derivs;
      }

    };

    template<typename UnperturbedHamiltonian>
    class UnperturbedDynamicSystem {

     public:
      using StateType = typename UnperturbedHamiltonian::StateType;

     private:
      UnperturbedHamiltonian h_;

      double dpdt (const FirstDerivatives& dh) const
      {
        return -dh.dq;
      };
      double dqdt (const FirstDerivatives& dh) const
      {
        return dh.dp;
      }
      double dphidt (const FirstDerivatives& dh) const
      {
        return dh.dF;
      }

      double oneFormTimeDerivative (const OneForm& of, const FirstDerivatives& dh) const
      {
        return of.p * dpdt(dh) + of.q * dqdt(dh);
      };

      OneForm dJ (const StateType& s) const
      {
        const double p = s[static_cast<unsigned>(CoordinateTag::p)];
        return OneForm{0, p};
      }

     public:

      explicit UnperturbedDynamicSystem (const UnperturbedHamiltonian& h)
          : h_{h}
      { }

      void operator() (const StateType& s, StateType& dsdt, const double /*t*/) const
      {

        const auto p = s[static_cast<unsigned >(CoordinateTag::p)];
        const auto dh = h_.first_derivatives(s);
        const auto d2h =h_.second_derivatives(s);

        const auto f_and_df =caluclate_translation_field_and_first_derivatives(dh,d2h);
        const auto beta = calculate_beta(p,f_and_df);
        const auto gamma = calculate_gamma(p,f_and_df,dh,d2h);


        dsdt[static_cast<unsigned>(CoordinateTag::p)] = dpdt(dh);
        dsdt[static_cast<unsigned>(CoordinateTag::q)] = dqdt(dh);
        dsdt[static_cast<unsigned>(CoordinateTag::F)] = 0;
        dsdt[static_cast<unsigned>(CoordinateTag::phi)] = dphidt(dh);
        dsdt[static_cast<unsigned>(CoordinateTag::J)] = oneFormTimeDerivative(dJ(s), dh);
        dsdt[static_cast<unsigned>(CoordinateTag::t)] = 1;
        dsdt[static_cast<unsigned>(CoordinateTag::beta)] = oneFormTimeDerivative(beta,dh);
        dsdt[static_cast<unsigned>(CoordinateTag::gamma)] = oneFormTimeDerivative(gamma,dh);

      }
    };

    template<typename UnperturbedHamiltonian>
    UnperturbedDynamicSystem<UnperturbedHamiltonian> makeUnperturbedDynamicSystem (const UnperturbedHamiltonian& h)
    {
      return UnperturbedDynamicSystem<UnperturbedHamiltonian>(h);
    }
}
#endif //ODE_INTEGRATORS_HAMILTONIAN_DYNAMIC_SYSTEM_HPP
