//
// Created by Panagiotis Zestanakis on 23/04/19.
//

#ifndef ODE_INTEGRATORS_HAMILTONIAN_DYNAMIC_SYSTEM_HPP
#define ODE_INTEGRATORS_HAMILTONIAN_DYNAMIC_SYSTEM_HPP
#include <cmath>
#include "fields_and_brackets.hpp"
#include "phase_space_state_types.hpp"

namespace DS
{



    inline double dpdt (const FirstDerivatives& dh) noexcept
    {
      return -dh.dq;
    }
    inline double dqdt (const FirstDerivatives& dh) noexcept
    {
      return dh.dp;
    }
    inline double dphidt (const FirstDerivatives& dh) noexcept
    {
      return dh.dF;
    }

    inline void get_time_derivatives(const FirstDerivatives& dh,
                                    PhaseSpaceState& dsdt) noexcept
    {

      dsdt[static_cast<unsigned>(CoordinateTag::p)] = dpdt(dh);
      dsdt[static_cast<unsigned>(CoordinateTag::q)] = dqdt(dh);
      dsdt[static_cast<unsigned>(CoordinateTag::F)] = 0;
      dsdt[static_cast<unsigned>(CoordinateTag::phi)] = dphidt(dh);
    }


    template<typename UnperturbedHamiltonian>
    class UnperturbedDynamicSystem {

     private:
      UnperturbedHamiltonian h_;

     public:
      using StateType = PhaseSpaceState;

      explicit UnperturbedDynamicSystem (const UnperturbedHamiltonian& h)
          : h_{h}
      { }

      void operator() (const PhaseSpaceState& s,
                       PhaseSpaceState& dsdt,
                       const double /*t*/) const
      {
        const auto dh = h_.first_derivatives(s);
        get_time_derivatives(dh, dsdt);
      }
    };

    template<typename UnperturbedHamiltonian>
    UnperturbedDynamicSystem<UnperturbedHamiltonian> makeUnperturbedDynamicSystem (const UnperturbedHamiltonian& h)
    {
      return UnperturbedDynamicSystem<UnperturbedHamiltonian>(h);
    }

    template<typename UnperturbedHamiltonian>
    class ActionDynamicSystem {

     public:
      using StateType = ExtendedSpaceState;

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

      explicit ActionDynamicSystem (const UnperturbedHamiltonian& h)
          : h_{h}
      { }

      void operator() (const StateType& s, StateType& dsdt, const double /*t*/) const
      {

        const auto p = s[static_cast<unsigned >(CoordinateTag::p)];
        const auto dh = h_.first_derivatives(s);
        const auto d2h = h_.second_derivatives(s);
        const auto d3h = h_.third_derivatives(s);

        const auto f_and_df = caluclate_translation_field_and_derivatives(dh, d2h, d3h);
        const auto beta_and_dbeta = calculate_beta(p, f_and_df);
        const auto& beta = beta_and_dbeta.g;
        const auto gamma_and_dgamma = calculate_gamma(p, f_and_df, dh, d2h, d3h);
        const auto& gamma = gamma_and_dgamma.g;
        const auto beta1 = calculate_beta1(f_and_df, beta_and_dbeta);
        const auto gamma1 = calculate_gamma1(f_and_df, gamma_and_dgamma);
        const auto beta2 = calculate_beta2(f_and_df, beta_and_dbeta, dh, d2h);
        const auto gamma2 = calculate_gamma2(f_and_df, gamma_and_dgamma, dh, d2h);

        dsdt[static_cast<unsigned>(CoordinateTag::p)] = dpdt(dh);
        dsdt[static_cast<unsigned>(CoordinateTag::q)] = dqdt(dh);
        dsdt[static_cast<unsigned>(CoordinateTag::F)] = 0;
        dsdt[static_cast<unsigned>(CoordinateTag::phi)] = dphidt(dh);
        dsdt[static_cast<unsigned>(CoordinateTag::J)] = oneFormTimeDerivative(dJ(s), dh);
        dsdt[static_cast<unsigned>(CoordinateTag::t)] = 1;
        dsdt[static_cast<unsigned>(CoordinateTag::beta)] = oneFormTimeDerivative(beta, dh);
        dsdt[static_cast<unsigned>(CoordinateTag::gamma)] = oneFormTimeDerivative(gamma, dh);
        dsdt[static_cast<unsigned>(CoordinateTag::beta1)] = oneFormTimeDerivative(beta1, dh);
        dsdt[static_cast<unsigned>(CoordinateTag::gamma1)] = oneFormTimeDerivative(gamma1, dh);
        dsdt[static_cast<unsigned>(CoordinateTag::beta2)] = oneFormTimeDerivative(beta2, dh);
        dsdt[static_cast<unsigned>(CoordinateTag::gamma2)] = oneFormTimeDerivative(gamma2, dh);

      }
    };

    template<typename UnperturbedHamiltonian>
    ActionDynamicSystem<UnperturbedHamiltonian> makeActionDynamicSystem (const UnperturbedHamiltonian& h)
    {
      return ActionDynamicSystem<UnperturbedHamiltonian>(h);
    }
}
#endif //ODE_INTEGRATORS_HAMILTONIAN_DYNAMIC_SYSTEM_HPP
