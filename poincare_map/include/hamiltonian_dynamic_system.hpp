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

    using PhaseSpaceState = armadillo_state<PHASE_SPACE_VARIABLES>;
    using ExtendedSpaceState = armadillo_state<EXTENDED_SPACE_VARIABLES>;

    ExtendedSpaceState phase_to_extended_space_state(const PhaseSpaceState &pss);

    class UnperturbedExtendedOscillatorHamiltonian {

     private:
      double M_;

     public:
      using StateType=ExtendedSpaceState;

      explicit UnperturbedExtendedOscillatorHamiltonian (double M)
          : M_{M}
      { };

      template<typename ST>
      double value(const ST& s) const noexcept
      {
        const auto& p = s[static_cast<unsigned>(CoordinateTag::p)];
        const auto& q = s[static_cast<unsigned>(CoordinateTag::q)];
        const auto& F = s[static_cast<unsigned>(CoordinateTag::F)];
        return M_ * p * p  + F * q * q;
      }

      template<typename ST>
      double operator() (const ST& s) const noexcept
      {
        return value(s);
      }

      template<typename ST>
      FirstDerivatives first_derivatives (const ST& s) const noexcept
      {
        FirstDerivatives derivs{};
        const auto& p = s[static_cast<unsigned>(CoordinateTag::p)];
        const auto& q = s[static_cast<unsigned>(CoordinateTag::q)];
        const auto& F = s[static_cast<unsigned>(CoordinateTag::F)];

        derivs.dp = 2 * M_ * p;
        derivs.dq = 2 * F * q;
        derivs.dF =  q * q;
        return derivs;
      }

      template<typename ST>
      SecondDerivatives second_derivatives (const ST& s) const noexcept
      {
        SecondDerivatives second_derivs{};

        const auto& q = s[static_cast<unsigned>(CoordinateTag::q)];
        const auto& F = s[static_cast<unsigned>(CoordinateTag::F)];

        second_derivs.dp2 = 2 * M_;
        second_derivs.dp_dq = 0;
        second_derivs.dp_dF = 0;
        second_derivs.dq2 = 2 * F;
        second_derivs.dq_dF = 2 * q;
        second_derivs.dF2 = 0;

        return second_derivs;
      }

      template<typename ST>
      ThirdDerivatives third_derivatives (const ST& /*s*/) const noexcept
      {
        ThirdDerivatives third_derivs{};


        third_derivs.dp3 = 0;
        third_derivs.dp2_dq = 0;
        third_derivs.dp2_dF = 0;
        third_derivs.dp_dq2 = 0;
        third_derivs.dp_dq_dF = 0;
        third_derivs.dp_dF2 = 0;

        third_derivs.dq3 = 0;
        third_derivs.dq2_dF = 2;
        third_derivs.dq_dF2 = 0;

        third_derivs.dF3 = 0;

        return third_derivs;
      }

      double action(const PhaseSpaceState& s) const;

      double dKdJ(const PhaseSpaceState& s) const noexcept;

      double dKdF(const PhaseSpaceState& s) const;

      double d2KdJ2(const PhaseSpaceState& /*s*/) const noexcept;

      double d2KdJdF(const PhaseSpaceState& s) const;

      double d2KdF2(const PhaseSpaceState& s) const;


    };

    class UnperturbedExtendedPendulumHamiltonian {

     private:
      double M_;

     public:
      using StateType=PhaseSpaceState;

      explicit UnperturbedExtendedPendulumHamiltonian (double M)
          : M_{M}
      { };

      double operator() (const StateType& s) const
      {
        const auto& p = s[static_cast<unsigned>(CoordinateTag::p)];
        const auto& q = s[static_cast<unsigned>(CoordinateTag::q)];
        const auto& F = s[static_cast<unsigned>(CoordinateTag::F)];
        return M_ * p * p / 2 - F * cos(q);
      }

      FirstDerivatives first_derivatives (const StateType& s) const noexcept
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

      SecondDerivatives second_derivatives (const StateType& s) const noexcept
      {
        SecondDerivatives second_derivs{};

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

      ThirdDerivatives third_derivatives (const StateType& s) const noexcept
      {
        ThirdDerivatives third_derivs{};

        const auto& q = s[static_cast<unsigned>(CoordinateTag::q)];
        const auto& F = s[static_cast<unsigned>(CoordinateTag::F)];

        third_derivs.dp3 = 0;
        third_derivs.dp2_dq = 0;
        third_derivs.dp2_dF = 0;
        third_derivs.dp_dq2 = 0;
        third_derivs.dp_dq_dF = 0;
        third_derivs.dp_dF2 = 0;

        third_derivs.dq3 = -F * sin(q);
        third_derivs.dq2_dF = cos(q);
        third_derivs.dq_dF2 = 0;

        third_derivs.dF3 = 0;

        return third_derivs;
      }

    };

    template<typename UnperturbedHamiltonian>
    class UnperturbedDynamicSystem {

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

      explicit UnperturbedDynamicSystem (const UnperturbedHamiltonian& h)
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
    UnperturbedDynamicSystem<UnperturbedHamiltonian> makeUnperturbedDynamicSystem (const UnperturbedHamiltonian& h)
    {
      return UnperturbedDynamicSystem<UnperturbedHamiltonian>(h);
    }
}
#endif //ODE_INTEGRATORS_HAMILTONIAN_DYNAMIC_SYSTEM_HPP
