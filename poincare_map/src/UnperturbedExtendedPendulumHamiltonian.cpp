//
// Created by Panagiotis Zestanakis on Nov 01, 2019
//

#include "hamiltonian_dynamic_system.hpp"
#include "UnperturbedExtendedPendulumHamiltonian.hpp"

namespace DS
{

    using UE_P  = UnperturbedExtendedPendulumHamiltonian;
    using P_S = PhaseSpaceState;
    using E_S = ExtendedSpaceState;

    /********************************************
    *  UnperturbedExtendedPendulumHamiltonian  *
    ********************************************/

      template<typename ST>
      double UE_P::value(const ST& s) const noexcept
      {
        const auto& p = s[static_cast<unsigned>(CoordinateTag::p)];
        const auto& q = s[static_cast<unsigned>(CoordinateTag::q)];
        const auto& F = s[static_cast<unsigned>(CoordinateTag::F)];
        return M_ * p * p / 2 - F * cos(q);
      }
      template double UE_P::value<P_S>(const P_S& s) const noexcept;
      template double UE_P::value<E_S>(const E_S& s) const noexcept;

      template<typename ST>
      double UE_P::operator() (const ST& s) const
      {
        return value(s);
      }
      template double UE_P::operator() (const P_S& s) const;
      template double UE_P::operator() (const E_S& s) const;

      template<typename ST>
      FirstDerivatives UE_P::first_derivatives (const ST& s) const noexcept
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
      template FirstDerivatives UE_P::first_derivatives(const P_S& s) const noexcept;
      template FirstDerivatives UE_P::first_derivatives(const E_S& s) const noexcept;

      template<typename ST>
      SecondDerivatives UE_P::second_derivatives (const ST& s) const noexcept
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
      template SecondDerivatives UE_P::second_derivatives(const P_S& s) const noexcept;
      template SecondDerivatives UE_P::second_derivatives(const E_S& s) const noexcept;


      template<typename ST>
      ThirdDerivatives UE_P::third_derivatives (const ST& s) const noexcept
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

      template ThirdDerivatives UE_P::third_derivatives(const P_S& s) const noexcept;
      template ThirdDerivatives UE_P::third_derivatives(const E_S& s) const noexcept;
}
