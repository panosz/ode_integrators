//
// Created by Panagiotis Zestanakis on 23/04/19.
//

#include "hamiltonian_dynamic_system.hpp"

namespace DS
{

    using UE_O = UnperturbedExtendedOscillatorHamiltonian;
    using UE_P  = UnperturbedExtendedPendulumHamiltonian;
    using P_S = PhaseSpaceState;
    using E_S = ExtendedSpaceState;


    /**********************************************
    *  UnperturbedExtendedOscillatorHamiltonian  *
    **********************************************/


      template<typename ST>
      double UE_O::value(const ST& s) const noexcept
      {
        const auto& p = s[static_cast<unsigned>(CoordinateTag::p)];
        const auto& q = s[static_cast<unsigned>(CoordinateTag::q)];
        const auto& F = s[static_cast<unsigned>(CoordinateTag::F)];
        return M_ * p * p  + F * q * q;
      }
      template double UE_O::value<P_S>(const P_S&) const noexcept;
      template double UE_O::value<E_S>(const E_S&) const noexcept;



      template<typename ST>
      double  UE_O::operator() (const ST& s) const noexcept
      {
        return value(s);
      }
      template double  UE_O::operator()(const P_S& ) const noexcept;
      template double  UE_O::operator()(const E_S& ) const noexcept;


      template<typename ST>
      FirstDerivatives UE_O::first_derivatives (const ST& s) const noexcept
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
      template FirstDerivatives UE_O::first_derivatives(const P_S& s) const noexcept;
      template FirstDerivatives UE_O::first_derivatives(const E_S& s) const noexcept;


      template<typename ST>
      SecondDerivatives UE_O::second_derivatives(const ST& s) const noexcept
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
      template SecondDerivatives UE_O::second_derivatives(const P_S&) const noexcept;
      template SecondDerivatives UE_O::second_derivatives(const E_S&) const noexcept;


      template<typename ST>
      ThirdDerivatives UE_O::third_derivatives(const ST& /*s*/) const noexcept
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
      template ThirdDerivatives UE_O::third_derivatives(const P_S&) const noexcept;
      template ThirdDerivatives UE_O::third_derivatives(const E_S&) const noexcept;



    double UE_O::action
      (const P_S& s) const
    {
      const auto& F = s[static_cast<unsigned>(CoordinateTag::F)];
      return value(s) / (2 * std::sqrt(M_ * F));
    }

    double UE_O::dKdJ
      (const P_S& s) const noexcept
    {
      const auto& F = s[static_cast<unsigned>(CoordinateTag::F)];

      return 2 * std::sqrt(M_ * F);
    }

    double UE_O::dKdF (const P_S& s) const
    {
      const auto& F = s[static_cast<unsigned>(CoordinateTag::F)];

      return action(s) * std::sqrt(M_ / F);
    }

    double UE_O::d2KdJ2 (const P_S&) const noexcept
    {
      return 0;
    }

    double UE_O::d2KdJdF (const P_S& s) const
    {
      const auto& F = s[static_cast<unsigned>(CoordinateTag::F)];

      return std::sqrt(M_ / F);
    }

    double UE_O::d2KdF2 (const P_S& s) const
    {
      const auto& F = s[static_cast<unsigned>(CoordinateTag::F)];

      const auto J = action(s);

      return -J * std::sqrt(M_ / F) / (2 * F);
    }

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
