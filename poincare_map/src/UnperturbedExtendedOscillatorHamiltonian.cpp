//
//Created by Panagiotis Zestanakis on Nov 01, 2019
//

#include "UnperturbedExtendedOscillatorHamiltonian.hpp"

namespace DS
{
    using UE_O = UnperturbedExtendedOscillatorHamiltonian;
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

    PhaseSpaceState UE_O::propagate(const PhaseSpaceState& s,
                                    double dt) const

    {
      PhaseSpaceState propagated{};
      using namespace std;

      const auto& p0 = s[static_cast<unsigned>(CoordinateTag::p)];
      const auto& q0 = s[static_cast<unsigned>(CoordinateTag::q)];
      const auto& F = s[static_cast<unsigned>(CoordinateTag::F)];
      const auto& phi0 = s[static_cast<unsigned>(CoordinateTag::phi)];

      const double E = value(s);
      const double A = sqrt(E);
      const double sqrtF = sqrt(F);
      const double sqrtM = sqrt(M_);

      const double theta_0 = atan2(sqrtF * q0,
                                        sqrtM * p0);

      const double omega = dKdJ(s);
      const double omega_phi = dKdF(s);
      const double propagation_phase =  omega * dt + theta_0;

      const double p = A / sqrtM * cos(propagation_phase);
      const double q = A / sqrtF * sin(propagation_phase);

      const double phi = omega_phi * dt
                         - 0.5 * omega_phi/omega * sin(2 * (propagation_phase))
                         + phi0;

      propagated[static_cast<unsigned>(CoordinateTag::p)] = p;
      propagated[static_cast<unsigned>(CoordinateTag::q)] = q;
      propagated[static_cast<unsigned>(CoordinateTag::F)] = F;
      propagated[static_cast<unsigned>(CoordinateTag::phi)] = phi;

      return propagated;

    }




}
