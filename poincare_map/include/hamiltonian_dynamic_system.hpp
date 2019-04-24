//
// Created by Panagiotis Zestanakis on 23/04/19.
//

#ifndef ODE_INTEGRATORS_HAMILTONIAN_DYNAMIC_SYSTEM_HPP
#define ODE_INTEGRATORS_HAMILTONIAN_DYNAMIC_SYSTEM_HPP
#include <cmath>
#include "armadillo_state.hpp"
#include "system_and_poincare_surface.hpp"
#include <boost/math/special_functions/pow.hpp>
namespace DS
{

    using myState = armadillo_state<6>;

    struct OneForm {
        double p = 0;
        double q = 0;
        OneForm () = default;

        OneForm (double P, double Q)
            : p(P), q(Q)
        { }
    };

    struct Field {
        double p = 0;
        double q = 0;

        Field () = default;
        Field (double p, double q);
    };

    struct FirstDerivatives {
        double dp = 0;
        double dq = 0;
        double dF = 0;
    };

    struct SecondDerivatives {
        double dp2 = 0;
        double dp_dq = 0;
        double dp_dF = 0;

        double dq2 = 0;
        double dq_dF = 0;

        double dF2 = 0;
    };

    struct ThirdDerivatives {

        double dp3 = 0;
        double dp2_dq = 0;
        double dp2_dF = 0;
        double dp_dq2 = 0;
        double dp_dq_dF = 0;
        double dp_dF2 = 0;

        double dq3 = 0;
        double dq2_dF = 0;
        double dq_dF2 = 0;

        double dF3 = 0;

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

        const auto dh = h_.first_derivatives(s);

        dsdt[static_cast<unsigned>(CoordinateTag::p)] = dpdt(dh);
        dsdt[static_cast<unsigned>(CoordinateTag::q)] = dqdt(dh);
        dsdt[static_cast<unsigned>(CoordinateTag::F)] = 0;
        dsdt[static_cast<unsigned>(CoordinateTag::phi)] = dphidt(dh);
        dsdt[static_cast<unsigned>(CoordinateTag::J)] = oneFormTimeDerivative(dJ(s), dh);
        dsdt[static_cast<unsigned>(CoordinateTag::t)] = 1;
      }
    };

    template<typename UnperturbedHamiltonian>
    UnperturbedDynamicSystem<UnperturbedHamiltonian> makeUnperturbedDynamicSystem (const UnperturbedHamiltonian& h)
    {
      return UnperturbedDynamicSystem<UnperturbedHamiltonian>(h);
    }

    struct VelocitySqAndFirstDerivatives {
        double v_Sq;
        FirstDerivatives dv_Sq;
    };

    VelocitySqAndFirstDerivatives
    v_fromHamiltonianDerivatives (const FirstDerivatives& dh, const SecondDerivatives& d2h)
    {
      using boost::math::pow;
      VelocitySqAndFirstDerivatives output{};

      output.v_Sq = pow<2>(dh.dp) + pow<2>(dh.dq);

      output.dv_Sq.dp = 2 * (dh.dp * d2h.dp2 + dh.dq * d2h.dp_dq);
      output.dv_Sq.dq = 2 * (dh.dp * d2h.dp_dq + dh.dq * d2h.dq2);
      output.dv_Sq.dF = 2 * (dh.dp * d2h.dp_dF + dh.dq * d2h.dq_dF);

      return output;
    };



    struct FieldFirstDerivatives {
        FirstDerivatives p{};
        FirstDerivatives q{};
    };

    struct FieldAndFirstDerivatives {
        Field f{};
        FieldFirstDerivatives df{};
        FieldAndFirstDerivatives (const Field& F, const FieldFirstDerivatives& dF);
    };

    FieldAndFirstDerivatives f_fromHamiltonianDerivatives (const FirstDerivatives& dh, const SecondDerivatives& d2h)
    {
      using boost::math::pow;
      Field f{};
      FieldFirstDerivatives df{};

      const auto [v_sq, dv_sq] = v_fromHamiltonianDerivatives(dh,d2h);

      f.p = dh.dp / v_sq;
      f.q = dh.dq / v_sq;

      const auto pow_2_v_sq = pow<2>(v_sq);


      df.p.dp = (v_sq * d2h.dp2 - dh.dp * dv_sq.dp)/pow_2_v_sq;
      df.p.dq = (v_sq * d2h.dp_dq - dh.dp * dv_sq.dq)/pow_2_v_sq;
      df.p.dF = (v_sq * d2h.dp_dF - dh.dp * dv_sq.dF)/pow_2_v_sq;

      df.q.dp = (v_sq * d2h.dp_dq - dh.dq * dv_sq.dp)/pow_2_v_sq;
      df.q.dq = (v_sq * d2h.dq2 - dh.dq * dv_sq.dq)/pow_2_v_sq;
      df.q.dF = (v_sq * d2h.dq_dF - dh.dq * dv_sq.dF)/pow_2_v_sq;

      return FieldAndFirstDerivatives{f,df};
    }
}
#endif //ODE_INTEGRATORS_HAMILTONIAN_DYNAMIC_SYSTEM_HPP
