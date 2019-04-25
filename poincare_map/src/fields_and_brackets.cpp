//
// Created by Panagiotis Zestanakis on 25/04/19.
//

#include <boost/math/special_functions/pow.hpp>
#include <boost/math/constants/constants.hpp>
#include "fields_and_brackets.hpp"
namespace DS
{

    Field::Field (double P, double Q)
        : p(P), q(Q)
    { }

    FieldAndFirstDerivatives::FieldAndFirstDerivatives (const DS::Field& F, const DS::FieldFirstDerivatives& dF)
        : f{F}, df{dF}
    { }

    OneFormAndFirstDerivatives::OneFormAndFirstDerivatives (const OneForm& G, const OneFormFirstDerivatives& dG)
        : g{G}, dg{dG}
    { }

    VelocitySqAndFirstDerivatives
    calculate_vSq_and_first_derivatives (const FirstDerivatives& dh, const SecondDerivatives& d2h)
    {
      using boost::math::pow;
      VelocitySqAndFirstDerivatives output{};

      output.v_Sq = pow<2>(dh.dp) + pow<2>(dh.dq);

      output.dv_Sq.dp = 2 * (dh.dp * d2h.dp2 + dh.dq * d2h.dp_dq);
      output.dv_Sq.dq = 2 * (dh.dp * d2h.dp_dq + dh.dq * d2h.dq2);
      output.dv_Sq.dF = 2 * (dh.dp * d2h.dp_dF + dh.dq * d2h.dq_dF);

      return output;
    }

    FieldAndFirstDerivatives
    caluclate_translation_field_and_first_derivatives (const FirstDerivatives& dh, const SecondDerivatives& d2h)
    {
      using boost::math::pow;
      Field f{};
      FieldFirstDerivatives df{};

      const auto[v_sq, dv_sq] = calculate_vSq_and_first_derivatives(dh, d2h);

      f.p = dh.dp / v_sq;
      f.q = dh.dq / v_sq;

      const auto pow_2_v_sq = pow<2>(v_sq);

      df.p.dp = (v_sq * d2h.dp2 - dh.dp * dv_sq.dp) / pow_2_v_sq;
      df.p.dq = (v_sq * d2h.dp_dq - dh.dp * dv_sq.dq) / pow_2_v_sq;
      df.p.dF = (v_sq * d2h.dp_dF - dh.dp * dv_sq.dF) / pow_2_v_sq;

      df.q.dp = (v_sq * d2h.dp_dq - dh.dq * dv_sq.dp) / pow_2_v_sq;
      df.q.dq = (v_sq * d2h.dq2 - dh.dq * dv_sq.dq) / pow_2_v_sq;
      df.q.dF = (v_sq * d2h.dq_dF - dh.dq * dv_sq.dF) / pow_2_v_sq;

      return FieldAndFirstDerivatives{f, df};
    }

    OneForm calculate_beta (double p, const FieldAndFirstDerivatives& fieldAndFirstDerivatives)
    {
      OneForm beta{};
      using boost::math::double_constants::two_pi;

      const auto &[f, df] = fieldAndFirstDerivatives;


      //OneForm[2 p \[Pi] Dt[fq, p], 2 \[Pi] (fp + p Dt[fq, q])]

      beta.p = two_pi * p * df.q.dp;
      beta.q = two_pi * (f.p + p * df.q.dq);

      return beta;
    }
    OneForm
    calculate_gamma (double p,
                     const FieldAndFirstDerivatives& fieldAndFirstDerivatives,
                     const FirstDerivatives& dh,
                     const SecondDerivatives& d2h)
    {
      OneForm gamma{};

      const auto &[f, df] = fieldAndFirstDerivatives;


      //OneForm[Hf p Dt[fq, p] + fq p Dt[Hf, p], fp Hf + Hf p Dt[fq, q] + fq p Dt[Hf, q]]
      gamma.p = dh.dF * p * df.q.dp + f.q * p * d2h.dp_dF;
      gamma.q = f.p * dh.dF + dh.dF * p * df.q.dq + f.q * p * d2h.dq_dF;
      return gamma;

    }

}