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

    VelocitySqAndDerivatives
    calculate_vSq_and_derivatives (const FirstDerivatives& dh, const SecondDerivatives& d2h, const ThirdDerivatives& d3h)
    {
      using boost::math::pow;
      VelocitySqAndDerivatives output{};

      output.v_Sq = pow<2>(dh.dp) + pow<2>(dh.dq);

      output.dv_Sq.dp = 2 * (dh.dp * d2h.dp2 + dh.dq * d2h.dp_dq);
      output.dv_Sq.dq = 2 * (dh.dp * d2h.dp_dq + dh.dq * d2h.dq2);
      output.dv_Sq.dF = 2 * (dh.dp * d2h.dp_dF + dh.dq * d2h.dq_dF);


      //2 (Dt[hp, p]^2 + hp Dt[hp, {p, 2}] + Dt[hq, p]^2 + hq Dt[hq, {p, 2}])
      output.d2v_Sq.dp2 =
          2 * (pow<2>(d2h.dp2) + dh.dp * d3h.dp3 + pow<2>(d2h.dp_dq) + dh.dq * d3h.dp2_dq);

      //2 (Dt[hp, p] Dt[hp, q] + Dt[hq, p] Dt[hq, q] + hp Dt[hp, p, q] + hq Dt[hq, p, q])
      output.d2v_Sq.dp_dq =
          2 * (d2h.dp2 * d2h.dp_dq + d2h.dp_dq * d2h.dq2 + dh.dp * d3h.dp2_dq + dh.dq * d3h.dp_dq2);

      //2 (Dt[hp, F] Dt[hp, p]
      // + Dt[hq, F] Dt[hq, p]
      // + hp Dt[hp, F, p]
      // + hq Dt[hq, F, p])
      output.d2v_Sq.dp_dF =
          2 * (d2h.dp_dF * d2h.dp2
               + d2h.dq_dF * d2h.dp_dq
               + dh.dp * d3h.dp2_dF
               + dh.dq * d3h.dp_dq_dF);


      //2 (Dt[hp, q]^2
      // + hp Dt[hp, {q, 2}]
      // + Dt[hq, q]^2
      // + hq Dt[hq, {q, 2}])
      output.d2v_Sq.dq2 =
          2 * (pow<2>(d2h.dp_dq)
               + dh.dp * d3h.dp_dq2
               + pow<2>(d2h.dq2)
               + dh.dq * d3h.dq3);

      //2 (Dt[hp, F] Dt[hp, q]
      // + Dt[hq, F] Dt[hq, q]
      // + hp Dt[hp, F, q]
      //  + hq Dt[hq, F, q])
      output.d2v_Sq.dq_dF =
          2 * (d2h.dp_dF * d2h.dp_dq
               + d2h.dq_dF * d2h.dq2
               + dh.dp * d3h.dp_dq_dF
               + dh.dq * d3h.dq2_dF);

      //2 (Dt[hp, F]^2
      // + hp Dt[hp, {F, 2}]
      // + Dt[hq, F]^2
      // + hq Dt[hq, {F, 2}])
      output.d2v_Sq.dF2 =
          2 * (pow<2>(d2h.dp_dF)
               + dh.dp * d3h.dp_dF2
               + pow<2>(d2h.dq_dF)
               + dh.dq * d3h.dq_dF2);

      return output;
    }

    FieldAndDerivatives
    caluclate_translation_field_and_derivatives (const FirstDerivatives& dh,
                                                 const SecondDerivatives& d2h,
                                                 const ThirdDerivatives& d3h)
    {
      using boost::math::pow;
      FieldAndDerivatives output{};

      auto &[f, df, d2f] = output;

      const auto[v_sq, dv_sq, d2v_sq] = calculate_vSq_and_derivatives(dh, d2h, d3h);

      f.p = dh.dp / v_sq;
      f.q = dh.dq / v_sq;

      const auto pow_2_v_sq = pow<2>(v_sq);

      // first derivatives

      //(vsq Dt[hp, p] - hp Dt[vsq, p])/vsq^2
      df.p.dp = (v_sq * d2h.dp2 - dh.dp * dv_sq.dp) / pow_2_v_sq;

      //(vsq Dt[hp, q] - hp Dt[vsq, q])/vsq^2
      df.p.dq = (v_sq * d2h.dp_dq - dh.dp * dv_sq.dq) / pow_2_v_sq;

      //(vsq Dt[hp, F] - hp Dt[vsq, F])/vsq^2
      df.p.dF = (v_sq * d2h.dp_dF - dh.dp * dv_sq.dF) / pow_2_v_sq;

      //(vsq Dt[hq, p] - hq Dt[vsq, p])/vsq^2
      df.q.dp = (v_sq * d2h.dp_dq - dh.dq * dv_sq.dp) / pow_2_v_sq;

      //(vsq Dt[hq, q] - hq Dt[vsq, q])/vsq^2
      df.q.dq = (v_sq * d2h.dq2 - dh.dq * dv_sq.dq) / pow_2_v_sq;

      //(vsq Dt[hq, F] - hq Dt[vsq, F])/vsq^2
      df.q.dF = (v_sq * d2h.dq_dF - dh.dq * dv_sq.dF) / pow_2_v_sq;

      const auto pow_3_v_sq = pow<3>(v_sq);

      //second derivatives


      //(vsq^2 Dt[hp, {p, 2}]
      // - 2 vsq Dt[hp, p] Dt[vsq, p]
      // + 2 hp Dt[vsq, p]^2
      // - hp vsq Dt[vsq, {p, 2}])/vsq^3

      d2f.p.dp2 = (pow_2_v_sq * d3h.dp3
                   - 2 * v_sq * d2h.dp2 * dv_sq.dp
                   + 2 * dh.dp * pow<2>(dv_sq.dp)
                   - dh.dp * v_sq * d2v_sq.dp2) / pow_3_v_sq;


      //-(1/(vsq^3))(vsq Dt[hp, q] Dt[vsq, p]
      // + vsq Dt[hp, p] Dt[vsq, q]
      // - 2 hp Dt[vsq, p] Dt[vsq, q]
      // - vsq^2 Dt[hp, p, q]
      // + hp vsq Dt[vsq, p, q])

      d2f.p.dp_dq = -(v_sq * d2h.dp_dq * dv_sq.dp
                      + v_sq * d2h.dp2 * dv_sq.dq
                      - 2 * dh.dp * dv_sq.dp * dv_sq.dq
                      - pow_2_v_sq * d3h.dp2_dq
                      + dh.dp * v_sq * d2v_sq.dp_dq) / pow_3_v_sq;

      //-(1/(vsq^3))(vsq Dt[hp, p] Dt[vsq, F]
      // + vsq Dt[hp, F] Dt[vsq, p]
      // -  2 hp Dt[vsq, F] Dt[vsq, p]
      // - vsq^2 Dt[hp, F, p]
      // +  hp vsq Dt[vsq, F, p])

      d2f.p.dp_dF = -(v_sq * d2h.dp2 * dv_sq.dF
                      + v_sq * d2h.dp_dF * dv_sq.dp
                      - 2 * dh.dp * dv_sq.dF * dv_sq.dp
                      - pow_2_v_sq * d3h.dp2_dF
                      + dh.dp * v_sq * d2v_sq.dp_dF) / pow_3_v_sq;


      //(vsq^2 Dt[hp, {q, 2}]
      // - 2 vsq Dt[hp, q] Dt[vsq, q]
      // + 2 hp Dt[vsq, q]^2
      // - hp vsq Dt[vsq, {q, 2}])/vsq^3
      d2f.p.dq2 = (pow_2_v_sq * d3h.dp_dq2
                   - 2 * v_sq * d2h.dp_dq * dv_sq.dq
                   + 2 * dh.dp * pow<2>(dv_sq.dq)
                   - dh.dp * v_sq * d2v_sq.dq2) / pow_3_v_sq;



      //-(1/(vsq^3))(vsq Dt[hp, q] Dt[vsq, F]
      // + vsq Dt[hp, F] Dt[vsq, q]
      // - 2 hp Dt[vsq, F] Dt[vsq, q]
      // - vsq^2 Dt[hp, F, q]
      // + hp vsq Dt[vsq, F, q])
      d2f.p.dq_dF = -(v_sq * d2h.dp_dq * dv_sq.dF
                      + v_sq * d2h.dp_dF * dv_sq.dq
                      - 2 * dh.dp * dv_sq.dF * dv_sq.dq
                      - pow_2_v_sq * d3h.dp_dq_dF
                      + dh.dp * v_sq * d2v_sq.dq_dF) / pow_3_v_sq;


      //(vsq^2 Dt[hp, {F, 2}]
      // - 2 vsq Dt[hp, F] Dt[vsq, F]
      // + 2 hp Dt[vsq, F]^2
      // - hp vsq Dt[vsq, {F, 2}])/vsq^3

      d2f.p.dF2 = (pow_2_v_sq * d3h.dp_dF2
                   - 2 * v_sq * d2h.dp_dF * dv_sq.dF
                   + 2 * dh.dp * pow<2>(dv_sq.dF)
                   - dh.dp * v_sq * d2v_sq.dF2) / pow_3_v_sq;



      //(vsq^2 Dt[hq, {p, 2}]
      // - 2 vsq Dt[hq, p] Dt[vsq, p]
      // + 2 hq Dt[vsq, p]^2
      // - hq vsq Dt[vsq, {p, 2}])/vsq^3

      d2f.q.dp2 = (pow_2_v_sq * d3h.dp2_dq
                   - 2 * v_sq * d2h.dp_dq * dv_sq.dp
                   + 2 * dh.dq * pow<2>(dv_sq.dp)
                   - dh.dq * v_sq * d2v_sq.dp2) / pow_3_v_sq;


      //-(1/(vsq^3))(vsq Dt[hq, q] Dt[vsq, p]
      // + vsq Dt[hq, p] Dt[vsq, q]
      // - 2 hq Dt[vsq, p] Dt[vsq, q]
      // - vsq^2 Dt[hq, p, q]
      // + hq vsq Dt[vsq, p, q])

      d2f.q.dp_dq = -(v_sq * d2h.dq2 * dv_sq.dp
                      + v_sq * d2h.dp_dq * dv_sq.dq
                      - 2 * dh.dq * dv_sq.dp * dv_sq.dq
                      - pow_2_v_sq * d3h.dp_dq2
                      + dh.dq * v_sq * d2v_sq.dp_dq) / pow_3_v_sq;


      //-(1/(vsq^3)) (vsq Dt[hq, p] Dt[vsq, F]
      // + vsq Dt[hq, F] Dt[vsq, p]
      // - 2 hq Dt[vsq, F] Dt[vsq, p]
      // - vsq^2 Dt[hq, F, p]
      // +  hq vsq Dt[vsq, F, p])

      d2f.q.dp_dF = -(v_sq * d2h.dp_dq * dv_sq.dF
                      + v_sq * d2h.dq_dF * dv_sq.dp
                      - 2 * dh.dq * dv_sq.dF * dv_sq.dp
                      - pow_2_v_sq * d3h.dp_dq_dF
                      + dh.dq * v_sq * d2v_sq.dp_dF) / pow_3_v_sq;


      //(vsq^2 Dt[hq, {q, 2}]
      // - 2 vsq Dt[hq, q] Dt[vsq, q]
      // + 2 hq Dt[vsq, q]^2
      // - hq vsq Dt[vsq, {q, 2}])/vsq^3

      d2f.q.dq2 = (pow_2_v_sq * d3h.dq3
                   - 2 * v_sq * d2h.dq2 * dv_sq.dq
                   + 2 * dh.dq * pow<2>(dv_sq.dq)
                   - dh.dq * v_sq * d2v_sq.dq2) / pow_3_v_sq;


      //-(1/(vsq^3))(vsq Dt[hq, q] Dt[vsq, F]
      //  + vsq Dt[hq, F] Dt[vsq, q]
      //  - 2 hq Dt[vsq, F] Dt[vsq, q]
      //  - vsq^2 Dt[hq, F, q]
      //  +  hq vsq Dt[vsq, F, q])

      d2f.q.dq_dF = -(v_sq * d2h.dq2 * dv_sq.dF
                      + v_sq * d2h.dq_dF * dv_sq.dq
                      - 2 * dh.dq * dv_sq.dF * dv_sq.dq
                      - pow_2_v_sq * d3h.dq2_dF
                      + dh.dq * v_sq * d2v_sq.dq_dF) / pow_3_v_sq;


      //(vsq^2 Dt[hq, {F, 2}]
      //  - 2 vsq Dt[hq, F] Dt[vsq, F]
      //  + 2 hq Dt[vsq, F]^2
      // - hq vsq Dt[vsq, {F, 2}])/vsq^3

      d2f.q.dF2 = (pow_2_v_sq * d3h.dq_dF2
                   - 2 * v_sq * d2h.dq_dF * dv_sq.dF
                   + 2 * dh.dq * pow<2>(dv_sq.dF)
                   - dh.dq * v_sq * d2v_sq.dF2) / pow_3_v_sq;
      return output;
    }

    OneFormAndFirstDerivatives calculate_beta (double p, const FieldAndDerivatives& fieldAndDerivatives)
    {
      OneFormAndFirstDerivatives output{};

      using boost::math::double_constants::two_pi;

      auto &[beta, dbeta] = output;

      const auto &[f, df, d2f] = fieldAndDerivatives;

//      const auto& f = fieldAndDerivatives.f;
//      const auto& df = fieldAndDerivatives.df;
//      const auto& d2f = fieldAndDerivatives.d2f;

      //                    calculate beta
      //----------------------------------------------------

      //OneForm[2 p \[Pi] Dt[fq, p], 2 \[Pi] (fp + p Dt[fq, q])]
      beta.p = two_pi * p * df.q.dp;
      beta.q = two_pi * (f.p + p * df.q.dq);


      //                    calculate derivatives
      //-----------------------------------------------------

      // 2 \[Pi] (Dt[fq, p] + p Dt[fq, {p, 2}])
      dbeta.p.dp = two_pi * (df.q.dp + p * d2f.q.dp2);

      //2 p \[Pi] Dt[fq, p, q]
      dbeta.p.dq = two_pi * p * d2f.q.dp_dq;

      //2 p \[Pi] Dt[fq, F, p]
      dbeta.p.dF = two_pi * p * d2f.q.dp_dF;

      //2 \[Pi] (Dt[fp, p] + Dt[fq, q] + p Dt[fq, p, q])
      dbeta.q.dp = two_pi * (df.p.dp + df.q.dq + p * d2f.q.dp_dq);

      //2 \[Pi] (Dt[fp, q] + p Dt[fq, {q, 2}])
      dbeta.q.dq = two_pi * (df.p.dq + p * d2f.q.dq2);

      //2 \[Pi] (Dt[fp, F] + p Dt[fq, F, q])
      dbeta.q.dF = two_pi * (df.p.dF + p * d2f.q.dq_dF);

      return output;
    }

    OneFormAndFirstDerivatives
    calculate_gamma (double p,
                     const FieldAndDerivatives& fieldAndDerivatives,
                     const FirstDerivatives& dh,
                     const SecondDerivatives& d2h,
                     const ThirdDerivatives& d3h)
    {
      OneFormAndFirstDerivatives output{};

      using boost::math::double_constants::two_pi;

      auto &[gamma, dgamma] = output;

//      const auto &[f, df, d2f] = fieldAndDerivatives;

      const auto& f = fieldAndDerivatives.f;
      const auto& df = fieldAndDerivatives.df;
      const auto& d2f = fieldAndDerivatives.d2f;

      //                  calculate gamma
      //------------------------------------------------------

      //OneForm[Hf p Dt[fq, p] + fq p Dt[Hf, p], fp Hf + Hf p Dt[fq, q] + fq p Dt[Hf, q]]
      gamma.p = dh.dF * p * df.q.dp + f.q * p * d2h.dp_dF;
      gamma.q = f.p * dh.dF + dh.dF * p * df.q.dq + f.q * p * d2h.dq_dF;


      //                   calculate the derivatives
      //--------------------------------------------------------



      //Hf p Dt[fq, {p, 2}] + Dt[fq, p] (Hf + 2 p Dt[Hf, p])
      // + fq (Dt[Hf, p] + p Dt[Hf, {p, 2}])
      dgamma.p.dp = dh.dF * p * d2f.q.dp2 + df.q.dp * (dh.dF + 2 * p * d2h.dp_dF)
                    + f.q * (d2h.dp_dF + p * d3h.dp2_dF);


      //p (Dt[fq, q] Dt[Hf, p] + Dt[fq, p] Dt[Hf, q] + Hf Dt[fq, p, q]
      //   + fq Dt[Hf, p, q])
      dgamma.p.dq = p * (df.q.dq * d2h.dp_dF + df.q.dp * d2h.dq_dF + dh.dF * d2f.q.dp_dq
                         + f.q * d3h.dp_dq_dF);

      //p (Dt[fq, p] Dt[Hf, F] + Dt[fq, F] Dt[Hf, p] + Hf Dt[fq, F, p]
      //  + fq Dt[Hf, F, p])
      dgamma.p.dF = p * (df.q.dp * d2h.dF2 + df.q.dF * d2h.dp_dF + dh.dF * d2f.q.dp_dF
                         + f.q * d3h.dp_dF2);


      //Hf Dt[fp, p] + fp Dt[Hf, p] + Dt[fq, q] (Hf + p Dt[Hf, p])
      // + fq Dt[Hf, q] + p Dt[fq, p] Dt[Hf, q] + Hf p Dt[fq, p, q]
      // + fq p Dt[Hf, p, q]
      dgamma.q.dp = dh.dF * df.p.dp + f.p * d2h.dp_dF + df.q.dq * (dh.dF + p * d2h.dp_dF)
                    + f.q * d2h.dq_dF + p * df.q.dp * d2h.dq_dF + dh.dF * p * d2f.q.dp_dq
                    + f.q * p * d3h.dp_dq_dF;

      //Hf Dt[fp, q] + Hf p Dt[fq, {q, 2}] + (fp + 2 p Dt[fq, q]) Dt[Hf, q]
      // + fq p Dt[Hf, {q, 2}]
      dgamma.q.dq = dh.dF * df.p.dq + dh.dF * p * d2f.q.dq2 + (f.p + 2 * p * df.q.dq) * d2h.dq_dF
                    + f.q * p * d3h.dq2_dF;

      //Hf Dt[fp, F] + (fp + p Dt[fq, q]) Dt[Hf, F]
      // + p (Dt[fq, F] Dt[Hf, q] + Hf Dt[fq, F, q] + fq Dt[Hf, F, q])
      dgamma.q.dF = dh.dF * df.p.dF + (f.p + p * df.q.dq) * d2h.dF2
                    + p * (df.q.dF * d2h.dq_dF + dh.dF * d2f.q.dq_dF + f.q * d3h.dq_dF2);

      return output;

    }
    OneForm
    calculate_beta1 (const FieldAndDerivatives& f_and_df,
                     const OneFormAndFirstDerivatives& beta_and_dbeta)
    {
      OneForm beta1{};

      const auto& f = f_and_df.f;
      const auto& df = f_and_df.df;

      const auto& beta = beta_and_dbeta.g;
      const auto& dbeta = beta_and_dbeta.dg;



      //\[Beta]p Dt[fp, p] + \[Beta]q Dt[fq, p] + fp Dt[\[Beta]p, p]
      // + fq Dt[\[Beta]p, q]
      beta1.p = beta.p * df.p.dp + beta.q * df.q.dp + f.p * dbeta.p.dp
                + f.q * dbeta.p.dq;


      //\[Beta]p Dt[fp, q] + \[Beta]q Dt[fq, q] + fp Dt[\[Beta]q, p]
      // + fq Dt[\[Beta]q, q]
      beta1.q = beta.p * df.p.dq + beta.q * df.q.dq + f.p * dbeta.q.dp
                + f.q * dbeta.q.dq;

      return beta1;
    }
    OneForm
    calculate_beta2 (const FieldAndDerivatives& f_and_df,
                     const OneFormAndFirstDerivatives& beta_and_dbeta,
                     const FirstDerivatives& dh,
                     const SecondDerivatives& d2h)
    {
      OneForm beta2{};

      const auto& f = f_and_df.f;
      const auto& df = f_and_df.df;

      const auto& beta = beta_and_dbeta.g;
      const auto& dbeta = beta_and_dbeta.dg;


      //-Hf \[Beta]p Dt[fp, p] - Hf \[Beta]q Dt[fq, p]
      // - fp \[Beta]p Dt[Hf, p] + Dt[\[Beta]p, F] - fp Hf Dt[\[Beta]p, p]
      // - fq (\[Beta]q Dt[Hf, p] + Hf Dt[\[Beta]p, q])
      beta2.p = -dh.dF * beta.p * df.p.dp - dh.dF * beta.q * df.q.dp
                - f.p * beta.p * d2h.dp_dF + dbeta.p.dF - f.p * dh.dF * dbeta.p.dp
                - f.q * (beta.q * d2h.dp_dF + dh.dF * dbeta.p.dq);

      //-Hf \[Beta]p Dt[fp, q] - Hf \[Beta]q Dt[fq, q]
      // - fp \[Beta]p Dt[Hf, q] - fp Hf Dt[\[Beta]q, p]
      // - fq (\[Beta]q Dt[Hf, q] + Hf Dt[\[Beta]q, q])
      beta2.q = -dh.dF * beta.p * df.p.dq - dh.dF * beta.q * df.q.dq
                - f.p * beta.p * d2h.dq_dF - f.p * dh.dF * dbeta.q.dp
                - f.q * (beta.q * d2h.dq_dF + dh.dF * dbeta.q.dq);

      return beta2;
    }

    OneForm
    calculate_gamma1 (const FieldAndDerivatives& f_and_df,
                      const OneFormAndFirstDerivatives& gamma_and_dgamma)
    {
      OneForm gamma1{};

      const auto& f = f_and_df.f;
      const auto& df = f_and_df.df;

      const auto& gamma = gamma_and_dgamma.g;
      const auto& dgamma = gamma_and_dgamma.dg;


      //\[Gamma]p Dt[fp, p] + \[Gamma]q Dt[fq, p] + fp Dt[\[Gamma]p, p]
      // + fq Dt[\[Gamma]p, q]
      gamma1.p = gamma.p * df.p.dp + gamma.q * df.q.dp + f.p * dgamma.p.dp
                 + f.q * dgamma.p.dq;

      //\[Gamma]p Dt[fp, q] + \[Gamma]q Dt[fq, q] + fp Dt[\[Gamma]q, p]
      // + fq Dt[\[Gamma]q, q]
      gamma1.q = gamma.p * df.p.dq + gamma.q * df.q.dq + f.p * dgamma.q.dp
                 + f.q * dgamma.q.dq;

      return gamma1;

    }

    OneForm
    calculate_gamma2 (const FieldAndDerivatives& f_and_df,
                      const OneFormAndFirstDerivatives& gamma_and_dgamma,
                      const FirstDerivatives& dh,
                      const SecondDerivatives& d2h)
    {
      OneForm gamma2{};

      const auto& f = f_and_df.f;
      const auto& df = f_and_df.df;

      const auto& gamma = gamma_and_dgamma.g;
      const auto& dgamma = gamma_and_dgamma.dg;

      //-Hf \[Gamma]p Dt[fp, p] - Hf \[Gamma]q Dt[fq, p]
      // - fp \[Gamma]p Dt[Hf, p] + Dt[\[Gamma]p, F] - fp Hf Dt[\[Gamma]p, p]
      // - fq (\[Gamma]q Dt[Hf, p] + Hf Dt[\[Gamma]p, q])
      gamma2.p = -dh.dF * gamma.p * df.p.dp - dh.dF * gamma.q * df.q.dp
                 - f.p * gamma.p * d2h.dp_dF + dgamma.p.dF - f.p * dh.dF * dgamma.p.dp
                 - f.q * (gamma.q * d2h.dp_dF + dh.dF * dgamma.p.dq);

      //-Hf \[Gamma]p Dt[fp, q] - Hf \[Gamma]q Dt[fq, q]
      // - fp \[Gamma]p Dt[Hf, q] - fp Hf Dt[\[Gamma]q, p]
      // - fq (\[Gamma]q Dt[Hf, q] + Hf Dt[\[Gamma]q, q])
      gamma2.q = -dh.dF * gamma.p * df.p.dq - dh.dF * gamma.q * df.q.dq
                 - f.p * gamma.p * d2h.dq_dF - f.p * dh.dF * dgamma.q.dp
                 - f.q * (gamma.q * d2h.dq_dF + dh.dF * dgamma.q.dq);

      return gamma2;

    }

}