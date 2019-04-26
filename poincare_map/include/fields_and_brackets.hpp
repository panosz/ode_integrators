//
// Created by Panagiotis Zestanakis on 25/04/19.
//

#ifndef ODE_INTEGRATORS_FIELDS_AND_BRACKETS_HPP
#define ODE_INTEGRATORS_FIELDS_AND_BRACKETS_HPP
namespace DS
{

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

    struct VelocitySqAndDerivatives {

        double v_Sq = 0;
        FirstDerivatives dv_Sq{};
        SecondDerivatives d2v_Sq{};
    };

    VelocitySqAndDerivatives
    calculate_vSq_and_derivatives (const FirstDerivatives& dh, const SecondDerivatives& d2h, const ThirdDerivatives& d3h);

    struct FieldFirstDerivatives {
        FirstDerivatives p{};
        FirstDerivatives q{};
    };

    struct FieldSecondDerivatives
    {
        SecondDerivatives p{};
        SecondDerivatives q{};
    };

    ///TODO: Calculate Field Second Derivatives

    struct FieldAndDerivatives {
        Field f{};
        FieldFirstDerivatives df{};
        FieldSecondDerivatives d2f{};
        FieldAndDerivatives (const Field& f, const FieldFirstDerivatives& df, const FieldSecondDerivatives& d2f);
        FieldAndDerivatives() = default;
    };

    struct OneFormFirstDerivatives {
        FirstDerivatives p{};
        FirstDerivatives q{};
    };

    struct OneFormAndFirstDerivatives {
        OneForm g{};
        OneFormFirstDerivatives dg{};
        OneFormAndFirstDerivatives (const OneForm& G, const OneFormFirstDerivatives& dG);
    };

    FieldAndDerivatives
    caluclate_translation_field_and_derivatives (const FirstDerivatives& dh,
                                                 const SecondDerivatives& d2h,
                                                 const ThirdDerivatives& d3h);

    OneForm calculate_beta (double p, const FieldAndDerivatives& fieldAndDerivatives);

    OneForm calculate_gamma (double p,
                             const FieldAndDerivatives& fieldAndDerivatives,
                             const FirstDerivatives& dh,
                             const SecondDerivatives& d2h);

}
#endif //ODE_INTEGRATORS_FIELDS_AND_BRACKETS_HPP
