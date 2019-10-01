
import numpy as np
from fractions import Fraction
import collections

import dynamic_system as ds
import harmonics

ResonancePosition = collections.namedtuple('ResonancePosition', 'F p')

resonance_input = \
    [
        {
            'energy': '-0.5',
            'resonances':
                [
                    Fraction(-4, 3), Fraction(-5, 4),
                    Fraction(-1, 1), Fraction(-4, 5),
                    Fraction(-3, 4), Fraction(-2, 3),
                    Fraction(-3, 5), Fraction(-1, 2)
                ],
            'resonance_positions':
                [
                    ResonancePosition(F=0.5349761847062395, p=0.26448510243958745),
                    ResonancePosition(F=0.5765494982622348, p=0.39127866862949423),
                    ResonancePosition(F=0.7518374183554419, p=0.7097005260748254),
                    ResonancePosition(F=0.9924282625378186, p=0.9923993778089732),
                    ResonancePosition(F=1.0782243052829135, p=1.0753830064520393),
                    ResonancePosition(F=1.2585854173979147, p=1.2317348881946266),
                    ResonancePosition(F=1.4510101910081645, p=1.3791375500711773),
                    ResonancePosition(F=1.8736471872720744, p=1.6574964176565055)
                ],
            'hessians':
                [
                    np.array([[-0.126571, 0.683637], [0.683637, -0.0306803]]),
                    np.array([[-0.128276, 0.65864], [0.65864, -0.0580747]]),
                    np.array([[-0.134001, 0.577567], [0.577567, -0.114056]]),
                    np.array([[-0.139375, 0.503866], [0.503866, -0.12979]]),
                    np.array([[-0.140849, 0.483783], [0.483783, -0.129605]]),
                    np.array([[-0.143431, 0.448462], [0.448462, -0.125619]]),
                    np.array([[-0.145615, 0.418264], [0.418264, -0.119146]]),
                    np.array([[-0.149085, 0.369005], [0.369005, -0.104125]])
                ]

        },
        {
            'energy': '0.5',
            'resonances':
                [
                    Fraction(-1, 4),
                    Fraction(-1, 5),
                    Fraction(0, 1),
                    Fraction(1, 5),
                    Fraction(1, 4),
                    Fraction(1, 3),
                    Fraction(2, 5),
                    Fraction(1, 2),
                    Fraction(3, 5),
                    Fraction(3, 5),
                    Fraction(1, 2),
                    Fraction(2, 5),
                    Fraction(1, 3),
                    Fraction(1, 4),
                    Fraction(1, 5)
                ],
            'resonance_positions':
                [
                    ResonancePosition(F=1.406426524950831, p=1.952652823699508),
                    ResonancePosition(F=1.0927634731311728, p=1.784804456029384),
                    ResonancePosition(F=0.7666012890999908, p=1.5916037755044379),
                    ResonancePosition(F=0.6568018494201529, p=1.52105348322809),
                    ResonancePosition(F=0.6392548378070139, p=1.5094733106663487),
                    ResonancePosition(F=0.6150744426203811, p=1.4933683019405368),
                    ResonancePosition(F=0.5992585180259321, p=1.4827397061021412),
                    ResonancePosition(F=0.5799496838439424, p=1.4696596094633223),
                    ResonancePosition(F=0.5646776183925167, p=1.4592310429760715),
                    ResonancePosition(F=0.47702818405859365, p=1.3978756626099431),
                    ResonancePosition(F=0.46243209552304043, p=1.3873947495381698),
                    ResonancePosition(F=0.43753019394353193, p=1.3693284441240035),
                    ResonancePosition(F=0.41136257000852217, p=1.3500833826164385),
                    ResonancePosition(F=0.3610667306194962, p=1.3123008272644623),
                    ResonancePosition(F=0.31725878277166963, p=1.2784825245357636)
                ],
            'hessians':
                [
                    np.array([[-0.197221, 0.450264], [0.450264, -0.28819]]),
                    np.array([[-0.21261, 0.522504], [0.522504, -0.413183]]),
                    np.array([[-0.259943, 0.672601], [0.672601, -0.740351]]),
                    np.array([[-0.310688, 0.789832], [0.789832, -1.02051]]),
                    np.array([[-0.324847, 0.819294], [0.819294, -1.09016]]),
                    np.array([[-0.350156, 0.869997], [0.869997, -1.20809]]),
                    np.array([[-0.372159, 0.91262], [0.91262, -1.30513]]),
                    np.array([[-0.408567, 0.981213], [0.981213, -1.45726]]),
                    np.array([[-0.449738, 1.05695], [1.05695, -1.62004]]),
                    np.array([[0.13407, 0.00771191], [0.00771191, 0.0168597]]),
                    np.array([[0.0631353, 0.154334], [0.154334, -0.270189]]),
                    np.array([[0.0155928, 0.26686], [0.26686, -0.523033]]),
                    np.array([[-0.00648129, 0.332726], [0.332726, -0.709041]]),
                    np.array([[-0.0256894, 0.418346], [0.418346, -1.04633]]),
                    np.array([[-0.0334166, 0.481992], [0.481992, -1.40626]])
                ]

        }
    ]

ResonanceDetails = collections.namedtuple('ResonanceDetails', 'ratio position hessian')

resonance_dict = {}

for resonance_entry in resonance_input:
    _resonance_positions = resonance_entry['resonance_positions']
    _hessians = resonance_entry['hessians']
    _resonances = resonance_entry['resonances']
    _res_details = [ResonanceDetails(ratio=r_ratio, position=r_pos, hessian=hes) for r_ratio, r_pos, hes in
                    zip(_resonances, _resonance_positions, _hessians)]

    resonance_dict[resonance_entry['energy']] = _res_details


def resonance_rotation_Jacobian(resonance: Fraction) -> np.array:
    """
    :param resonance: Fraction, the resonance ratio expressed as a fration of coprime integers
    :return: np.array: The Jacobian of the transform to rotated coordinates that eliminate the resonance
    """
    r = resonance.numerator
    s = resonance.denominator

    return np.array([[r, 0],
                     [-s, 1]])


def to_resonance_hessian(resonance: Fraction, original_hessian: np.array) -> np.array:
    """
    Transforms the Hessian of an Action function in the original coordinates to the Hessian in the rotated coodrdinates
    :param resonance: Fraction, the resonance ratio expressed as a fration of coprime integers
    :param original_hessian: The Hessian of an Action Function, usually the Hamiltonian, in the original coordinates
    :return: The hessian in the rotated coordinates
    """
    Jac = resonance_rotation_Jacobian(resonance)

    return Jac.transpose() @ original_hessian @ Jac


def effective_mass(resonance, hessian):
    """
    returns the effective mass of the pendulum-like motion around a fixed point on a resonance
    :param resonance: a fraction of co-prime integers
    :param hessian: 2x2 symmetric np.array shaped as [[hjj hjf],[hjf,hff]]
    """

    return np.abs(to_resonance_hessian(resonance, hessian)[0, 0])


def k_on_resonance(res, m_harmonic, An_s):
    r = res.numerator
    s = res.denominator

    if m_harmonic % s:
        return 0
    else:
        p = m_harmonic // s
        return 2 * abs(An_s[-p * r])


def get_harmonics_on_resonant_orbit(F, p, nperturb, mperturb, n_samples=32):
    def standard_perturb(p, q, chi):
        return np.sin(nperturb * q + mperturb * chi)

    standard_perturb.m = mperturb
    orbit = ds.make_closed_orbit(ds.HarmonicOsillator(F), p0=p, t_integration_max=1000)
    return harmonics.orbital_harmonics_on_const_chi_bar(standard_perturb, orbit, n_samples)


NominalWidth = collections.namedtuple('NominalWidth', 'delta_J delta_F')
ResonanceOutput = collections.namedtuple('ResonanceOutput','position nominal_widths')


def get_nominal_resonance_widths(nperturb, mperturb, energy):
    if not energy in resonance_dict.keys():
        raise RuntimeError(f"No entries for Energy={energy}")

    resonance_widths = []

    for res_details in resonance_dict[energy]:
        ratio_res = res_details.ratio
        F_res = res_details.position.F
        p_res = res_details.position.p
        hessian_res = res_details.hessian

        har = get_harmonics_on_resonant_orbit(F=F_res, p=p_res, nperturb=nperturb, mperturb=mperturb)
        k_res = k_on_resonance(ratio_res, mperturb, har)

        if k_res:
            delta_J_hat_nominal = 2 * np.sqrt(k_res / effective_mass(ratio_res, hessian_res))

            r = ratio_res.numerator
            s = ratio_res.denominator

            delta_J = abs(r) * delta_J_hat_nominal
            delta_F = abs(s) * delta_J_hat_nominal

            resonance_widths.append(ResonanceOutput(position=res_details.position,
                                                    nominal_widths=NominalWidth(delta_J=delta_J, delta_F=delta_F)))

    return resonance_widths


# if __name__ == "__main__":
#
    # amplitude = 0.04
    # nperturb = -5
    # mperturb = 2
    #
    # for (res, nominal_widths) in get_nominal_resonance_widths(nperturb=nperturb, mperturb=mperturb).items():
    #     delta_J = np.sqrt(amplitude) * nominal_widths.delta_J
    #     delta_F = np.sqrt(amplitude) * nominal_widths.delta_F
    #     print("widths at {}: delta_J = {}, delta_F = {}".format(res,
    #                                                             delta_J,
    #                                                             delta_F))
    #
    #     F_min, F_max = resonance_dict[res].position.F - delta_F, resonance_dict[res].position.F + delta_F
