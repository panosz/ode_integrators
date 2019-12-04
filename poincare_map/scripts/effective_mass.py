
import numpy as np
from fractions import Fraction
import collections


ResonancePosition = collections.namedtuple('ResonancePosition', 'F p')

resonance_input = \
    [
        {
            'energy': '-0.5',
            'resonances':
                [Fraction(*fract)
                    for fract in ((-4, 3),
                                  (-5, 4),
                                  (-1, 1),
                                  (-4, 5),
                                  (-3, 4),
                                  (-2, 3),
                                  (-3, 5),
                                  (-1, 2))],
            'resonance_positions':
                [
                    ResonancePosition(F, p)
                    for F, p in
                    ((0.5349761847062395, 0.26448510243958745),
                     (0.5765494982622348, 0.39127866862949423),
                     (0.7518374183554419, 0.7097005260748254),
                     (0.9924282625378186, 0.9923993778089732),
                     (1.0782243052829135, 1.0753830064520393),
                     (1.2585854173979147, 1.2317348881946266),
                     (1.4510101910081645, 1.3791375500711773),
                     (1.8736471872720744, 1.6574964176565055))
                ],
            'hessians':
                [np.array(hes) for hes in
                    ([[-0.126571, 0.683637], [0.683637, -0.0306803]],
                     [[-0.128276, 0.65864],  [0.65864, -0.0580747]],
                     [[-0.134001, 0.577567], [0.577567, -0.114056]],
                     [[-0.139375, 0.503866], [0.503866, -0.12979]],
                     [[-0.140849, 0.483783], [0.483783, -0.129605]],
                     [[-0.143431, 0.448462], [0.448462, -0.125619]],
                     [[-0.145615, 0.418264], [0.418264, -0.119146]],
                     [[-0.149085, 0.369005], [0.369005, -0.104125]])]

        },
        {
            'energy': '0.5',
            'resonances':
                [Fraction(*fract)
                    for fract in ((-1, 4),
                                  (-1, 5),
                                  (0, 1),
                                  (1, 5),
                                  (1, 4),
                                  (1, 3),
                                  (2, 5),
                                  (1, 2),
                                  (3, 5),
                                  (3, 5),
                                  (1, 2),
                                  (2, 5),
                                  (1, 3),
                                  (1, 4),
                                  (1, 5))],
            'resonance_positions':
                [ResonancePosition(F, p) for F, p in
                    ((1.406426524950831,   1.952652823699508),
                     (1.0927634731311728,  1.784804456029384),
                     (0.7666012890999908,  1.5916037755044379),
                     (0.6568018494201529,  1.52105348322809),
                     (0.6392548378070139,  1.5094733106663487),
                     (0.6150744426203811,  1.4933683019405368),
                     (0.5992585180259321,  1.4827397061021412),
                     (0.5799496838439424,  1.4696596094633223),
                     (0.5646776183925167,  1.4592310429760715),
                     (0.47702818405859365, 1.3978756626099431),
                     (0.46243209552304043, 1.3873947495381698),
                     (0.43753019394353193, 1.3693284441240035),
                     (0.41136257000852217, 1.3500833826164385),
                     (0.3610667306194962,  1.3123008272644623),
                     (0.31725878277166963, 1.2784825245357636))],
            'hessians':
                [np.array(hes) for hes in
                    ([[-0.197221, 0.450264],   [0.450264, -0.28819]],
                     [[-0.21261, 0.522504],    [0.522504, -0.413183]],
                     [[-0.259943, 0.672601],   [0.672601, -0.740351]],
                     [[-0.310688, 0.789832],   [0.789832, -1.02051]],
                     [[-0.324847, 0.819294],   [0.819294, -1.09016]],
                     [[-0.350156, 0.869997],   [0.869997, -1.20809]],
                     [[-0.372159, 0.91262],    [0.91262, -1.30513]],
                     [[-0.408567, 0.981213],   [0.981213, -1.45726]],
                     [[-0.449738, 1.05695],    [1.05695, -1.62004]],
                     [[0.13407, 0.00771191],   [0.00771191, 0.0168597]],
                     [[0.0631353, 0.154334],   [0.154334, -0.270189]],
                     [[0.0155928, 0.26686],    [0.26686, -0.523033]],
                     [[-0.00648129, 0.332726], [0.332726, -0.709041]],
                     [[-0.0256894, 0.418346],  [0.418346, -1.04633]],
                     [[-0.0334166, 0.481992],  [0.481992, -1.40626]])]

        }
    ]

ResonanceDetails = collections.namedtuple('ResonanceDetails',
                                          'ratio position hessian')

resonance_dict = {}

_fields = ('resonances', 'resonance_positions', 'hessians')

for resonance_entry in resonance_input:
    _energy_data = (resonance_entry[f] for f in _fields)
    _energy = resonance_entry['energy']
    resonance_dict[_energy] = [ResonanceDetails(*data)
                               for data in zip(*_energy_data)]


def resonance_rotation_Jacobian(resonance: Fraction) -> np.array:
    """
    :param resonance: Fraction
        the resonance ratio expressed as a fraction of co-prime integers
    :return: np.array:
        The Jacobian of the transform to rotated coordinates
        that eliminate the resonance
    """
    r = resonance.numerator
    s = resonance.denominator

    return np.array([[r, 0],
                     [-s, 1]])


def to_resonance_hessian(resonance: Fraction,
                         original_hessian) -> np.array:
    """
    Transforms the Hessian of an Action function in the original coordinates
    to the Hessian in the rotated coordinates
    :param resonance: Fraction,
    the resonance ratio expressed as a fraction of co-prime integers
    :param original_hessian:  array-like
    The Hessian of an Action Function,
    usually the Hamiltonian, in the original coordinates
    :return: The hessian in the rotated coordinates
    """
    Jac = resonance_rotation_Jacobian(resonance)

    return Jac.transpose() @ original_hessian @ Jac


def effective_mass(resonance, hessian):
    """
    returns the effective mass of the pendulum-like motion
    around a fixed point on a resonance
    :param resonance:
        a fraction of co-prime integers
    :param hessian: array-like
        2x2 symmetric array shaped as [[hjj hjf],[hjf,hff]]
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


#  def get_harmonics_on_resonant_orbit(F, p, nperturb, mperturb, n_samples=32):
    #  def standard_perturb(p, q, chi):
        #  return np.sin(nperturb * q + mperturb * chi)

    #  standard_perturb.m = mperturb
    #  orbit = ds.make_closed_orbit(ds.HarmonicOsillator(F),
    #                               p0=p,
    #                               t_integration_max=1000)
    #  return harmonics.orbital_harmonics_on_const_chi_bar(standard_perturb,
    #                                                      orbit,
    #                                                      n_samples)


NominalWidth = collections.namedtuple('NominalWidth',
                                      'delta_J delta_F')
ResonanceOutput = collections.namedtuple('ResonanceOutput',
                                         'position nominal_widths')


# def get_nominal_resonance_widths(nperturb, mperturb, energy):
#     if energy not in resonance_dict.keys():
#         raise RuntimeError(f"No entries for Energy={energy}")
#
#     resonance_widths = []
#
#     for res_details in resonance_dict[energy]:
#         ratio_res = res_details.ratio
#         F_res = res_details.position.F
#         p_res = res_details.position.p
#         hessian_res = res_details.hessian
#
#         har = get_harmonics_on_resonant_orbit(F=F_res,
#                                               p=p_res,
#                                               nperturb=nperturb,
#                                               mperturb=mperturb)
#         k_res = k_on_resonance(ratio_res, mperturb, har)
#
#         if k_res:
#             eff_mass = effective_mass(ratio_res, hessian_res)
#             delta_J_hat_nominal = 2 * np.sqrt(k_res / eff_mass)
#
#             r = ratio_res.numerator
#             s = ratio_res.denominator
#
#             delta_J = abs(r) * delta_J_hat_nominal
#             delta_F = abs(s) * delta_J_hat_nominal
#             position = res_details.position
#             nominal_width = NominalWidth(delta_J=delta_J, delta_F=delta_F)
#             res_width = ResonanceOutput(position=position,
#                                         nominal_widths=nominal_width)
#             resonance_widths.append(res_width)
#
#     return resonance_widths
#

# if __name__ == "__main__":
#
#     amplitude = 0.04
#     nperturb = -5
#     mperturb = 2
#
#     res_n_width = get_nominal_resonance_widths(nperturb=nperturb,
#                                                mperturb=mperturb)
#
#     for (res, nominal_widths) in res_n_width.items():
#         delta_J = np.sqrt(amplitude) * nominal_widths.delta_J
#         delta_F = np.sqrt(amplitude) * nominal_widths.delta_F
#         print("widths at {}: delta_J = {}, delta_F = {}".format(res,
#                                                                 delta_J,
#                                                                 delta_F))
#
#     F_min = resonance_dict[res].position.F - delta_F
#
#     F_max = resonance_dict[res].position.F + delta_F
