import matplotlib.pyplot as plt
import numpy as np
import action_integration_ext as ai
from fractions import Fraction
from panos_utilities import roots
from pprint import pprint

p_vals = np.linspace(1.98, 2.03, 500)
F = 1.
M = 1.

list_of_points = [np.array([p, 0, F, 0]) for p in p_vals]


@np.vectorize
def get_g_factor(p):
    s = np.array([p, 0, F, 0])
    numer_ho = ai.integrate_E_Pendulum(s, mass=M)
    return numer_ho.g_factor()


#  for s in list_of_points:
    #  numer_ho = ai.integrate_E_H_O(s, mass=M)
    #  print('G = {}'.format(numer_ho.g_factor()))

@np.vectorize
def get_omega_theta(p):
    s = np.array([p, 0, F, 0])
    numer_ho = ai.integrate_E_Pendulum(s, mass=M)
    return numer_ho.omega_theta()


@np.vectorize
def get_omega_phi(p):
    s = np.array([p, 0, F, 0])
    numer_ho = ai.integrate_E_Pendulum(s, mass=M)
    return numer_ho.omega_phi()


ratio_set = {Fraction(x, y) for x in range(-5, 0) if x != 0 for y in range(1, 5) if y >= abs(x)}
ratio_set_float = [float(ratio) for ratio in ratio_set]

crossing_list = roots.find_iso_points(get_g_factor,
                                      ratio_set_float,
                                      window=(0.0001, 2.6), n_samples=25)

for crossing in crossing_list:
    pprint(type(crossing))
pprint(crossing_list)


iso_crosses = dict(zip(ratio_set, crossing_list))

for level, crossings in iso_crosses.items():
    root_string = ','.join([str(sol.root) for sol in crossings.solutions])
    print('crossing {} at {}'.format(level, root_string))
plt.plot(p_vals, get_omega_phi(p_vals))
plt.plot(p_vals, get_omega_theta(p_vals))
plt.show()
