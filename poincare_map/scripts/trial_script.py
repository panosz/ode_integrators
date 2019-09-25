import matplotlib.pyplot as plt
import numpy as np
import action_integration_ext as ai

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


plt.plot(p_vals, get_omega_phi(p_vals))
plt.plot(p_vals, get_omega_theta(p_vals))
plt.show()
