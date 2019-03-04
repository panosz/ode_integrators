import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

from panos_poincare_plot import poincare as poinc


def wrap2pi(x):
    return x - 2 * np.pi * np.floor_divide(x, 2 * np.pi)


a = sp.loadtxt('cross_exact.txt')

fig, ax = plt.subplots(1, 1)

def energy(orb):
    p=a[:,0]
    q=a[:,1]
    F=a[:,2]
    chi=a[:,3]
    return 0.5*p**2 - F* np.cos(q) + np.sin(q+2*chi)

poinc.plot_poincare_data(ax, wrap2pi(a[:, 3]),  a[:,2], x_label='$\chi$', y_label='$p$')

ax.set_xlim(-.01, 2 * np.pi + 0.1)

plt.show()
