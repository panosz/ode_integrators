
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

from panos_poincare_plot import poincare as poinc

def wrap2pi( x):
    return x - 2*np.pi*np.floor_divide(x,2*np.pi)


a = sp.loadtxt('../cmake-build-release/cross_exact.txt')


fig, ax = plt.subplots(1, 1)

poinc.plot_poincare_data(ax, wrap2pi(a[:,2]), a[:,0], x_label='$\chi$', y_label='$p$')

ax.set_xlim(-.01,2*np.pi+0.1)

plt.show()
