import numpy as np
import h5py
import matplotlib.pyplot as plt

from panos_poincare_plot import poincare as poinc


def wrap2pi(x):
    return x - 2 * np.pi * np.floor_divide(x, 2 * np.pi)


hf = h5py.File('cross.hdf5')

fig, ax = plt.subplots(1, 1)


list_array = [ np.array(dataset['cross_points']) for (_,dataset) in hf.items()]

a=np.concatenate(list_array,axis=0)

poinc.plot_poincare_data(ax, wrap2pi(a[:,3]),  a[:,2], x_label='$\chi$', y_label='$p$')

ax.set_xlim(-.01, 2 * np.pi + 0.1)
plt.show()

