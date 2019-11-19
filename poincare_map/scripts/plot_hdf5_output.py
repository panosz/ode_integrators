import h5py
from pprint import pprint
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
from panos_poincare_plot import poincare as poinc
from panos_utilities import wrap


def is_orbit_key(k):
    return k.startswith('orbit')


def parse_hdf5_output(filename):
    d = defaultdict(list)

    with h5py.File(filename + '.hdf5', 'r') as f:
        for k in filter(is_orbit_key, f.keys()):
            d['initial_point'].append(f[k]['initial_point'][...])
            d['cross_points'].append((f[k]['cross_points'][...]))

    init = np.concatenate(d['initial_point'])
    cross = np.concatenate(d['cross_points'])
    pprint(init)

    return init, cross


def plot_hdf5_output(axes, filename, plot_initial_points):

    init, cross = parse_hdf5_output(filename)

    poinc.plot_poincare_data(axes, wrap.wrap_2pi(cross[:, 3]), cross[:, 2], x_label='$\chi$', y_label='$p$',
                             axes_to_data=False,alpha=0.8)

    if plot_initial_points:
        poinc.plot_poincare_data(axes, wrap.wrap_2pi(init[:, 3]), init[:, 2], marker='x', markersize=2,
                                 x_label='$\chi$', y_label='$p$', color='r', axes_to_data=False)
