import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import host_subplot, make_axes_locatable
from panos_utilities import coprimes
from panos_utilities import roots as panos_roots
from fractions import Fraction
from pprint import pprint
import pandas as pd
import action_integration_ext as ai

_colors = ['#1b9e77', '#d95f02', '#7570b3']

M = 1.


def get_p(f, energy):
    p = np.sqrt(2 * (energy + f))
    return p


def g_of_f(f, energy):
    p = get_p(f, energy)
    s = np.array([p, 0, f, 0])
    numer_ho = ai.integrate_E_Pendulum(s, mass=M)
    return numer_ho.g_factor()


def g_of_f_vectorized(f, energy):
    if np.isscalar(f):
        return g_of_f(f, energy)
    else:
        return np.array([g_of_f(myf, energy) for myf in f])


def initial_points_on_energy_surface(energy, n_points, f_window):
    fmin, fmax = f_window
    f = np.linspace(fmin, fmax, n_points)
    g = g_of_f_vectorized(f, energy)

    initial_points = np.zeros([n_points, 4])
    initial_points[:, 0] = get_p(f, energy)
    initial_points[:, 2] = f
    return initial_points, g


def make_ratio_set(max_int):
    """
    returns a list of positive and negative ratios constructed from integers in the interval [0 max_int]
    """
    positive_ratio_set = list(coprimes.coprime_fractions(max_int))
    negative_ratio_set = [- x for x in positive_ratio_set]
    ratio_set = sorted(negative_ratio_set + [Fraction(0, 1)] + positive_ratio_set)
    return ratio_set


def resonances_on_energy_surface(energy, ratio_set, f_window, n_samples=15):
    crossings = panos_roots.find_roots(lambda x: g_of_f_vectorized(x, energy),
                                       ratio_set,
                                       window=f_window, n_samples=n_samples)

    return pd.DataFrame((crossing for crossing in crossings if crossing.solutions)
                        , columns=['ratio', 'F_list'])


def plot_resonances(ax, resonance_data):
    def relative_x(x):
        x0, x1 = ax.get_xlim()
        return (x - x0) / (x1 - x0)

    def relative_y(y):
        y0, y1 = ax.get_ylim()
        return (y - y0) / (y1 - y0)

    x_ticks = []
    y_ticks = []
    for _, resonance_row in resonance_data.iterrows():
        resonance_ratio = resonance_row['ratio']
        resonance_positions = sorted(resonance_row['F_list'])

        x_ticks.append(resonance_ratio)
        y_ticks += resonance_positions

        resonance_position_max = resonance_positions[-1]
        ax.axvline(resonance_ratio, color='black', alpha=0.5, ymax=relative_y(resonance_position_max), linestyle='--')

        for position in resonance_positions:
            ax.axhline(position, color='black', alpha=0.5, xmin=relative_x(resonance_ratio), linestyle='--')

        ax.plot([resonance_ratio] * len(resonance_positions), resonance_positions, 'x', color=_colors[1])

    y_ticks = sorted(y_ticks)
    y_tick_labels = ["{:.2f}".format(x) for x in y_ticks]
    x_tick_labels = [str(x_tick) for x_tick in x_ticks]

    ax2 = ax.twin()
    ax2.set_yticks(y_ticks)
    ax2.set_yticklabels(y_tick_labels)
    ax2.axis["top"].major_ticklabels.set_visible(False)
    ax2.axis["top"].major_ticks.set_visible(False)
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_tick_labels)
    return ax2


def resonance_data_from_json(filename):
    resonance_data = pd.read_json(filename)
    resonance_data['ratio'] = \
        resonance_data['ratio'].transform(lambda x: Fraction(x['numerator'], x['denominator']))
    return resonance_data


def cache_resonance_data(filename, energy, ratio_list, f_window, n_samples=15):
    try:
        resonance_data = resonance_data_from_json(filename)
        print(f'read resonance data from {filename}')
    except ValueError:
        resonance_data = resonances_on_energy_surface(energy, ratio_list, f_window, n_samples=n_samples)
        resonance_data.to_json(filename)
        print("had to calculate a")

    return resonance_data


def run_example(ax, energy, n_points, f_window, y_labels_no_print_interval, filenames, ratio_max_int):
    initial_points, g = initial_points_on_energy_surface(energy, n_points, f_window)

    np.savetxt(filenames['initial_points_filename'], initial_points, delimiter=" ")

    ratio_list = make_ratio_set(ratio_max_int)

    resonance_data = cache_resonance_data(filenames['resonance_filename'], energy, ratio_list, f_window, n_samples=15)

    pprint(resonance_data)

    ax.plot(g, initial_points[:, 2], color=_colors[0])

    if float(energy) > 0:  # if there is a separatrix
        # limit x axis range when there is a separatrix
        left, right = plt.xlim()
        right = min(right, 0.8)
        ax.set_xlim(left, right)

        # plot the separatrix
        ax.axhline(float(energy), color='black', linestyle='-', linewidth=3)

    ax.set_xlabel(r"""$G$""", fontsize=12)
    ax.set_ylabel(r"""$F$""", fontsize=12)

    ax2 = plot_resonances(ax, resonance_data)

    # avoid overlapping resonance y labels near the separatrix, if there is a separatrix
    if float(energy) > 0 and y_labels_no_print_interval:
        label_list = ax2.get_yticklabels()
        for label in label_list:
            if y_labels_no_print_interval[0] < float(label.get_text()) < y_labels_no_print_interval[1]:
                label.set_visible(False)


positive_energy_params = \
    {
        'energy': 0.5,
        'n_points': 40,
        'f_window': (0.3001, 1.5),
        'y_labels_no_print_interval': (0.3, 0.67),
        'ratio_max_int': 5
    }

negative_energy_params = \
    {
        'energy': -0.5,
        'n_points': 40,
        'f_window': (0.5001, 2),
        'y_labels_no_print_interval': None,
        'ratio_max_int': 5
    }

if __name__ == "__main__":
    Energy = 'P'

    params = positive_energy_params if Energy is 'P' else negative_energy_params

    amplitudes_negative = [0.004, 0.008, 0.01, 0.022]
    amplitudes_positive = [0.0008, 0.003, 0.0045]

    for ampl in amplitudes_positive:
        poincare_params = \
            {
                'time': 40000,
                'amplitude': ampl,
                'n': -5,
                'm': 4,
                'energy': str(params['energy']),
                'plot_initial_points': False
            }

        file_names = \
            {
                'initial_points_filename': f'initial_points@energy={params["energy"]}.txt',
                'figure_filename': 'figures/resonances_E_{}A_{}n_{}m_{}.pdf'.format(
                    params["energy"],
                    poincare_params["amplitude"],
                    poincare_params['n'],
                    poincare_params['m']),
                'resonance_filename': f'resonances@energy={params["energy"]},nmax={params["ratio_max_int"]}.json'
            }

        sns.set_style("ticks")

        fig = plt.gcf()

        ax = host_subplot(121)
        ax_right = host_subplot(122)

        #  run_and_plot.run_and_plot_poincare(ax_right, **poincare_params)

        ax.set_ylim(ax_right.get_ylim())
        run_example(ax, filenames=file_names, **params)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=30)

        ax_right.get_yaxis().set_visible(False)

        fig.set_size_inches(12.5, 6.5, forward=True)
        plt.savefig(file_names['figure_filename'], quality=95, transparent=False, dpi=200)
        plt.close(fig)
    # plt.show()
