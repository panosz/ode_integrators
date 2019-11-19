import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import host_subplot
from panos_utilities import coprimes
from panos_utilities import roots as panos_roots
from fractions import Fraction
from pprint import pprint
import pandas as pd
import action_integration_ext as ai

_colors = ['#1b9e77', '#d95f02', '#7570b3']

M = 1.
dynamic_system = ai.PendulumDynamicSystem(M)
option_dict = {'abs_err': 1e-12, 'rel_err': 1e-12, 'dt': 1e-2}
options = ai.IntegrationOptions(**option_dict)


def get_p(f, energy):
    p = np.sqrt(2 * (energy + f))
    return p


def points_on_energy_surface(*, f, energy):
    p = get_p(f, energy)
    if np.isscalar(f):
        return np.array([p, 0, f, 0])
    else:
        points = np.zeros([f.size, 4])
        points[:, 0] = get_p(f, energy)
        points[:, 2] = f
        return points


@np.vectorize
def g_of_f(f, energy):
    s = points_on_energy_surface(f=f, energy=energy)
    action_integrals = dynamic_system.action_integrals(s=s,
                                                       time=1000,
                                                       options=options)
    return action_integrals.g_factor()


def initial_points_and_g_on_energy_surface(energy, n_points, f_window):
    fmin, fmax = f_window
    f = np.linspace(fmin, fmax, n_points)
    g = g_of_f(f, energy)

    initial_points = points_on_energy_surface(f=f, energy=energy)
    return initial_points, g


def make_ratio_set(max_int):
    """
    returns a list of positive and negative ratios constructed from integers in
    the interval [0 max_int]
    """
    positive_ratio_set = list(coprimes.coprime_fractions(max_int))
    negative_ratio_set = [- x for x in positive_ratio_set]
    ratio_set = sorted(negative_ratio_set +
                       [Fraction(0, 1)] +
                       positive_ratio_set)
    return ratio_set


def resonances_on_energy_surface(energy, ratio_set, f_window, n_samples=15):
    crossings = panos_roots.find_roots(lambda x: g_of_f(x, energy),
                                       ratio_set,
                                       window=f_window, n_samples=n_samples)

    valid_crossings = (crossing for crossing in crossings
                       if crossing.solutions)

    output_list = []

    for vc in valid_crossings:
        level = vc.level
        f_values = vc.solutions
        init_points = [points_on_energy_surface(f=f, energy=energy)
                       for f in f_values]
        action_integrals = [dynamic_system.action_integrals(s=s,
                                                            time=1000,
                                                            options=options)
                            for s in init_points]

        output_list.append((level, f_values, init_points, action_integrals))

    return pd.DataFrame(output_list,
                        columns=['ratio',
                                 'F_list',
                                 'init_points',
                                 'action_integrals'])


def relative_x(ax, x):
    '''
        calculates the relative value of x with respect to the ax axes
    '''
    x0, x1 = ax.get_xlim()
    return (x - x0) / (x1 - x0)


def relative_y(ax, y):
    '''
        calculates the relative value of y with respect to the ax axes
    '''
    y0, y1 = ax.get_ylim()
    return (y - y0) / (y1 - y0)


def plot_single_resonance_level(ax,
                                resonance_ratio,
                                sorted_resonance_positions,
                                cut_line_style,
                                resonance_position_style):
    '''
        plots the positions of a resonance ratio
    '''
    resonance_position_max = sorted_resonance_positions[-1]
    relative_y_max = relative_y(ax, resonance_position_max)
    relative_x_min = relative_x(ax, resonance_ratio)

    ax.axvline(resonance_ratio,
               ymax=relative_y_max,
               **cut_line_style)

    for position in sorted_resonance_positions:
        ax.axhline(position,
                   xmin=relative_x_min,
                   **cut_line_style)

    ax.plot([resonance_ratio] * len(sorted_resonance_positions),
            sorted_resonance_positions, **resonance_position_style)


def plot_resonances(ax, resonance_data):
    resonance_line_style = {'color': 'black',
                            'alpha': 0.5,
                            'linestyle': '--'}

    resonance_position_style = {'linestyle': '',
                                'marker': 'x',
                                'color': _colors[1]}
    x_ticks = []
    y_ticks = []
    for _, resonance_row in resonance_data.iterrows():
        resonance_ratio = resonance_row['ratio']
        resonance_positions = sorted(resonance_row['F_list'])

        x_ticks.append(resonance_ratio)
        y_ticks += resonance_positions

        plot_single_resonance_level(ax,
                                    resonance_ratio,
                                    resonance_positions,
                                    resonance_line_style,
                                    resonance_position_style)

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


def run_example(ax,
                energy,
                n_points,
                f_window,
                y_labels_no_print_interval,
                filenames, ratio_max_int):
    initial_points, g = initial_points_and_g_on_energy_surface(energy,
                                                               n_points,
                                                               f_window)

    np.savetxt(filenames['initial_points_filename'],
               initial_points,
               delimiter=" ")

    ratio_list = make_ratio_set(ratio_max_int)

    resonance_data = resonances_on_energy_surface(energy,
                                                  ratio_list,
                                                  f_window,
                                                  n_samples=15)

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

    # if there is a separatrix,
    # avoid overlapping resonance y labels near the separatrix.
    if float(energy) > 0 and y_labels_no_print_interval:
        label_list = ax2.get_yticklabels()
        y_no_print_min, y_no_print_max = y_labels_no_print_interval
        for label in label_list:
            y = float(label.get_text())
            if y_no_print_min < y < y_no_print_max:
                label.set_visible(False)


if __name__ == "__main__":
    Energy = 'P'

    params = {
                'energy': 0.5,
                'n_points': 40,
                'f_window': (0.3001, 1.5),
                'y_labels_no_print_interval': (0.5, 0.67),
                'ratio_max_int': 4,
                'n': -5,
                'm': 4,
                'plot_initial_points': False
              }

    init_points, g = initial_points_and_g_on_energy_surface(params['energy'],
                                                            params['n_points'],
                                                            params['f_window'])

    sns.set_style("ticks")

    fig = plt.gcf()

    ax = host_subplot(121)
    ax_right = host_subplot(122)

    #  run_and_plot.run_and_plot_poincare(ax_right, **poincare_params)

    #  ax.set_ylim(ax_right.get_ylim())
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30)

    ax_right.get_yaxis().set_visible(False)

    fig.set_size_inches(12.5, 6.5, forward=True)
    ax.plot(g, init_points[:, 2])

    ratio_list = make_ratio_set(params['ratio_max_int'])

    resonance_data = resonances_on_energy_surface(params['energy'],
                                                  ratio_list,
                                                  params['f_window'],
                                                  n_samples=15)
    pprint([[pos_int.hessian() for pos_int in action_integrals]
           for action_integrals in resonance_data['action_integrals'].tolist()])

    ax2 = plot_resonances(ax, resonance_data)

    plt.show()
