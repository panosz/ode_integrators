
from pathlib import Path
import subprocess
from collections import ChainMap
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from plot_hdf5_output import plot_hdf5_output
import effective_mass_on_resonance as effective_mass
import seaborn as sns


def stringify_list(input_list):
    return [str(x) for x in input_list]


def make_arg_list(command, input_file, integration_time, perturbation_amplitude, q_harmonic, chi_harmonic, output_file):
    return stringify_list(
        [command, input_file, integration_time, perturbation_amplitude, q_harmonic, chi_harmonic, output_file])


def add_dictionary_info_to_hdf5(file, info_set_name, **kwargs):
    file[info_set_name] = info_set_name
    for key, value in kwargs.items():
        file[info_set_name].attrs[key] = value


command_path = str(Path("bin/generalized_pendulum").resolve())


def run_my_command(**kwargs):
    default = {'command': command_path,
               'input_file': "my_foo.txt",
               'output_file': 'cross'}

    arguments = ChainMap(kwargs, default)

    arg_list = make_arg_list(**arguments)

    proc = subprocess.run(args=arg_list)

    init_points = np.loadtxt(arguments['input_file'])

    output_file_name = arguments['output_file'] + '.hdf5'

    with h5py.File(output_file_name, 'r+') as f:
        add_dictionary_info_to_hdf5(f, 'command_arguments', **arguments)
        f['init_points'] = init_points

    return proc


def run_and_plot_poincare(ax, time, amplitude, n, m, energy, plot_initial_points):
    output_file = f'data/cross={amplitude}_n={n}_m={m}'
    run_my_command(input_file=f'initial_points@energy={energy}.txt', integration_time=time,
                   perturbation_amplitude=amplitude,
                   q_harmonic=n, chi_harmonic=m,
                   output_file=output_file)

    sns.set_style("ticks")

    nominal_widths_output = effective_mass.get_nominal_resonance_widths(nperturb=n, mperturb=m, energy=energy)

    colors = sns.color_palette("Paired", len(nominal_widths_output))

    for (res, color) in zip(nominal_widths_output, colors):
        nominal_widths = res.nominal_widths
        delta_J = np.sqrt(amplitude) * nominal_widths.delta_J
        delta_F = np.sqrt(amplitude) * nominal_widths.delta_F
        print("widths at {}: delta_J = {}, delta_F = {}".format(res,
                                                                delta_J,
                                                                delta_F))
        F = res.position.F

        F_min, F_max = F - delta_F, F + delta_F

        plt.axhspan(F_min, F_max, facecolor=color, alpha=0.8, edgecolor=color)

    if float(energy) > 0:
        ax.axhline(float(energy), color='black', linestyle='-', linewidth=3)

    plot_hdf5_output(ax, output_file, plot_initial_points)
    ax.set_xlabel(r'$\bar{\phi}$')
    ax.set_ylabel(r'$F$')


if __name__ == "__main__":
    time = 200000
    amplitude = 3e-3
    n = -5
    m = 4
    energy = "0.5"
    fig, ax = plt.subplots(1, 1)

    run_and_plot_poincare(ax, time, amplitude, n, m, energy, plot_initial_points=True)

    figure_file = f'figures/poincare_E={energy}_A={amplitude}_n={n}_m={m}.jpeg'

    # plt.savefig(figure_file)

    plt.show()
