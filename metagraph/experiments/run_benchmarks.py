#!/usr/bin/env python3

"""
Run this script from build.

Example:
../experiments/run_benchmarks.py flat norepl 20
../experiments/run_benchmarks.py brwt_custom norepl 20 --arity 5

"""


import numpy as np
import sys
from glob import glob
import os
import glob
import subprocess
from multiprocessing import Pool


import sys
import numpy as np
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
from scipy.special import comb

mpl.rcParams['figure.dpi'] = 300
rcParams.update({'errorbar.capsize' : 2})
rcParams['text.antialiased'] = True
rcParams['font.family'] = "Helvetica"
plot_ratio = (1 + np.sqrt(5)) / 2
markers = ['^', 'v', '<', '>', 'o', 's', 'P', '*']

DIR = "./simulate"

command = (os.path.dirname(os.path.realpath(__file__))) \
        + "/../build/run_experiments matrices simulate {simulation} --anno-type {annotator}" \
        + " -d {density} -m {num_columns} -n {num_rows}"

query_command = (os.path.dirname(os.path.realpath(__file__))) \
              + "/../build/run_experiments query {type} --anno-type {annotator}" \
              + " -n {num_samples}"

num_rows = 10**6
num_unique_count = 5
num_init_fraction = 10

num_unique_rows = int(num_rows / num_unique_count)
num_init_rows = int(num_rows / num_init_fraction)
decay = 0.2

max_density = 0.01
d_array = np.linspace(0, max_density, num=13)
methods = [
    "column",
    "flat",
    "rbfish",
    "brwt",
    "bin_rel_wt_sdsl",
    "bin_rel_wt"
]

methods_to_plot = [
    "column",
    "flat",
    "rbfish",
    "bin_rel_wt_sdsl",
    "bin_rel_wt",
    "brwt_arity_2",
    "brwt_greedy_relax_7",
]

methods_name_map = {
    "column": "Column",
    "flat": "Flat",
    "rbfish": "Rainbowfish",
    "bin_rel_wt_sdsl": "BinRel-WT (SDSL)",
    "bin_rel_wt": "BinRel-WT",
    "brwt_arity_2": "BRWT2",
    "brwt_arity_7": "BRWT7",
    "brwt_greedy": "BRWT2 (greedy)",
    "brwt_greedy_relax_7": "BRWT7 (greedy)",
}

# simulations = [ "norepl", "uniform", "weighted" ]
# simulation_names = {
#     "norepl": "Random columns",
#     "uniform": "Uniformly distributed",
#     "weighted": "Weighted"
#}

simulations = [ "norepl", "uniform_rows", "uniform_columns",
#"weighted_rows"
]
simulation_names = {
    "norepl": "Random columns",
    "uniform_rows": "Duplicated rows",
    "uniform_columns": "Duplicated columns",
    #"weighted_rows": "Weighted duplicate rows",
    #"weighted_columns": "Weighted duplicate columns",
}

queries = [ "point", "row", "column" ]
num_samples = {
    "point": 100000,
    "row": 1000,
    "column": 10
}


def run_experiment(command, cwd = DIR):
    print(command)
    completed_process = subprocess.run(command.split(), cwd=cwd)
    if completed_process.returncode:
        print("ERROR {} occured".format(completed_process.returncode))


def get_size(signature):
    f = glob.glob(signature)
    if len(f) != 1:
        return 0

    return os.path.getsize(f[0])

def plot_query(query_files):
    lines = np.array([line.rstrip().split(" ")
        for query_file in query_files for line in open(query_file,"r")])
    data = lines[:,[2,3,4,6,7]].astype(np.double)
    sorter = np.argsort(data[:,0])
    data = data[sorter,:]
    lines = lines[sorter,:]

    exp_types = np.unique(lines[:,1])
    for exp_type in exp_types:
        exp_selector = lines[:,1] == exp_type
        query_types = np.unique(lines[:,5])
        # TODO: space between subplots
        f, axes = plt.subplots(ncols=len(query_types), sharex=False, figsize=(16,6))
        if len(query_types) == 1:
            axes = [axes]

        for i, (ax, query_type) in enumerate(zip(axes, queries)):
            if query_type not in query_types:
                continue

            query_selector = (lines[:,5] == query_type) & exp_selector
            num_query = 0
            compressor_types = np.unique(lines[query_selector,0])
            for j, ctype in enumerate(methods):
                if ctype not in compressor_types:
                    continue

                #print("Plotting",exp_type, query_type, ctype)
                type_selector = (lines[:,0] == ctype) & query_selector
                type_selector = type_selector & (data[:,0] <= 0.02)
                num_query = data[query_selector,3].astype(np.int)[0]
                ax.plot(data[type_selector,0] * 100, data[type_selector,4],
                        label=ctype,
                        marker=markers[(i * len(queries) + j) % len(markers)],
                        ms=3)
            ax.set_title("{}, {} query, n = {}".format(
                exp_type, query_type, num_query))
            ax.set_ylabel("Query time (ms)")
            ax.set_xlabel("Density (%)")
            ax.grid(True)
            ax.legend()

        plt.tight_layout()
        plt.savefig(
            'query_time_{}.pdf'.format(exp_type),
            fmt='pdf')
        plt.show()


def create_plots(methods, num_columns, d_array):
    simulation_arrays = []
    for simulation in simulations:
        size_array = []
        for method in methods_to_plot:
            if method == "brwt":

                params = ["greedy"]
                for arity in range(2, 12):
                    params.append("arity_{}".format(arity))

                partitionings = params[:]
                for relax in range(2, 12):
                    for part in partitionings:
                        params.append(part + "_relax_{}".format(relax))

                for param in params:
                    sizes = [get_size(DIR + "/simulate.{simulation}_{d}*_{m}_brwt_{param}.annomat".format(
                                simulation=simulation,
                                d="{:.6f}".format(d)[:7],
                                n=num_rows,
                                m=num_columns,
                                param=param
                            )) for d in d_array]
                    if np.unique(sizes).size > 1 and (method + "_" + param) in methods:
                        size_array.append((method + "_" + param, np.array(sizes)))
            else:
                sizes = [get_size(DIR + "/simulate.{simulation}_{d}*_{m}_{method}.annomat".format(
                            simulation=simulation,
                            d="{:.6f}".format(d)[:7],
                            n=num_rows,
                            m=num_columns,
                            method=method
                        )) for d in d_array]
                if method in methods:
                    size_array.append((method, np.array(sizes)))

        simulation_arrays.append((simulation, size_array))
    f, axes = plt.subplots(ncols=len(simulation_arrays),
                           sharex=False,
                           sharey=False,
                           figsize=(len(simulation_arrays) * plot_ratio * 3, 4))
    for i, (ax, (simulation, size_array)) in enumerate(zip(axes, simulation_arrays)):
        for j, (method, method_array) in enumerate(size_array):
            ax.plot(d_array * 100,
                    method_array * 8 / num_rows / num_columns,
                    label=method,
                    marker=markers[(i * len(simulation_arrays) + j) % len(markers)],
                    ms=5)
        ax.set_xlabel('Density (%)')
        ax.set_ylabel('Bits per matrix element')
        ax.set_title('{}: n={:.1e}, m={:.1e}'.format(
            simulation_names[simulation],
            num_rows,
            num_columns))
        ax.grid(True)
        ax.legend(loc='best', fontsize=8)

    plt.tight_layout()
    plt.savefig(
        'rrr_mat_shape_{}_{}_drange_{}_{}.pdf'.format(
            min(d_array), max(d_array),
            num_rows, num_columns
        ),
        fmt='pdf')
    plt.show()



def main():
    if len(sys.argv) > 1 and sys.argv[1][:4] == "plot":
        if len(sys.argv) < 3:
            print("Usage:\n{} plot <num_columns> [methods_to_plot ...]".format(sys.argv[0]))
            exit(1)

        print("Start plotting...")

        print("Compression")
        num_columns = int(sys.argv[2])
        methods_to_plot = sys.argv[3:]
        if len(methods_to_plot) == 0:
            methods_to_plot = methods[:]
        create_plots(methods_to_plot, num_columns, d_array)

        print("Query")
        plot_query(glob.glob(DIR + "/simulate.*.query_stats"))

        print("Done")

        return

    if len(sys.argv) < 5:
        print("Usage:\n{} <compression_method> <simulation_method> <num_columns> <num_threads> [method_params ...]".format(sys.argv[0]))
        print("{} plot <num_columns> [methods_to_plot ...]".format(sys.argv[0]))
        exit(1)

    input_methods, input_simulations, num_columns, num_threads = sys.argv[1:5]
    num_columns = int(num_columns)
    num_unique_columns = int(num_columns / num_unique_count)
    method_params = ' '.join(sys.argv[5:])

    if not os.path.exists(DIR):
        os.mkdir(DIR)

    if input_methods == ".":
        input_methods = methods
    else:
        input_methods = [input_methods]

    if input_simulations == ".":
        input_simulations = simulations
    else:
        input_simulations = [input_simulations]

    commands_to_run = []
    for method in input_methods:
        if method in methods:
            for simulation in input_simulations:
                if simulation in simulations:
                    if simulation == "norepl":
                        commands_to_run += [ command.format(
                            simulation=simulation,
                            annotator=method,
                            density=d,
                            num_rows=num_rows,
                            num_columns=num_columns) +
                                " " +
                                method_params for d in d_array
                        ]
                    elif simulation == "uniform_rows":
                        commands_to_run += [ command.format(
                            simulation=simulation,
                            annotator=method,
                            density=d,
                            num_rows=num_rows,
                            num_columns=num_columns) +
                                " -u {}".format(num_unique_rows) +
                                " " + method_params for d in d_array
                        ]
                    elif simulation == "uniform_columns":
                        commands_to_run += [ command.format(
                            simulation=simulation,
                            annotator=method,
                            density=d,
                            num_rows=num_rows,
                            num_columns=num_columns) +
                                " -u {}".format(num_unique_columns) +
                                " " + method_params for d in d_array
                        ]
                    elif simulation == "weighted_rows":
                        commands_to_run += [ command.format(
                            simulation=simulation,
                            annotator=method,
                            density=d,
                            num_rows=num_rows,
                            num_columns=num_columns) +
                                " -u {} -c {} -D {}".format(num_unique_rows,
                                                            num_init_rows,
                                                            decay) +
                                " " + method_params for d in d_array
                        ]
                elif simulation in queries:
                    if method != "brwt":
                        files = [file[(len(DIR) + 1):]
                            for sim in simulations
                            for file in glob.glob(
                                DIR + "/simulate.{}_*_{}.annomat".format(sim, method)
                            )
                        ]
                    else:
                        files = [file[(len(DIR) + 1):]
                            for sim in simulations
                            for file in glob.glob(
                                DIR + "/simulate.{}_*_{}_*.annomat".format(sim, method)
                            )
                        ]
                    commands_to_run += [ query_command.format(
                        type=simulation,
                        annotator=method,
                        num_samples=num_samples[simulation]) +
                            " " + mat for mat in files
                    ]
                else:
                    print("Simulation type {} unknown, skipping experiments.".format(simulation))
        else:
            print("Method {} unknown, skipping experiments.".format(method))

    pool = Pool(processes=int(num_threads))
    pool.imap_unordered(run_experiment, commands_to_run)
    pool.close()
    pool.join()

    print("Done!")


if __name__ == '__main__':
    main()
