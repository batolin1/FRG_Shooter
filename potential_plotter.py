import os
import numpy as np
import matplotlib.pyplot as plt


def parse_label_from_filename(file_name):
    """
    Extract parameters from filename and build a LaTeX-style label.
    """

    name = file_name.replace(".txt", "")
    parts = name.split("_")

    params = {}
    for part in parts:
        key, value = part.split("=")
        params[key] = value

    # Match your convention (eigenperturbation-like case)
    label = (
        f"d={params['dimension']}, "
        f"$\\eta$={params['anomalous-dimension']}, "
        f"s={params['s-factor']}, "
        f"$\\sigma$={params['sigma']}"
    )

    return label

def plot_trajectory(file_path):
    """
    Given a trajectory file (produced by the C++ code), this function reads
    the data and plots potential_0prime as a function of the wavefunction.

    @param file_path    path to the trajectory file
    @return             (fig, ax) matplotlib objects
    """

    # Read file (each line is a dataset)
    with open(file_path, "r") as f:
        lines = [line.strip() for line in f.readlines() if line.strip()]

    # Convert lines into numpy arrays (ignore trailing commas)
    data = [np.array([float(x) for x in line.split(",") if x]) for line in lines]

    # Assign variables according to C++ structure
    wavefunction      = data[0]
    potential_0prime  = data[1]
    potential_1prime  = data[2]
    potential_2prime  = data[3]

    # Extract parameters from filename for labeling
    filename = os.path.basename(file_path)
    label = parse_label_from_filename(filename)

    # Create plot
    fig, ax = plt.subplots(figsize=(8, 6))

    ax.plot(
        wavefunction,
        potential_1prime,
        linestyle='-',
        color='k',
        alpha=0.8,
        label=label
    )

    # Labels and styling
    ax.set_xlabel(r"$\tilde{\rho}$", fontsize=14)
    ax.set_ylabel(r"$\tilde{U}^{(1)}$", fontsize=14)
    ax.set_title("Potential Derivative as a function of Field", fontsize=14)

    ax.grid()
    ax.legend()

    return fig, ax


def plot_all_trajectories(folder="output-files/trajectories"):
    """
    Reads all trajectory files in a folder and plots them one by one.

    @param folder    directory containing trajectory files
    """

    for filename in os.listdir(folder):
        if filename.endswith(".txt"):
            file_path = os.path.join(folder, filename)

            fig, ax = plot_trajectory(file_path)

            # Show each plot separately
            plt.show()

plot_all_trajectories (folder="output-files/trajectories")