import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict


def parse_label(dim, eta, s, sigma):
    return (
        f"d={dim}, "
        f"$\\eta$={eta}, "
        f"s={s}, "
        f"$\\sigma$={sigma}"
    )


def load_trajectories(file_path):
    """
    Reads ONE file containing MULTIPLE trajectories.
    Groups rows by (dim, eta, s, sigma).
    """

    data = np.loadtxt(file_path, delimiter=",")

    if data.ndim == 1:
        data = data.reshape(1, -1)

    trajectories = defaultdict(lambda: {
        "field": [], 
        "U0": [],
        "U1": [],
        "U2": [],
        "denominator": [],
        "real_denominator": []})

    for row in data:
        dim, eta, s, sigma = row[:4]

        key = (dim, eta, s, sigma)

        trajectories[key]["field"].append(row[-6])
        trajectories[key]["U0"].append(row[-5])
        trajectories[key]["U1"].append(row[-4])
        trajectories[key]["U2"].append(row[-3])
        trajectories[key]["denominator"].append(row[-2])
        trajectories[key]["real_denominator"].append(row[-1])

    return trajectories


def plot_file(file_path):
    trajectories = load_trajectories(file_path)

 

    for (dim, eta, s, sigma), t in trajectories.items():
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.plot(
            t["field"],
            t["U2"],
            ".",
            label=parse_label(dim, eta, s, sigma),
            alpha=0.8
        )

        ax.set_xlabel(r"$\tilde{\rho}$")
        ax.set_ylabel(r"$\tilde{U}^{(0)}$")
        ax.set_title("Potential Derivative vs Field")
        ax.grid()
        ax.legend()
        plt.show ()

plot_file("output_files/output_eigenperturbation_trajectory.txt")
