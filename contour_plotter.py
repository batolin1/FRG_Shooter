from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
from scipy.stats import rankdata


def process_file(filepath):
    print(f"Processing {filepath}")

    df = pd.read_csv(
        filepath,
        header=None,
        names=[
            "dimension",
            "s_factor",
            "anomalous_dimension",
            "sigma",
            "asymptotic_field",
        ],
    )

    # # -----------------------
    # # Plot 1
    # # -----------------------
    # plt.figure(figsize=(8, 6))

    # sc = plt.scatter(
    #     df["sigma"],
    #     df["asymptotic_field"],
    #     c=df["anomalous_dimension"],
    #     cmap="YlGn",
    #     s=10,
    # )

    # plt.xlabel("Sigma")
    # plt.ylabel("Asymptotic Field")
    # plt.colorbar(sc, label="Anomalous Dimension")
    # plt.tight_layout()
    # plt.show()

    z = df["asymptotic_field"]
    z_norm = rankdata(z) / len(z)
    # z_norm = z

    # -----------------------
    # Plot 2
    # -----------------------
    plt.figure(figsize=(8, 6))

    sc = plt.scatter(
        df["sigma"],
        df["anomalous_dimension"],
        c=z_norm,
        cmap="YlGn",
        s=10,
    )

    plt.xlabel("Sigma")
    plt.ylabel("Anomalous Dimension")
    plt.colorbar(sc, label="Asymptotic Field")
    plt.tight_layout()
    plt.show()

    # -----------------------
    # 3D scatter
    # -----------------------
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")

    sc = ax.scatter(
        df["sigma"],
        df["anomalous_dimension"],
        df["asymptotic_field"],
        c=df["asymptotic_field"],
        cmap="viridis",
        s=10,
    )

    ax.set_xlabel("Sigma")
    ax.set_ylabel("Anomalous Dimension")
    ax.set_zlabel("Asymptotic Field")

    fig.colorbar(sc, label="Asymptotic Field")
    plt.show()

    # # -----------------------
    # # Interpolated surface
    # # -----------------------
    # sigma_grid = np.linspace(df["sigma"].min(), df["sigma"].max(), 500)
    # eta_grid = np.linspace(
    #     df["anomalous_dimension"].min(),
    #     df["anomalous_dimension"].max(),
    #     500,
    # )

    # X, Y = np.meshgrid(sigma_grid, eta_grid)

    # Z = griddata(
    #     (df["sigma"], df["anomalous_dimension"]),
    #     df["asymptotic_field"],
    #     (X, Y),
    #     method="linear",
    # )

    # fig = plt.figure(figsize=(10, 8))
    # ax = fig.add_subplot(111, projection="3d")

    # ax.scatter(
    #     X.ravel(),
    #     Y.ravel(),
    #     Z.ravel(),
    #     c=Z.ravel(),
    #     cmap="viridis",
    #     s=10,
    # )

    # sc = ax.scatter(
    #     df["sigma"],
    #     df["anomalous_dimension"],
    #     df["asymptotic_field"],
    #     c=df["asymptotic_field"],
    #     cmap="viridis",
    #     s=10,
    # )

    # ax.set_xlabel("Sigma")
    # ax.set_ylabel("Anomalous Dimension")

    # fig.colorbar(sc, ax=ax, label="Asymptotic Field")
    # plt.show()


# Process every .txt file in the folder
folder = Path("output_files/testers/annealing_fractional_dimensions")

for filepath in folder.glob("*.txt"):
    process_file(filepath)