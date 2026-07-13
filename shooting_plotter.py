import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from collections import defaultdict
from itertools import cycle

shooting_data = "output_files/output_shooting.txt"
eigenperturbation_data = "output_files/output_eigenperturbation.txt"

# ---------------------------------------------------------------------------
# Global aesthetic configuration: clean white background, a soft creamy
# grid, and a richer/darker color palette. Curves are drawn with a
# semi-transparent fill so that overlapping regions naturally read as
# darker -- this block only affects *appearance*; none of the data
# processing logic below is touched.
# ---------------------------------------------------------------------------
mpl.rcParams.update({
    "figure.dpi": 100,
    "savefig.dpi": 100,
    "savefig.bbox": "tight",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    #"font.family": "serif",
    "font.size": 12,
    "axes.titlesize": 16,
    "axes.titleweight": "bold",
    "axes.titlecolor": "#2B2B2B",
    "axes.labelsize": 14,
    "axes.labelcolor": "#333333",
    "axes.edgecolor": "#999999",
    "axes.linewidth": 0.9,
    "axes.grid": True,
    "grid.color": "#EDE0C8",
    "grid.linestyle": "-",
    "grid.linewidth": 0.9,
    "grid.alpha": 0.9,
    "legend.frameon": True,
    "legend.framealpha": 0.92,
    "legend.fancybox": True,
    "legend.fontsize": 10,
    "legend.edgecolor": "#CCCCCC",
    "legend.facecolor": "white",
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "lines.linewidth": 2.2,
    "lines.solid_capstyle": "round",
    "lines.dash_capstyle": "round",
    "axes.prop_cycle": mpl.cycler(color=[
        "#1B4F72",  # deep blue
        "#922B21",  # deep red
        "#1D6E50",  # deep green
        "#6C3483",  # deep purple
        "#9A6D00",  # deep mustard
        "#0E6E6E",  # deep teal
        "#7D3C98",  # plum
        "#4A4A4A",  # charcoal
    ]),
})

# Toggle this to turn the peak-pointing arrows on/off everywhere.
SHOW_PEAK_ARROWS = True

# Fill opacity for the area under each curve. Kept moderate so that
# overlapping curves visibly darken where they cross.
FILL_ALPHA = 0.0

# Accent colors for the spike/peak markers and the intercept line.
SPIKE_COLOR = "#C0392B"
SPIKE_EDGE = "white"
INTERCEPT_COLOR = "#555555"



def add_spikes(x, y, ax, show_arrows=SHOW_PEAK_ARROWS):
    return 0

    """
    Helper method to add spikes to a plot. The spikes are identified namely by
    finding the gradient dy/dx and marking all the gradients whose numerical 
    value is above the n-th percentile of the prescribed dataset. 
    @param x            the x-values of the dataset in question. 
    @param y            the y-values of the dataset in question. 
    @param ax           the plot to which the spikes will be added.
    @param show_arrows  if True, draws an arrow pointing at each marked peak.
    @return      a plot with the spikes added to it.
    """

    # Prefixed parameters for identifying the spikes. 
    gradient_threshold = 20
    value_percentile = 68
    lower_threshold = 0.02
    upper_threshold = 2.95
    
    # Sorts the dataset
    index = np.argsort(x)
    x = x[index]
    y = y[index]

    # Calculates gradient
    gradient_y = np.gradient(y, x)

    # Finds the threshold based on percentiles.
    value_threshold = np.percentile(gradient_y, value_percentile)
    
    # Sets the mask up.
    spike_mask = (np.abs(gradient_y) > gradient_threshold) #& (y > value_threshold)

    # Placeholders for the data to be scattered.
    scatter_x, scatter_y = [], []

    # Loops over the mask. When mask applies, adds to the points to scatter.
    in_spike = False
    for i, value in enumerate(spike_mask):
        # Also adds an arbitrary threshold beyond which spikes are to be 
        # regarded as spurious.
        if (value and not in_spike and 
            abs (x[i]) > lower_threshold and 
            abs (x[i]) < upper_threshold):
            if (x[i] < 0):
                # If the gradient was in [i], mark the [i+1] point.
                scatter_x.append (x [i+1])
                scatter_y.append (y [i+1])
            else:
                # The other way around when signs invert. 
                scatter_x.append (x [i])
                scatter_y.append (y [i])
            in_spike = True
        elif not value:
            in_spike = False

    ax.plot(
        scatter_x, scatter_y,
        marker="o", linestyle="None",
        markersize=10, markerfacecolor=SPIKE_COLOR,
        markeredgecolor=SPIKE_EDGE, markeredgewidth=1.6,
        zorder=6,
    )

    # Optionally draws an arrow pointing straight at each marked peak,
    # offset slightly above it so the marker itself stays clearly visible.
    if show_arrows and len(scatter_x) > 0:
        y_range = (max(y) - min(y)) or 1.0
        for px, py in zip(scatter_x, scatter_y):
            ax.annotate(
                "",
                xy=(px, py),
                xytext=(px, py + 0.18 * y_range),
                arrowprops=dict(
                    arrowstyle="->",
                    color=SPIKE_COLOR,
                    lw=1.6,
                    shrinkA=0,
                    shrinkB=4,
                ),
                zorder=6,
            )

    return ax

def generate_output (file_name, show_arrows=SHOW_PEAK_ARROWS):

    """
    Given a filename (either the output from the shooting of from the 
    eigenperturbation method), this method will read and process the data on
    the file in question, and plot the data. 
    @param file_name    the file name. 
    @param show_arrows  if True, draws arrows pointing at the detected peaks
                         (only relevant for the shooting plots).
    @return             the plot with the data and labels attached to it. 
    """

    # Read data from file, remove empty lines if they emerge, and assign 
    # each row to their respective groups, using the first four columns
    # as label. 
    data = np.genfromtxt(file_name, delimiter=',')
    data = data[~np.isnan(data).any(axis=1)]
    groups = defaultdict(list)
    for row in data:
        key = tuple (row [:4])  
        groups [key].append (row)

    # Creates instance of the plot.
    fig, ax = plt.subplots(figsize=(8.5, 6.2))
    
    # Some markers to distinguish the lines ... 
    line_cycler = cycle([':', '--',  (0, (3, 1, 1, 1, 1, 1)), '-'])

    # Loops over the different groups and plots (x, y) for each.
    for key, rows in groups.items():

        rows = np.array(rows)
        x = rows[:, -2] 
        y = rows[:, -1] / max (rows [:,-1])
        #y = 1/ rows[:, -1] * (2 - rows [:, -3]) / rows [:,0] * 4
        # Distinguish the label depending on what we are plotting.
        label = ""

        if "shooting" in file_name:
            label = f"d={key[0]}, $\\eta$={key[1]}, s={key[2]}"
            
        elif "eigenperturbation" in file_name:
            label = f"d={key[0]}, $\\eta$={key[1]}, s={key[2]}, $\\sigma$={key[3]}"
    
        # Actually plots: draws the line on top of a semi-transparent fill,
        # so overlapping curves visibly darken where they cross.
        line, = ax.plot(
            x, y,
            # ".",
            linestyle=next(line_cycler),
            label=label,
            linewidth=2.2,
            solid_capstyle="round",
            zorder=3,
        )
        ax.fill_between(x, 0, y, color=line.get_color(), alpha=FILL_ALPHA,
                          zorder=2, linewidth=0)

        # For the shooting graph, also want to add the spikes. 
        if "shooting" in file_name:
            add_spikes (x, y, ax, show_arrows=show_arrows)

    #Also distinguishes labels and title
    xlabel, ylabel, title = "", "", ""

    if "shooting" in file_name:
        title = "Potential fixed point solutions"
        xlabel = r"$\sigma$"
        ylabel = r"$\rho_\infty ^{(\sigma)}$"

    else:
        title = "RG Eigendirections"
        xlabel = r"Eigenvalue $y$"
        ylabel = r"Asymptotic eigenvector $\nu_\infty$"

        # Also includes the intercept in this case.
        ax.plot([-3, 3], [0, 0], color=INTERCEPT_COLOR, linewidth=2.2,
                 linestyle="-", zorder=1, alpha=0.7)
        ax.plot([0,0], [-1, 1], color=INTERCEPT_COLOR, linewidth=2.2,
                 linestyle="-", zorder=1, alpha=0.7)

    # Adds labels, titles and grid to plot, and return.
    ax.set_ylabel(ylabel, fontsize=14, labelpad=10)
    ax.set_xlabel(xlabel, fontsize=14, labelpad=10)
    ax.set_title(title, fontsize=16, pad=14)

    # Creamy grid, kept behind the data and fills.
    ax.grid(True, which="major", zorder=0)
    ax.set_axisbelow(True)

    # Clean white-background frame.
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=5, width=0.9, colors="#444444")

    # Legend kept inside the plot area.
    ax.legend(loc="best")

    fig.tight_layout()

    return fig, ax


from pathlib import Path
import matplotlib.pyplot as plt

# output_dir = Path("output_files/shooting-eff-dim-changed-formula")
output_dir = Path("output_files/eigenperturbations-long-range/")

for file_path in output_dir.glob("*.txt"):
    print(f"Processing {file_path}")

    # if not "shooting" in str(file_path):
    #     continue

    fig, ax = generate_output(str(file_path), show_arrows=SHOW_PEAK_ARROWS)

    # Optional: save the figure
    fig.savefig(file_path.with_suffix(".png"), dpi=300)

    plt.show()      # or remove this if you only want saved figures
    plt.close(fig)  # prevents memory buildup