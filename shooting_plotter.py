import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from itertools import cycle

shooting_data = "output-files/output_shooting.txt"
eigenperturbation_data = "output-files/output_eigenperturbation.txt"



def add_spikes(x, y, ax):

    """
    Helper method to add spikes to a plot. The spikes are identified namely by
    finding the gradient dy/dx and marking all the gradients whose numerical 
    value is above the n-th percentile of the prescribed dataset. 
    @param x     the x-values of the dataset in question. 
    @param y     the y-values of the dataset in question. 
    @param ax    the plot to which the spikes will be added.
    @return      a plot with the spikes added to it.
    """

    # Prefixed parameters for identifying the spikes. 
    gradient_threshold = 1000
    value_percentile = 68
    lower_threshold = 0.025
    upper_threshold = 0.68
    
    # Sorts the dataset
    index = np.argsort(x)
    x = x[index]
    y = y[index]

    # Calculates gradient
    gradient_y = np.gradient(y, x)

    # Finds the threshold based on percentiles.
    value_threshold = np.percentile(gradient_y, value_percentile)
    
    # Sets the mask up.
    spike_mask = (np.abs(gradient_y) > gradient_threshold) & (y > value_threshold)

    # Placeholders for the data to be scattered.
    scatter_x, scatter_y = [], []

    # Loops over the mask. When mask applies, adds to the points to scatter.
    in_spike = False
    for i, value in enumerate(spike_mask):
        # Also adds an arbitrary threshold beyond which spikes are to be 
        # regarded as spurious.
        if (
            value and not in_spike and 
            abs (x[i]) > lower_threshold and 
            abs (x[i]) < upper_threshold):
            # If the gradient was in [i], mark the [i+1] point.
            scatter_x.append (x[i+1])
            scatter_y.append (y[i+1])
            in_spike = True
        elif not value:
            in_spike = False

    ax.plot (scatter_x, scatter_y, "k.", markersize=10)

    return ax

def generate_output (file_name):

    """
    Given a filename (either the output from the shooting of from the 
    eigenperturbation method), this method will read and process the data on
    the file in question, and plot the data. 
    @param file_name    the file name. 
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
    fig, ax = plt.subplots (figsize = (8,6))
    
    # Some markers to distinguish the lines ... 
    line_cycler = cycle([':', '--',  (0, (3, 1, 1, 1, 1, 1)), '-'])

    # Loops over the different groups and plots (x, y) for each.
    for key, rows in groups.items():

        rows = np.array(rows)
        x = rows[:, -2] 
        y = rows[:, -1] 

        # Distinguish the label depending on what we are plotting.
        label = ""

        if "shooting" in file_name:
            label = f"d={key[0]}, $\\eta$={key[1]}, s={key[2]}, N={key[3]}"
            
        elif "eigenperturbation" in file_name:
            label = f"d={key[0]}, $\\eta$={key[1]}, s={key[2]}, $\\sigma$={key[3]}"
    
        # Actually plots.
        ax.plot(x, y, color="k",  linestyle=next (line_cycler), label=label, alpha=0.7)

        # For the shooting graph, also want to add the spikes. 
        if "shooting" in file_name:
            add_spikes (x, y, ax)

    #Also distinguishes labels and title
    xlabel, ylabel, title = "", "", ""

    if "shooting" in file_name:
        title = "Potential fixed point solutions"
        xlabel = r"$\sigma$"
        ylabel = r"$\rho_\infty ^{(\sigma)}$"

    elif "eigenperturbation" in file_name:
        title = "RG Eigendirections"
        xlabel = r"Eigenvalue $y$"
        ylabel = r"Asymptotic eigenvector $\nu_\infty$"

        # Also includes the intercept in this case.
        ax.plot ([-5,0], [0,0], "k-", linewidth=2.5)

    # Adds labels, titles and grid to plot, and return.
    ax.set_ylabel (ylabel, fontsize=14)
    ax.set_xlabel (xlabel, fontsize=14)
    ax.set_title (title,fontsize=14)
    ax.grid ()
    ax.legend ()
    return fig, ax

# Generates the plot(s).
ax, fig = generate_output (shooting_data)
plt.show ()
ax, fig = generate_output (eigenperturbation_data)
plt.show ()