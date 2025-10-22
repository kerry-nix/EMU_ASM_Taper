from matplotlib import colors
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.collections import LineCollection
import pandas as pd


def plot_iEEG_data(data, t, fig=None, ax=None, linecolor="k"):
    """ "
    2021.06.23. Python 3.8
    Akash Pattnaik
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Purpose:
    To plot iEEG data
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Input
        data: iEEG data in pandas.DataFrame or numpy.array
        time: time array
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Output:
        Returns figure handle
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    """

    if fig is None and ax is None:
        fig, ax = plt.subplots()

    # Show only bottom and left axis for visibility, format tick labels
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.ticklabel_format(useOffset=False)

    if data.ndim == 2:
        n_rows = data.shape[1]
    else:
        n_rows = 1

    data = pd.DataFrame(data)
    
    ticklocs = []
    ax.set_xlim(t[0], t[-1])
    dmin = data.min().min()
    dmax = data.max().min()

    dr = (dmax - dmin) * 0.8  # Crowd them a bit.

    y0 = dmin
    y1 = (n_rows - 1) * dr + dmax
    ax.set_ylim(y0, y1)

    segs = []
    for i in range(n_rows):
        if isinstance(data, pd.DataFrame):
            segs.append(np.column_stack((t, data.iloc[:, i])))
        elif isinstance(data, np.ndarray):
            segs.append(np.column_stack((t, data[:, i])))
        else:
            print("Data is not in valid format")

    for i in reversed(range(n_rows)):
        ticklocs.append(i * dr)

    offsets = np.zeros((n_rows, 2), dtype=float)
    offsets[:, 1] = ticklocs

    lines = LineCollection(
        segs, offsets=offsets, transOffset=None, colors=linecolor, linewidth=0.2
    )
    ax.add_collection(lines)

    # # Set the yticks to use axes coordinates on the y axis
    ax.set_yticks(ticklocs)
    if isinstance(data, pd.DataFrame):
        ax.set_yticklabels(data.columns)

    ax.set_xlabel("Time (s)")
    return fig, ax
