"""
Functions for calculating 3D migration of RFs in cartesian coordinates.

Note: Based on codes originally written by Matteo Scarponi.

Location: Chavannes-pres-renens, CH
Date: Mar 2022
Author: Konstantinos Michailos
"""

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from collections import OrderedDict
import matplotlib.patches as patches


def add_legend(ax, fontsize=12):
    fontsize = 10
    h, l = ax.get_legend_handles_labels()
    by_l = OrderedDict(zip(l, h))
    legend = ax.legend(by_l.values(), by_l.keys(), loc="best", fontsize=fontsize)
    frame = legend.get_frame()
    frame.set_facecolor("white")
    frame.set_edgecolor("black")
    frame.set_linewidth(0.5)
    frame.set_alpha(1)

    return legend


def add_colorbar(ax, m, title=False, ticks=False, ticks_Flag=False, visible=True):
    # colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.05)
    if ticks_Flag:
        cbar = plt.colorbar(m, cax=cax, ticks=ticks)
    else:
        cbar = plt.colorbar(m, cax=cax)
    cbar.ax.set_visible(visible)
    if title:
        cbar.ax.set_title(title, fontsize=fontsize)
    cbar.ax.tick_params(axis="both", which="major", labelsize=fontsize)
    return cbar


fontsize = 12
markersize = 100


def plot_migration_profile(Gp, migration_param_dict, sta, work_directory, filename=False):


    # Read migration parameters
    minx = migration_param_dict['minx']
    maxx = migration_param_dict['maxx']
    pasx = migration_param_dict['pasx']
    miny = migration_param_dict['miny']
    maxy = migration_param_dict['maxy']
    pasy = migration_param_dict['pasy']
    minz = migration_param_dict['minz']
    maxz = migration_param_dict['maxz']
    pasz = migration_param_dict['pasz']
    inc = migration_param_dict['inc']
    zmax = migration_param_dict['zmax']

    zz = np.arange(minz, maxz + pasz, pasz)
    xx = np.arange(minx, maxx + pasx, pasx)

    XX, ZZ = np.meshgrid(xx, zz)

    # COLOR PALETTE AND COLORMAP
    # cork, bam, broc, vik
    pal_col = work_directory + "/data/colormaps/vik.txt"
    pal_col = pd.read_csv(pal_col, header=None, index_col=False, sep="\s+", names=["R", "G", "B"])
    cm = LinearSegmentedColormap.from_list("blue2red", pal_col.values, len(pal_col))
    c = np.min([np.max(Gp), 0.1])
    c = 0.06
    CL = 2

    # PLOT
    f = plt.figure(1, figsize=[10, 8])
    gs0 = gridspec.GridSpec(nrows=1, ncols=1, figure=f,
                            hspace=0.08, right=0.91, left=0.09, bottom=0.08, top=0.96,)

    ax = f.add_subplot(gs0[0])  # Ray tracing
    # bwr, seismic, coolwarm, RdBu
    m = ax.pcolormesh(XX, ZZ, Gp.T, cmap=cm, vmin=-c / CL, vmax=c / CL,
                      zorder=1, shading="auto")
    add_colorbar(ax, m)
    ax.scatter(sta["XSTA"].values, sta["ZSTA"].values,
               markersize, facecolors="grey", edgecolors="k",
               marker="v", lw=0.95, zorder=3, clip_on=False,
               label="Seismic stations", )

    ax.set_aspect("equal")
    ax.set_xlabel("x [km]", fontsize=fontsize)
    ax.set_ylabel("z [km]", fontsize=fontsize)

    majorLocator = MultipleLocator(10)
    minorLocator = MultipleLocator(2.5)
    ax.xaxis.set_major_locator(majorLocator)
    ax.xaxis.set_minor_locator(minorLocator)
    # ax.set_xticks(np.arange(0, 140, 10))

    majorLocator = MultipleLocator(10)
    minorLocator = MultipleLocator(2.5)
    ax.yaxis.set_major_locator(majorLocator)
    ax.yaxis.set_minor_locator(minorLocator)
    ax.set_yticks(np.arange(10, 140, 10))
    ax.set_ylim([100, 0])

    ax.tick_params(axis="both", which="major", labelsize=fontsize)
    ax.tick_params(axis="both", which="minor", labelsize=fontsize)
    if filename:
        f.savefig(filename, dpi=200)

    plt.show()
    plt.close()

    return




