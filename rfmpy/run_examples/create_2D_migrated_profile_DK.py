"""
Location: Chavannes-pres-renens, CH
Date: April 2022
Author: Konstantinos Michailos
"""

import rfmpy.core.migration_sphr as rf_mig
import rfmpy.utils.migration_plots_spher as plot_migration_sphr
import os
import platform
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from obspy.geodetics import degrees2kilometers, kilometers2degrees




def plot_migration_profile(Gp, xx, zz, migration_param_dict, sta, work_directory, filename=False, plot_title=None):
    """

    :param Gp:
    :param xx:
    :param zz:
    :param migration_param_dict:
    :param sta:
    :param work_directory:
    :param filename:
    :return:
    """

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

    XX, ZZ = np.meshgrid(xx, zz)
    fontsize = 12
    markersize = 100
    # zz_corr = []
    # R = 6371.009
    # for i, dist in enumerate(prof_dist):
    #     correction = np.sqrt(dist**2 + R)
    #     print(correction)
    #     zz_corr = zz + correction



    # COLOR PALETTE AND COLORMAP
    # cork, bam, broc, vik
    pal_col = work_directory + "/data/colormaps/vik.txt"
    pal_col = pd.read_csv(pal_col, header=None, index_col=False, sep="\s+", names=["R", "G", "B"])
    cm = LinearSegmentedColormap.from_list("blue2red", pal_col.values, len(pal_col))
    c = np.min([np.max(Gp), 0.1])
    c = 0.15
    CL = 2

    # PLOT
    f = plt.figure(1, figsize=[8,15])
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

    majorLocator = MultipleLocator(20)
    minorLocator = MultipleLocator(5)
    ax.xaxis.set_major_locator(majorLocator)
    ax.xaxis.set_minor_locator(minorLocator)
    # ax.set_xticks(np.arange(0, 140, 2))

    majorLocator = MultipleLocator(20)
    minorLocator = MultipleLocator(5)
    ax.yaxis.set_major_locator(majorLocator)
    ax.yaxis.set_minor_locator(minorLocator)
    ax.set_yticks(np.arange(20, zz[-1], 20))
    ax.set_ylim([800, 0])

    ax.tick_params(axis="both", which="major", labelsize=fontsize)
    ax.tick_params(axis="both", which="minor", labelsize=fontsize)
    if plot_title:
        ax.set_title(plot_title)

    # plt.tight_layout()
    if filename:
        f.savefig(filename, dpi=200)
    plt.show()
    plt.close()

    return



# Set up paths
if platform.node().startswith('kmichailos-laptop'):
    data_root_dir = '/media/kmichailos/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/home/kmichailos/Desktop/codes/github'
    desktop_dir = '/home/kmichailos/Desktop'
    hard_drive_dir = '/media/kmichailos/SEISMIC_DATA/'
else:
    data_root_dir = '/media/konstantinos/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/github'
    desktop_dir = '/home/konstantinos/Desktop'
    hard_drive_dir = '/media/konstantinos/SEISMIC_DATA/'

# Define paths
work_dir = os.getcwd()
# Example RFs from a couple of teleseismic events
# path = work_dir + "/data/RF/RF/"
# Path to RFs in the hard drive
# path='/media/kmichailos/SEISMIC_DATA/RF_calculations/RF/'
# Path to RFs in the Desktop
# path = desktop_dir + "/all_rfs/RF/"
path = desktop_dir + "/RF_test/"

# Define MIGRATION parameters
# Ray-tracing parameters
inc = 10
zmax = 700
# Determine study area (x -> perpendicular to the profile)
minx = 0.0
maxx = 30.0
pasx = 0.5
miny = 30.0
maxy = 60.0
pasy = 0.5
minz = -5
# maxz needs to be >= zmax
maxz = 700
pasz = 15
# Pass all the migration parameters in a dictionary to use them in functions called below
m_params = {'minx': minx, 'maxx': maxx,
            'pasx': pasx, 'pasy': pasy, 'miny': miny, 'maxy': maxy,
            'minz': minz, 'maxz': maxz, 'pasz': pasz, 'inc': inc, 'zmax': zmax}
################
#############################################
# read stations
sta = rf_mig.read_stations_from_sac(path2rfs=path)
# Read the 3D numpy array of the RF amplitudes
# with open('/home/kmichailos/Desktop/All_EPcrust.npy', 'rb') as f:
#     mObs_ep = np.load(f)
with open(desktop_dir + '/All_zmodel_m60_vel.npy', 'rb') as f:
    mObs = np.load(f)

profile_A = np.array([[7.5, 45], [7.5, 48]])
prof_name = 'Cross-section_test'
G2_, sta_, xx, zz = plot_migration_sphr.create_2d_profile(mObs, m_params, profile_A, sta, swath=250, plot=True)
################
# Smoothing    #
################
# mObs = rf_mig.ccp_smooth(G2_, m_params)
# mObs[np.abs(mObs) < np.max(np.abs(mObs)) * 15 / 100] = 0
# ################
# # Plotting     #
# ################
plot_migration_profile(Gp=G2_, xx=xx, zz=zz, migration_param_dict=m_params, sta=sta_,
                                           work_directory=work_dir, filename='DK', plot_title='DK')
######################################################################################
######################################################################################
# File for creating cross-sections with GMT
for i, x in enumerate(xx):
    for j, z in enumerate(zz):
        print(kilometers2degrees(x), z, mObs[i,j])
        with open('/home/kmichailos/Desktop/codes/github/rfmpy/rfmpy/visualisation/gmt/cross_sections/appendix_moho_picks/cc_files/' + prof_name + 'TEST.txt', 'a') as of:
            of.write('{} {} {} \n'.
                     format(kilometers2degrees(x), z, mObs[i, j]))
######################################################################################
######################################################################################