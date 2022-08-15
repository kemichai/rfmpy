import matplotlib
matplotlib.use('TkAgg')
import rfmpy.core.migration_sphr as rf_mig
import rfmpy.utils.migration_plots_spher as plot_migration_sphr
import numpy as np
import platform
import os
import matplotlib.pyplot as plt
import time
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from collections import OrderedDict
import matplotlib.patches as patches
from scipy.interpolate import RegularGridInterpolator
from pyproj import Geod
import pyproj
from obspy.geodetics import degrees2kilometers, kilometers2degrees



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

# Set up paths
if platform.node().startswith('kmichailos-laptop'):
    data_root_dir = '/media/kmichailos/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/home/kmichailos/Desktop/codes/github'
    desktop_dir = '/home/kmichailos/Desktop'
    hard_drive_dir = '/media/kmichailos/SEISMIC_DATA/'
else:
    data_root_dir = '/media/kmichall/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/github'
    desktop_dir = '/home/kmichall/Desktop'
    hard_drive_dir = '/media/kmichall/SEISMIC_DATA/'

# Define paths
work_dir = os.getcwd()
# path = work_dir + "/data/RF/RF/"
# path = desktop_dir + "/RF_test/RF/"
# path='/media/kmichailos/SEISMIC_DATA/RF_calculations/RF/'
path = desktop_dir + "/all_rfs/RF/"

# Define MIGRATION parameters
# min lat=45.0
# max lat=50.0
# min lon= 5.0
# max lon = 10.0
# Ray-tracing parameters
inc = 0.25
zmax = 100
# Determine study area (x -> perpendicular to the profile)
minx = 0.0
maxx = 30.0
pasx = 0.05

miny = 30.0
maxy = 60.0
pasy = 0.05

minz = -5
# maxz needs to be >= zmax
maxz = 100
pasz = 0.5
# Pass all the migration parameters in a dictionary to use them in functions called below
m_params = {'minx': minx, 'maxx': maxx,
            'pasx': pasx, 'pasy': pasy, 'miny': miny, 'maxy': maxy,
            'minz': minz, 'maxz': maxz, 'pasz': pasz, 'inc': inc, 'zmax': zmax}

# read stations
sta = rf_mig.read_stations_from_sac(path2rfs=path)
# COPY PASTE FROM HERE |
#                     /|\


with open('/home/kmichailos/Desktop/All_EPcrust.npy', 'rb') as f:
    mObs_ep = np.load(f)

with open('/home/kmichailos/Desktop/All_iasp91.npy', 'rb') as f:
    mObs_ia = np.load(f)



# 3D to 2D
# NOTE
# profile_A = np.array([[8, 46], [8, 48]])
profile_A = np.array([[13.35, 50.6], [13.35, 45.6]])
profile_A = np.array([[3., 43.5], [19., 48]])

# test for picking
profile_A = np.array([[10., 40], [10., 50]])
# profile_A = np.array([[5., 43], [5., 50]])
# profile_A = np.array([[15., 43], [15., 50]])
# profile_A = np.array([[20., 43], [20., 50]])
# #
# # profile_A = np.array([[5., 45], [20., 45]])
# # profile_A = np.array([[5., 47.5], [20., 47.5]])
# # profile_A = np.array([[5., 50], [20., 50]])
# profile_A = np.array([[5., 49], [20., 49]])
#
# profile_A = np.array([[15., 49], [20., 49]])
# profile_A = np.array([[20., 47.5], [25., 47.5]])
# profile_A = np.array([[15., 46], [25., 46]])
# profile_A = np.array([[15., 46.5], [20., 46.5]])


G2_, sta, xx, zz = plot_migration_sphr.create_2d_profile_4_moho_picker(mObs_ep, m_params, profile_A, sta, swath=50, plot=True)


def ccp_smooth(G2, migration_param_dict):
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
    zbegin_lisse = -5
    # pasx is in degrees so we modify the line below
    l0 = 1./111.11
    # dl = 1000000
    # dl = 100
    with np.errstate(divide="warn"):
        G3 = G2
        for iz in range(G2.shape[1]):
            if zz[iz] < zbegin_lisse:
                G3[:, iz] = G2[:, iz]
            else:
                # sigmal = (zz[iz] / dl + l0) / (pasx*111.11)
                sigmal = (l0) / pasx
                # print(sigmal)
                nbml = G2.shape[0]
                mm = int(np.round(nbml / 2))
                C = (1 / (sigmal * np.sqrt(2 * np.pi)) * np.exp(-0.5 * np.square(np.arange(nbml) - mm) / (sigmal * sigmal))           )
                C = C / np.sum(C)
                temp = np.convolve(G2[:, iz], C)
                G3[:, iz] = temp[mm : mm + G2.shape[0]]
    return G3
# The value of sigma roughly
# corresponds to the number of neighbouring cells on either side of a cell over which
# a signal is smoothed. Pick a place on the raw image where the adjacent cells are
# 0, and see the effect when sigma is 0.5, 1 or 2. 2 is too much, I think Matteo used
# it as he had shorter horizontal cells -- so you either shorten as well, or reduce sigma.
# I would say that your sigma=1 is visually most appealing to me, but it does not include
# smoothing yet. I would keep this little filter as harmless as possible, and tune the smoothing function.
#
# Looking at the code again, it seems that nbm, b and a in the ccpFilter function
# relate to number of cells, not true distances.
#
# Ultimately, a single trace's single Moho conversion signal should roughly ' \
#                           'correspond to the Fresnel zone radius at that depth. Let's say
# that is is ca. 10 km, then a signal should be spread over 20 km distance, that is 4 cells in your current setting
# -- we could reduce it a bit.

from scipy import signal

def ccpFilter(G2):
    """
    Convolution with a Gaussian bell for local smooth.

    :type G2:
    :param G2:

    :returns:
    """

    nbm = 5
    b, a =  3, 1.5
    sigma = 1.0
    # sigma = 0.1
    C = np.zeros((nbm, nbm))
    mm = np.floor(nbm / 2)
    for i in range(nbm):
        for j in range(nbm):
            r = np.sqrt(((i - mm) ** 2) / (a ** 2) + ((j - mm) ** 2 / (b ** 2)))
            if r < mm:
                C[i, j] = np.exp(-0.5 * (r / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))
    C = C / np.sum(C)
    miniG2 = signal.convolve2d(G2, C, "same")
    return miniG2


G2 = ccp_smooth(G2_, m_params)
# G2[np.abs(G2) < np.max(np.abs(G2)) * 15 / 100] = 0
G2 = ccpFilter(G2)

# ################
# # Plotting     #
# ################

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

    XX, ZZ = np.meshgrid(xx, zz)

    # COLOR PALETTE AND COLORMAP
    # cork, bam, broc, vik
    pal_col = work_directory + "/data/colormaps/vik.txt"
    pal_col = pd.read_csv(pal_col, header=None, index_col=False, sep="\s+", names=["R", "G", "B"])
    cm = LinearSegmentedColormap.from_list("blue2red", pal_col.values, len(pal_col))
    c = np.min([np.max(Gp), 0.1])
    c = 0.1
    CL = 1

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

    majorLocator = MultipleLocator(50)
    minorLocator = MultipleLocator(10)
    ax.xaxis.set_major_locator(majorLocator)
    ax.xaxis.set_minor_locator(minorLocator)
    # ax.set_xticks(np.arange(0, 140, 10))

    majorLocator = MultipleLocator(10)
    minorLocator = MultipleLocator(5)
    ax.yaxis.set_major_locator(majorLocator)
    ax.yaxis.set_minor_locator(minorLocator)
    ax.set_yticks(np.arange(0, zz[-1], 10))
    ax.set_ylim([100, -10])
    ax.set_xlim([xx[0] - 10 , xx[-1] + 10])

    ax.tick_params(axis="both", which="major", labelsize=fontsize)
    ax.tick_params(axis="both", which="minor", labelsize=fontsize)
    if plot_title:
        ax.set_title(plot_title)

    plt.tight_layout()
    if filename:
        f.savefig(filename, dpi=200)
    plt.show()
    plt.close()

    return


# plot_migration_profile(Gp=G2, xx=xx, zz=zz, migration_param_dict=m_params, sta=sta,
#                        work_directory=work_dir, filename='epcrust',
#                        plot_title='epcrust')

# Manually pick moho deps
# IMPORTANT NOTE: only works with cross-sections the have S-N and E-W directions!!!
plot_migration_sphr.moho_picker(Gp=G2, xx=xx, zz=zz, migration_param_dict=m_params,
                                sta=sta, work_directory=work_dir, profile=profile_A)


# for i, x in enumerate(xx):
#     for j, z in enumerate(zz):
#         print(kilometers2degrees(x), z, G2[i,j])
#         with open('/home/kmichailos/Desktop/codes/github/rfmpy/rfmpy/visualisation/gmt/cross_sections/xyz_smoothed_test.txt', 'a') as of:
#             of.write('{} {} {} \n'.
#                      format(kilometers2degrees(x), z, G2[i, j]))

######################################################################################
######################################################################################

