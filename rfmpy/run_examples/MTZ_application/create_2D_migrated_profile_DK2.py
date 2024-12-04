"""
Location: Chavannes-pres-renens, CH
Date: April 2022
Author: Konstantinos Michailos
"""

#import rfmpy.core.migration_sphr as rf_mig
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

import numpy as np
from obspy.taup import TauPyModel
from scipy.interpolate import RegularGridInterpolator
from scipy import interpolate
import obspy
import glob
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from scipy import signal
from scipy.interpolate import interp1d
from obspy.geodetics.base import gps2dist_azimuth as gps2dist
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
    #pal_col = work_directory + "/data/colormaps/vik.txt"
    pal_col = '/home/kalmard/rfmpy/data/colormaps/vik.txt'
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

    majorLocator = MultipleLocator(50)
    minorLocator = MultipleLocator(5)
    ax.xaxis.set_major_locator(majorLocator)
    ax.xaxis.set_minor_locator(minorLocator)
    # ax.set_xticks(np.arange(0, 140, 2))

    majorLocator = MultipleLocator(50)
    minorLocator = MultipleLocator(5)
    ax.yaxis.set_major_locator(majorLocator)
    ax.yaxis.set_minor_locator(minorLocator)
    ax.set_yticks(np.arange(50, zz[-5], 50))
    ax.set_ylim([950, 0])

    ax.tick_params(axis="both", which="major", labelsize=fontsize)
    ax.tick_params(axis="both", which="minor", labelsize=fontsize)
    if plot_title:
        ax.set_title(plot_title)

    # plt.tight_layout()
    if filename:
        f.savefig(filename, dpi=300)
    plt.show()
    plt.close()

    return

def ccp_smooth(G2, migration_param_dict):
    """
    Apply smoothing to the migrated RF image.

    :type G2:
    :param G2:
    :type migration_param_dict: dict
    :param migration_param_dict: Dictionary of grid points for the migration.

    :returns: Smoothed amplitudes in 2D grid.
    """
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
    # NOTE: l0, dl are parameters to decide from
    #       which depth the smoothing becomes important and how much

    zbegin_lisse = -2
    # pasx is in degrees so we modify the line below
    #l0 = 0.05
    l0 = 1./111.11
    #dl = 1000000
    # dl = 100

    with np.errstate(divide="warn"):
        G3 = G2
        for iz in range(G2.shape[1]):
            if zz[iz] < zbegin_lisse:
                G3[:, iz] = G2[:, iz]
            else:
                #sigmal = (zz[iz] / dl + l0) / (pasx *111)
                sigmal = (l0) / pasx
                print(sigmal)

                nbml = G2.shape[0]

                mm = int(np.round(nbml / 2))
                C = (
                    1
                    / (sigmal * np.sqrt(2 * np.pi))
                    * np.exp(-0.5 * np.square(np.arange(nbml) - mm) / (sigmal * sigmal))
                )
                C = C / np.sum(C)
                temp = np.convolve(G2[:, iz], C)
                G3[:, iz] = temp[mm : mm + G2.shape[0]]

    return G3


def ccpFilter(G2):
    """
    Convolution with a Gaussian bell for local smooth.

    :type G2:
    :param G2:

    :returns:
    """

    nbm = 5
    b, a = 3, 1.5
    sigma = 1.0
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

def read_stations_from_sac(path2rfs):
    """
    Read all SAC files to get the seismic site details.

    :type path2rfs: str
    :param path2rfs: Path to the stored RF SAC files.

    :return: Pandas DataFrame of the seismic site details.
    """
    import glob
    import obspy
    import pandas as pd

    all_rfs = glob.glob(path2rfs + '*.SAC')
    sta_names = []
    sta_lats = []
    sta_lons = []
    sta_eles = []
    for rf in all_rfs:
        trace = obspy.read(rf)
        # This way we only append items for unique stations
        if trace[0].stats.sac.kstnm not in sta_names:
            sta_names.append(trace[0].stats.sac.kstnm)
            sta_lats.append(trace[0].stats.sac.stla)
            sta_lons.append(trace[0].stats.sac.stlo)
            sta_eles.append(trace[0].stats.sac.stel)

    d_ = {}
    d_['NAMESTA'] = sta_names
    d_['LATSTA'] = sta_lats
    d_['LONSTA'] = sta_lons
    d_['ALTSTA'] = sta_eles

    sta = pd.DataFrame(d_)
    # Add elevations
    sta["ZSTA"] = (-1) * sta["ALTSTA"].values / 1000
    print(sta)

    return sta



# Set up paths
if platform.node().startswith('kalmard-laptop'):
    data_root_dir = '/media/kalmard/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/home/kalmard/Desktop/codes/github'
    mig_file_dir = '/home/kalmard/Post_doc/Prf_mtz/Python_3Dmigration'
    hard_drive_dir = '/media/kalmard/SEISMIC_DATA/'
else:
    data_root_dir = '/media/kalmard/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/github'
    mig_file_dir = '/home/kalmard/Post_doc/Prf_mtz/Python_3Dmigration'
    hard_drive_dir = '/media/kalmard/SEISMIC_DATA/'

# Define paths
work_dir = os.getcwd()
# Example RFs from a couple of teleseismic events
# path = work_dir + "/data/RF/RF/"
# Path to RFs in the hard drive
# path='/media/kmichailos/SEISMIC_DATA/RF_calculations/RF/'
# Path to RFs in the Desktop
# path = desktop_dir + "/all_rfs/RF/"
path = '/home/kalmard/Post_doc/Prf_mtz/matlab_deconvolution/test2/'

## Define MIGRATION parameters
# Ray-tracing parameters
inc = 2  # km
zmax = 800 # km
# Determine study area (x -> perpendicular to the profile)
minx = -13.0 # degrees 2.optional:2 or -4
maxx = 46.0 # degrees 2.optional:30 or 38
pasx = 0.26 # degrees oldest 0.38
miny = 30.0 # degrees 2.optional:41 or 38
maxy = 64.0 # degrees 2.optional:51 or 54
pasy = 0.18 # degrees oldest 0.27
minz = -5 # km
# maxz needs to be >= zmax
maxz = 800 # km
pasz = 2 # km
# Pass all the migration parameters in a dictionary to use them in functions called below
m_params = {'minx': minx, 'maxx': maxx,
            'pasx': pasx, 'pasy': pasy, 'miny': miny, 'maxy': maxy,
            'minz': minz, 'maxz': maxz, 'pasz': pasz, 'inc': inc, 'zmax': zmax}
################
#############################################
# read stations
#sta = rf_mig.read_stations_from_sac(path2rfs=path)
sta = read_stations_from_sac(path2rfs=path)
# Read the 3D numpy array of the RF amplitudes
# with open('/home/kmichailos/Desktop/All_EPcrust.npy', 'rb') as f:
#     mObs_ep = np.load(f)
with open('/home/kalmard/Post_doc/Prf_mtz/Python_3Dmigration/All_zhu_model_800.npy', 'rb') as f:
    mObs = np.load(f)

#Intresting sections to the paper
##EU04
#profile_A = np.array([[2, 43.8], [18.2, 51]])
#prof_name = 'Cross-section_EU04'
##EU14
#profile_A = np.array([[2, 46.2], [30, 46.2]])
#prof_name = 'Cross-section_EU14'
##EU17
#profile_A = np.array([[8, 50.5], [30, 45.2]])
#prof_name = 'Cross-section_EU17'
##EU30
#profile_A = np.array([[5.7, 50], [18, 44]])
#prof_name = 'Cross-section_EU30'
##Hetényi et al., 2009 a; Direction: N-S
#profile_A = np.array([[19.9, 51], [19.6, 44]])
#prof_name = 'Cross-section_hetenyi_2009a'
##Hetényi et al., 2009 b; Direction: N-S
#profile_A = np.array([[20.8, 51], [20.4, 44]])
#prof_name = 'Cross-section_hetenyi_2009b'
##Hetényi et al., 2009 c; Direction: N-S
#profile_A = np.array([[21.7, 51], [21.2, 44]])
#prof_name = 'Cross-section_hetenyi_2009c'
##Hetényi et al., 2009 d; Direction: W-E
#profile_A = np.array([[12, 47.3], [24.5, 47.3]])
#prof_name = 'Cross-section_hetenyi_2009d'
##Hetényi et al., 2009 e; Direction: W-E
#profile_A = np.array([[12, 48.2], [24.5, 48.2]])
#prof_name = 'Cross-section_hetenyi_2009e'
##Hetényi et al., 2009 f; Direction: W-E
#profile_A = np.array([[12, 49.1], [24.5, 49.1]])
#prof_name = 'Cross-section_hetenyi_2009f'
##Kind et al., 2017 Only ~410 dis 6a; Direction: N-S
#profile_A = np.array([[13, 58], [13, 42]])
#prof_name = 'Cross-section_kind_2017_6a'
##Kind et al., 2017 Only ~410 dis 6b; Direction: N-S
#profile_A = np.array([[16, 58], [16, 42]])
#prof_name = 'Cross-section_kind_2017_6b'
##Kind et al., 2017 Only ~410 dis 6c; Direction: N-S
#profile_A = np.array([[19, 58], [19, 42]])
#prof_name = 'Cross-section_kind_2017_6c'
##Kind et al., 2017 Only ~410 dis 6d; Direction: N-S
#profile_A = np.array([[22, 58], [22, 42]])
#prof_name = 'Cross-section_kind_2017_6d'
##Kind et al., 2017 Only ~410 dis 7a; Direction: N-S
#profile_A = np.array([[24, 58], [24, 42]])
#prof_name = 'Cross-section_kind_2017_7a'
##Kind et al., 2017 Only ~410 dis 7b; Direction: N-S
#profile_A = np.array([[26, 58], [26, 42]])
#prof_name = 'Cross-section_kind_2017_7b'
##Kind et al., 2017 Only ~410 dis 8a; Direction: W-E
#profile_A = np.array([[2, 52], [27, 52]])
#prof_name = 'Cross-section_kind_2017_8a'
##Kind et al., 2017 Only ~410 dis 8b; Direction: W-E
#profile_A = np.array([[2, 49], [27, 49]])
#prof_name = 'Cross-section_kind_2017_8b'
##Kind et al., 2017 Only ~410 dis 8c; Direction: W-E
#profile_A = np.array([[2, 46], [27, 46]])
#prof_name = 'Cross-section_kind_2017_8c'
##Kind et al., 2017 Only ~410 dis 8d; Direction: W-E
#profile_A = np.array([[2, 44], [27, 44]])
#prof_name = 'Cross-section_kind_2017_8d'
##Lombardi et al., 2009 a-c; Direction: W-E
#profile_A = np.array([[0, 49], [20, 49]])
#prof_name = 'Cross-section_lombardi_2009_ac'
##Lombardi et al., 2009 d-f; Direction: W-E
#profile_A = np.array([[0, 45], [20, 45]])
#prof_name = 'Cross-section_lombardi_2009_df'
##Lombardi et al., 2009 g-i; Direction: N-S
#profile_A = np.array([[9, 52], [9, 40]])
#prof_name = 'Cross-section_lombardi_2009_gi'
##Lombardi et al., 2009 j-l; Direction: N-S
#profile_A = np.array([[14, 52], [14, 40]])
#prof_name = 'Cross-section_lombardi_2009_jl'
##Lu et al., 2023 Fig.3a; Direction: N-S
#profile_A = np.array([[5, 48], [16, 40]])
#prof_name = 'Cross-section_lu_2023_3a'
##Laura 1; Direction: 
#profile_A = np.array([[15, 50], [25, 43]])
#prof_name = 'Cross-section_laura1'
##Laura 2; Direction:
#profile_A = np.array([[20, 50], [25, 43]])
#prof_name = 'Cross-section_laura2'
##Laura 3; Direction:
#profile_A = np.array([[23, 50.5], [16, 44]])
#prof_name = 'Cross-section_laura3'
##Laura 4; Direction:
profile_A = np.array([[8, 43], [26, 48]])
prof_name = 'Cross-section_A_and_A'

#Non paralell cross-sectinos
#EU01
#profile_A = np.array([[2, 45.9], [12.4, 50.3]])
#prof_name = 'Cross-section_EU01'
##EU02
#profile_A = np.array([[2, 45.5], [15, 51]])
#prof_name = 'Cross-section_EU02'
##EU03
#profile_A = np.array([[2, 44.5], [16, 51]])
#prof_name = 'Cross-section_EU03'
##EU04
#profile_A = np.array([[2, 43.8], [18.2, 51]])
#prof_name = 'Cross-section_EU04'
##EU04 test
#profile_A = np.array([[5, 44.8], [10, 47.9]])
#prof_name = 'Cross-section_EU04_test'
##EU05
#profile_A = np.array([[2, 42.9], [20, 51]])
#prof_name = 'Cross-section_EU05'
##EU06
#profile_A = np.array([[2, 42], [22, 51]])
#prof_name = 'Cross-section_EU06'
##EU07
#profile_A = np.array([[4, 42], [23.4, 51]])
#prof_name = 'Cross-section_EU07'
##EU08
#profile_A = np.array([[7, 42.5], [24, 51]])
#prof_name = 'Cross-section_EU08'
##EU09
#profile_A = np.array([[8, 42.5], [27.5, 51]])
#prof_name = 'Cross-section_EU09'
##EU10
#profile_A = np.array([[8, 41.9], [25, 49]])
#prof_name = 'Cross-section_EU10'
##EU11
#profile_A = np.array([[2, 48.2], [30, 48.2]])
#prof_name = 'Cross-section_EU11'
##EU12
#profile_A = np.array([[2, 47.8], [30, 47.8]])
#prof_name = 'Cross-section_EU12'
##EU13
#profile_A = np.array([[2, 47.2], [30, 47.2]])
#prof_name = 'Cross-section_EU13'
##EU14
#profile_A = np.array([[2, 46.2], [30, 46.2]])
#prof_name = 'Cross-section_EU14'
##EU15
#profile_A = np.array([[2, 46], [30, 46]])
#prof_name = 'Cross-section_EU15'
##EU16
#profile_A = np.array([[2, 45.3], [30, 45.3]])
#prof_name = 'Cross-section_EU16'
##EU17
#profile_A = np.array([[8, 50.5], [30, 45.2]])
#prof_name = 'Cross-section_EU17'
##EU18
#profile_A = np.array([[6, 50.5], [30, 44.5]])
#prof_name = 'Cross-section_EU18'
##EU19
#profile_A = np.array([[6, 50], [30, 43.9]])
#prof_name = 'Cross-section_EU19'
##EU20
#profile_A = np.array([[17, 51], [30, 45.1]])
#prof_name = 'Cross-section_EU20'
##EU21
#profile_A = np.array([[15.8, 51], [30, 44.3]])
#prof_name = 'Cross-section_EU21'
##EU22
#profile_A = np.array([[14.4, 51], [30, 44]])
#prof_name = 'Cross-section_EU22'
##EU23
#profile_A = np.array([[13, 51], [30, 43.6]])
#prof_name = 'Cross-section_EU23'
##EU24
#profile_A = np.array([[11.7, 51], [30, 43.3]])
#prof_name = 'Cross-section_EU24'
##EU25
#profile_A = np.array([[10.3, 51], [28, 43.5]])
#prof_name = 'Cross-section_EU25'
##EU26
#profile_A = np.array([[10, 50.2], [27, 43.4]])
#prof_name = 'Cross-section_EU26'
##EU27
#profile_A = np.array([[8, 50.5], [26, 42.5]])
#prof_name = 'Cross-section_EU27'
##EU28
#profile_A = np.array([[7, 50.5], [24, 43]])
#prof_name = 'Cross-section_EU28'
##EU29
#profile_A = np.array([[6, 50.2], [20, 44]])
#prof_name = 'Cross-section_EU29'
##EU30
#profile_A = np.array([[5.7, 50], [18, 44]])
#prof_name = 'Cross-section_EU30'
##EU31
#profile_A = np.array([[5, 49.5], [16, 44.1]])
#prof_name = 'Cross-section_EU31'
##EU32
#profile_A = np.array([[4, 49], [14, 45]])
#prof_name = 'Cross-section_EU32'
##EU33
#profile_A = np.array([[4, 48.5], [14, 44.3]])
#prof_name = 'Cross-section_EU33'
##EU34
#profile_A = np.array([[3, 48.3], [14, 43.8]])
#prof_name = 'Cross-section_EU34'
##EU35
#profile_A = np.array([[3.2, 48], [13, 43.5]])
#prof_name = 'Cross-section_EU35'
##EU36
#profile_A = np.array([[2, 47.8], [12, 43.6]])
#prof_name = 'Cross-section_EU36'
##EU37
#profile_A = np.array([[2, 47.5], [12, 43.7]])
#prof_name = 'Cross-section_EU37'
##EU38
#profile_A = np.array([[2, 46.8], [12, 42.8]])
#prof_name = 'Cross-section_EU38'
##EU39
#profile_A = np.array([[2, 46.2], [11.7, 42]])
#prof_name = 'Cross-section_EU39'
##EU40
#profile_A = np.array([[2, 46], [10.4, 42]])
#prof_name = 'Cross-section_EU40'
##EU41
#profile_A = np.array([[2, 45.4], [9, 42]])
#prof_name = 'Cross-section_EU41'
###################################################
###################################################
##Paralell cross-sectinos
#EU01
#profile_A = np.array([[2, 51], [2, 41.5]])
#prof_name = 'Cross-section_EU01'
##EU02
#profile_A = np.array([[3, 51], [3, 41.5]])
#prof_name = 'Cross-section_EU02'
##EU03
#profile_A = np.array([[4, 51], [4, 41.5]])
#prof_name = 'Cross-section_EU03'
##EU04
#profile_A = np.array([[5, 51], [5, 41.5]])
#prof_name = 'Cross-section_EU04'
##EU05
#profile_A = np.array([[6, 51, [6, 41.5]])
#prof_name = 'Cross-section_EU05'
##EU06
#profile_A = np.array([[7, 51], [7, 41.5]])
#prof_name = 'Cross-section_EU06'
##EU07
#profile_A = np.array([[8, 51], [8, 41.5]])
#prof_name = 'Cross-section_EU07'
##EU08
#profile_A = np.array([[9, 51], [9, 41.5]])
#prof_name = 'Cross-section_EU08'
##EU09
#profile_A = np.array([[10, 51], [10, 41.5]])
#prof_name = 'Cross-section_EU09'
##EU10
#profile_A = np.array([[11, 51], [11, 41.5]])
#prof_name = 'Cross-section_EU10'
##EU11
#profile_A = np.array([[12, 51], [12, 41.5]])
#prof_name = 'Cross-section_EU11'
##EU12
#profile_A = np.array([[13, 51], [13, 41.5]])
#prof_name = 'Cross-section_EU12'
##EU13
#profile_A = np.array([[14, 51], [14, 41.5]])
#prof_name = 'Cross-section_EU13'
##EU14
#profile_A = np.array([[15, 51], [15, 41.5]])
#prof_name = 'Cross-section_EU14'
##EU15
#profile_A = np.array([[16, 51], [16, 41.5]])
#prof_name = 'Cross-section_EU15'
##EU16
#profile_A = np.array([[17, 51], [17, 41.5]])
#prof_name = 'Cross-section_EU16'
##EU17
#profile_A = np.array([[18, 51], [18, 41.5]])
#prof_name = 'Cross-section_EU17'
##EU18
#profile_A = np.array([[19, 51], [19, 41.5]])
#prof_name = 'Cross-section_EU18'
##EU19
#profile_A = np.array([[20, 51], [20, 41.5]])
#prof_name = 'Cross-section_EU19'
##EU20
#profile_A = np.array([[21, 51], [21, 41.5]])
#prof_name = 'Cross-section_EU20'
##EU21
#profile_A = np.array([[22, 51], [22, 41.5]])
#prof_name = 'Cross-section_EU21'
##EU22
#profile_A = np.array([[23, 51], [23, 41.5]])
#prof_name = 'Cross-section_EU22'
##EU23
#profile_A = np.array([[24, 51], [24, 41.5]])
#prof_name = 'Cross-section_EU23'
##EU24
#profile_A = np.array([[25, 51], [25, 41.5]])
#prof_name = 'Cross-section_EU24'
##EU25
#profile_A = np.array([[26, 51], [26, 41.5]])
#prof_name = 'Cross-section_EU25'
##EU26
#profile_A = np.array([[27, 51], [27, 41.5]])
#prof_name = 'Cross-section_EU26'
##EU27
#profile_A = np.array([[28, 51], [28, 41.5]])
#prof_name = 'Cross-section_EU27'
##EU28
#profile_A = np.array([[29, 51], [29, 41.5]])
#prof_name = 'Cross-section_EU28'
##EU29
#profile_A = np.array([[30, 51], [30, 41.5]])
#prof_name = 'Cross-section_EU29'
##EU30
#profile_A = np.array([[2, 51], [30, 51]])
#prof_name = 'Cross-section_EU30'
##EU31
#profile_A = np.array([[2, 50.5], [30, 50.5]])
#prof_name = 'Cross-section_EU31'
##EU32
#profile_A = np.array([[2, 50], [30, 50]])
#prof_name = 'Cross-section_EU32'
##EU33
#profile_A = np.array([[2, 49.5], [30, 49.5]])
#prof_name = 'Cross-section_EU33'
##EU34
#profile_A = np.array([[2, 49], [30, 49]])
#prof_name = 'Cross-section_EU34'
##EU35
#profile_A = np.array([[2, 48.5], [30, 48.5]])
#prof_name = 'Cross-section_EU35'
##EU36
#profile_A = np.array([[2, 48], [30, 48]])
#prof_name = 'Cross-section_EU36'
##EU37
#profile_A = np.array([[2, 47.5], [30, 47.5]])
#prof_name = 'Cross-section_EU37'
##EU38
#profile_A = np.array([[2, 47, [30, 47]])
#prof_name = 'Cross-section_EU38'
##EU39
#profile_A = np.array([[2, 46.5], [30, 46.5]])
#prof_name = 'Cross-section_EU39'
##EU40
#profile_A = np.array([[2, 46], [30, 46]])
#prof_name = 'Cross-section_EU40'
##EU41
#profile_A = np.array([[2, 45.5], [30, 45.5]])
#prof_name = 'Cross-section_EU41'
##EU42
#profile_A = np.array([[2, 45], [30, 45]])
#prof_name = 'Cross-section_EU42'
##EU43
#profile_A = np.array([[2, 44.5], [30, 44.5]])
#prof_name = 'Cross-section_EU43'
##EU44
#profile_A = np.array([[2, 44], [30, 44]])
#prof_name = 'Cross-section_EU44'
##EU45
#profile_A = np.array([[2, 43.5], [30, 43.5]])
#prof_name = 'Cross-section_EU45'
##EU46
#profile_A = np.array([[2, 43], [30, 43]])
#prof_name = 'Cross-section_EU46'
##EU47
#profile_A = np.array([[2, 42.5], [30, 42.5]])
#prof_name = 'Cross-section_EU47'
##EU48
#profile_A = np.array([[2, 42], [30, 42]])
#prof_name = 'Cross-section_EU48'
##EU49
#profile_A = np.array([[2, 41.5], [30, 41.5]])
#prof_name = 'Cross-section_EU49'


G2_, sta_, xx, zz = plot_migration_sphr.create_2d_profile(mObs, m_params, profile_A, sta, swath=300, plot=True)
################
# Smoothing    #
################
#mObs = rf_mig.ccp_smooth(G2_, m_params)
mObs = ccp_smooth(G2_, m_params)
#mObs[np.abs(mObs) < np.max(np.abs(mObs)) * 15 / 100] = 0
#mObs = rf_mig.ccpFilter(mObs)
mObs = ccpFilter(mObs)
# ################
# # Plotting     #
# ################
plot_migration_profile(Gp=mObs, xx=xx, zz=zz, migration_param_dict=m_params, sta=sta_,
                                           work_directory=work_dir, filename='test', plot_title='test')
######################################################################################
######################################################################################
# File for creating cross-sections with GMT
for i, x in enumerate(xx):
    for j, z in enumerate(zz):
        print(kilometers2degrees(x), z, mObs[i,j])
        with open('/home/kalmard/Post_doc/Prf_mtz/Python_3Dmigration/' + prof_name + '.txt', 'a') as of:
            of.write('{} {} {} \n'.
                     format(kilometers2degrees(x), z, mObs[i, j]))
######################################################################################
######################################################################################
