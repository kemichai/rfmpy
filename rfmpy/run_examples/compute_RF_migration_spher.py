"""
Code for calculating 3D receiver function migrations...
-
=============================================
Requirements:
    * obspy
    * scipy
    * pandas
    * ...
=============================================

Location: Chavannes-pres-renens, CH
Date: Mar 2022
Author: Konstantinos Michailos
"""

# for picking moho deps
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import matplotlib.pyplot as plt


import rfmpy.core.migration_sphr as rf_mig
import rfmpy.utils.migration_plots_spher as plot_migration_sphr
import numpy as np
import platform
import os
import matplotlib.pyplot as plt


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
path = work_dir + "/data/RF/RF/"
# path='/media/kmichailos/SEISMIC_DATA/RF_calculations/RF/'
#################
# Read stations #
#################
# Read station coordinates from the rfs (sac files) in a pandas dataframe
sta = rf_mig.read_stations_from_sac(path2rfs=path)

plt.scatter(sta["LONSTA"], sta["LATSTA"],
            c='r', marker='v', edgecolor='k', s=100)
plt.show()



################
# Read RFs     #
################
stream = rf_mig.read_traces_sphr(path2rfs=path, sta=sta)

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
maxx = 25.0
pasx = 0.5

miny = 40.0
maxy = 55.0
pasy = 0.5

minz = -5
# maxz needs to be >= zmax
maxz = 100
pasz = 0.5
# Pass all the migration parameters in a dictionary to use them in functions called below
m_params = {'minx': minx, 'maxx': maxx,
            'pasx': pasx, 'pasy': pasy, 'miny': miny, 'maxy': maxy,
            'minz': minz, 'maxz': maxz, 'pasz': pasz, 'inc': inc, 'zmax': zmax}

################
# Ray tracing  #
################
stream_ray_trace = rf_mig.tracing_3D_sphr(stream=stream, migration_param_dict=m_params,
                                          zmoho=50)
# Plot ray tracing...
plot_migration_sphr.plot_ray_tracing(stream_ray_trace)

piercing_lon = []
piercing_lat = []
for i, tr in enumerate(stream_ray_trace):
    tr.stats.station
    for j, z in enumerate(tr.Z):
        if z > 34 and z < 35:
            # print(tr.Xp[j], tr.Yp[j])
            piercing_lon.append(tr.Xp[j])
            piercing_lat.append(tr.Yp[j])
        # elif z > 49 and z < 51:
        #     # print(tr.Xp[j], tr.Yp[j])
        #     piercing_lon.append(tr.Xp[j])
        #     piercing_lat.append(tr.Yp[j])

plt.scatter(piercing_lon, piercing_lat, alpha=.3,
            c='gray', marker='x', edgecolor='gray', s=50)
plt.scatter(sta["LONSTA"], sta["LATSTA"],
            c='r', marker='v', edgecolor='k', s=100)
plt.show()

wav_p_lon = []
wav_p_lat = []
wav_p_dep = []
for i, tr in enumerate(stream_ray_trace):
    tr.stats.station
    for j, z in enumerate(tr.Z):
            # print(tr.Xp[j], tr.Yp[j])
            wav_p_lon.append(tr.Xp[j])
            wav_p_lat.append(tr.Yp[j])
            wav_p_dep.append(z)



plt.scatter(wav_p_lon, wav_p_lat, alpha=0.5,
            c=wav_p_dep, marker='.', edgecolor=None, s=1)
plt.scatter(sta["LONSTA"], sta["LATSTA"],
            c='r', marker='v', edgecolor='k', s=100)
plt.show()

################
# Migration    #
################
mObs = rf_mig.ccpm_3d(stream_ray_trace, m_params, output_file="all_G3", phase="PS")



# 3D to 2D
profile_A = np.array([[8, 46.5], [8, 47.7]])
G2_, sta, xx, zz = plot_migration_sphr.create_2d_profile(mObs, m_params, profile_A, sta, swath=4, plot=True)

################
# Smoothing    #
################
G2 = rf_mig.ccp_smooth(G2_, m_params)
G2[np.abs(G2) < np.max(np.abs(G2)) * 15 / 100] = 0
G2 = rf_mig.ccpFilter(G2)

# ################
# # Plotting     #
# ################
# plot_migration_sphr.plot_migration_profile(Gp=G2, xx=xx, zz=zz, migration_param_dict=m_params, sta=sta,
#                                       work_directory=work_dir, filename=False)

######################################################################################
######################################################################################
# WIP
# for picking moho deps

from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
fontsize = 12
markersize = 100

# for picking moho deps
migration_param_dict=m_params
work_directory=work_dir
filename=False
XX, ZZ = np.meshgrid(xx, zz)
Gp=G2

pal_col = work_directory + "/data/colormaps/vik.txt"
pal_col = pd.read_csv(pal_col, header=None, index_col=False, sep="\s+", names=["R", "G", "B"])
cm = LinearSegmentedColormap.from_list("blue2red", pal_col.values, len(pal_col))
c = np.min([np.max(Gp), 0.1])
c = 0.06
CL = 2


plt.close('all')
# PLOT
f = plt.figure(1, figsize=[15, 8])
gs0 = gridspec.GridSpec(nrows=1, ncols=1, figure=f,
                        hspace=0.08, right=0.91, left=0.09, bottom=0.08, top=0.96,)
ax = f.add_subplot(gs0[0])  # Ray tracing
ax.scatter(XX, ZZ, c=Gp.T, cmap=cm, s=50, vmin=-c / CL, vmax=c / CL, alpha=.5,
                  zorder=1, picker=True)

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
ax.set_yticks(np.arange(10, zz[-1], 10))
ax.set_ylim([50, 0])

ax.tick_params(axis="both", which="major", labelsize=fontsize)
ax.tick_params(axis="both", which="minor", labelsize=fontsize)

def onpick(event):
    index = event.ind
    index = index[0]
    xy = event.artist.get_offsets()
    print('Dist:', xy[index][0], 'Moho:', xy[index][1])
    # TODO: print lon or lat and knowing the profile ...

f.canvas.mpl_connect('pick_event', onpick)
plt.show()
######################################################################################
######################################################################################
######################################################################################
######################################################################################







####
# Create files for gmt
for i, lon in enumerate(piercing_lon):
    with open('/home/kmichailos/Desktop/codes/github/rfmpy/rfmpy/visualisation/gmt/maps/files/pp.txt', 'a') as f:
            f.write('{} {} \n'.format(lon, piercing_lat[i]))

for i, lon in enumerate(wav_p_lon):
    with open('/home/kmichailos/Desktop/codes/github/rfmpy/rfmpy/visualisation/gmt/maps/files/rp.txt', 'a') as f:
            f.write('{} {} {} \n'.format(lon, wav_p_lat[i], wav_p_dep[i]))