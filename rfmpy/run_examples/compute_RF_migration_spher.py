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
pasz = 2
# Pass all the migration parameters in a dictionary to use them in functions called below
m_params = {'minx': minx, 'maxx': maxx,
            'pasx': pasx, 'pasy': pasy, 'miny': miny, 'maxy': maxy,
            'minz': minz, 'maxz': maxz, 'pasz': pasz, 'inc': inc, 'zmax': zmax}

################
# Ray tracing  #
################
stream_ray_trace = rf_mig.tracing_3D_sphr(stream=stream, migration_param_dict=m_params,
                                          zMoho=50)
# Plot ray tracing...
plot_migration_sphr.plot_ray_tracing(stream_ray_trace)

piercing_lon = []
piercing_lat = []
for i, tr in enumerate(stream_ray_trace):
    tr.stats.station
    for j, z in enumerate(tr.Z):
        if z > 29 and z < 31:
            # print(tr.Xp[j], tr.Yp[j])
            piercing_lon.append(tr.Xp[j])
            piercing_lat.append(tr.Yp[j])
        elif z > 49 and z < 51:
            # print(tr.Xp[j], tr.Yp[j])
            piercing_lon.append(tr.Xp[j])
            piercing_lat.append(tr.Yp[j])

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
            print(tr.Xp[j], tr.Yp[j])
            wav_p_lon.append(tr.Xp[j])
            wav_p_lat.append(tr.Yp[j])
            wav_p_dep.append(z)



plt.scatter(wav_p_lon, wav_p_lat, alpha=0.5,
            c=wav_p_dep, marker='.', edgecolor=None, s=100)
plt.scatter(sta["LONSTA"], sta["LATSTA"],
            c='r', marker='v', edgecolor='k', s=100)
plt.show()

################
# Migration    #
################
mObs = rf_mig.ccpm_3d(stream_ray_trace, m_params, output_file="all_G3", phase="PS")



################
# Smoothing    #
################
# mObs = rf_mig.ccp_smooth(mObs, m_params)
# mObs[np.abs(mObs) < np.max(np.abs(mObs)) * 15 / 100] = 0
# mObs = rf_mig.ccpFilter(mObs)
#
# ################
# # Plotting     #
# ################
plot_migration_sphr.plot_migration_profile(Gp=mObs, migration_param_dict=m_params, sta=sta,
                                      work_directory=work_dir, filename=False)

