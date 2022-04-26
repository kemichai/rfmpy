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

Note: Based on codes originally written by Matteo Scarponi.

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
path = work_dir + "/data/RF/"

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
minx = 5.0
maxx = 15.0
pasx = 0.5

miny = 45.0
maxy = 55.0
pasy = 0.5

minz = -2
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
################
# Migration    #
################
mObs = rf_mig.ccpm_3d(stream_ray_trace, m_params, phase="PS")



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

