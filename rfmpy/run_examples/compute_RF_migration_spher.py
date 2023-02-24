"""
Code for calculating 3D receiver function migrations...
-
=============================================
Requirements:
    * obspy
    * scipy
    * pandas

=============================================

Location: Chavannes-pres-renens, CH
Date: Mar 2022
Author: Konstantinos Michailos
"""
import matplotlib
matplotlib.use('TkAgg')
import rfmpy.core.migration_sphr as rf_mig
import rfmpy.utils.migration_plots_spher as plot_migration_sphr
import numpy as np
import platform
import os
import matplotlib.pyplot as plt
import time

t_beg = time.time()
# Set up paths
if platform.node().startswith('kmichailos-laptop'):
    data_root_dir = '/media/kmichailos/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/home/kmichailos/Desktop/codes/github'
    desktop_dir = '/home/kmichailos/Desktop'
    hard_drive_dir = '/media/kmichailos/SEISMIC_DATA/'
else:
    data_root_dir = '/media/konstantinos/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/home/Desktop/codes'
    desktop_dir = '/home/konstantinos/Desktop'
    hard_drive_dir = '/media/konstantinos/SEISMIC_DATA/'

# Define paths
work_dir = os.getcwd()
# Example RFs from a couple of teleseismic events
# path = work_dir + "/data/RF/RF/"
# Path to RFs in the hard drive
# path='/media/kmichailos/SEISMIC_DATA/RF_calculations/RF/'
# Path to RFs in the Desktop
path = desktop_dir + "/all_rfs/RF/"
#################
# Read stations #
#################
# Read station coordinates from the rfs (sac files) in a pandas dataframe
sta = rf_mig.read_stations_from_sac(path2rfs=path)

################
# Read RFs     #
################
stream = rf_mig.read_traces_sphr(path2rfs=path, sta=sta)

# Define MIGRATION parameters
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
################
# Ray tracing  #
################
# Pick one of the two velocity models
# 'EPcrust' or 'iasp91'
stream_ray_trace = rf_mig.tracing_3D_sphr(stream=stream, migration_param_dict=m_params,
                                          velocity_model='EPcrust')
# Write piercing points in a file
plot_migration_sphr.write_files_4_piercing_points_and_raypaths(stream_ray_trace, sta, piercing_depth=35, plot=True)
################
# Migration    #
################
mObs = rf_mig.ccpm_3d(stream_ray_trace, m_params, output_file="/home/konstantinos/Desktop/All_EPcrust_new_mantle_vel", phase="PS")


total_time = time.time() - t_beg
print('Ray tracing took ' + str(round(total_time)/60) + ' minutes in total.')
