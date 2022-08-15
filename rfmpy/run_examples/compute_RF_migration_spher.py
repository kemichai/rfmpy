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
    data_root_dir = '/media/kmichall/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/github'
    desktop_dir = '/home/kmichall/Desktop'
    hard_drive_dir = '/media/kmichall/SEISMIC_DATA/'

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
plot_migration_sphr.write_files_4_piercing_points_and_raypaths(stream_ray_trace)
################
# Migration    #
################
mObs = rf_mig.ccpm_3d(stream_ray_trace, m_params, output_file="/home/kmichailos/Desktop/ALL_EPcrust", phase="PS")
total_time = time.time() - t_beg
print('Ray tracing took ' + str(round(total_time)/60) + ' minutes in total.')






# 3D to 2D
profile_A = np.array([[8, 46], [8, 48]])
G2_, sta, xx, zz = plot_migration_sphr.create_2d_profile(mObs, m_params, profile_A, sta, swath=25, plot=True)

################
# Smoothing    #
################
G2 = rf_mig.ccp_smooth(G2_, m_params)
# G2[np.abs(G2) < np.max(np.abs(G2)) * 15 / 100] = 0
G2 = rf_mig.ccpFilter(G2)

# ################
# # Plotting     #
# ################
plot_migration_sphr.plot_migration_profile(Gp=G2, xx=xx, zz=zz, migration_param_dict=m_params, sta=sta,
                                      work_directory=work_dir, filename='EPcrust', plot_title='EPcrust')

######################################################################################
######################################################################################


# Manually pick moho deps
# IMPORTANT NOTE: only works with cross-sections the have S-N and E-W directions!!!
plot_migration_sphr.moho_picker(Gp=G2, xx=xx, zz=zz, migration_param_dict=m_params,
                                sta=sta, work_directory=work_dir, profile=profile_A)
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
