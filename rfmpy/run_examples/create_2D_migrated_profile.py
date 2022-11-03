"""
Location: Chavannes-pres-renens, CH
Date: April 2022
Author: Konstantinos Michailos
"""

import rfmpy.core.migration_sphr as rf_mig
import rfmpy.utils.migration_plots_spher as plot_migration_sphr
import numpy as np
import os
import matplotlib.pyplot as plt
from obspy.geodetics import degrees2kilometers, kilometers2degrees
import platform


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
#############################################
# read stations
sta = rf_mig.read_stations_from_sac(path2rfs=path)
# Read the 3D numpy array of the RF amplitudes
# with open('/home/kmichailos/Desktop/All_EPcrust.npy', 'rb') as f:
#     mObs_ep = np.load(f)
with open('/home/kmichailos/Desktop/All_EPcrust_new_mantle_vel.npy', 'rb') as f:
    mObs_ep = np.load(f)
# with open('/home/kmichailos/Desktop/All_iasp91.npy', 'rb') as f:
#     mObs_ia = np.load(f)

# 3D to 2D
profile_A = np.array([[10., 40], [10., 50]])
# EASI test
profile_A = np.array([[13.35, 50.6], [13.35, 45.6]])
# prof_name = 'EPcrust'
# prof_name = 'iasp91'
# A - A'
profile_A = np.array([[3, 44.1], [9, 44.8]])
prof_name = 'A_A'
# # B - B'
# profile_A = np.array([[6, 49], [11.5, 44]])
# prof_name = 'B_B'
# # C - C'
# profile_A = np.array([[11, 45.5], [22, 50]])
# prof_name = 'C_C'




# # 3D to 2D
# 0
# profile_A = np.array([[3, 43], [3, 50]])
# prof_name = 'Cross-section_0'
# 1
profile_A = np.array([[4, 43], [4, 50]])
prof_name = 'Cross-section_1'
# 2
profile_A = np.array([[5, 43], [5, 50]])
prof_name = 'Cross-section_2'
# 3
profile_A = np.array([[6, 43], [6, 50]])
prof_name = 'Cross-section_3'
# 4
profile_A = np.array([[7, 43], [7, 50]])
prof_name = 'Cross-section_4'
# 5
profile_A = np.array([[8, 43], [8, 50]])
prof_name = 'Cross-section_5'
# 6
profile_A = np.array([[9, 43], [9, 52]])
prof_name = 'Cross-section_6'
# 7
profile_A = np.array([[10, 43], [10, 52]])
prof_name = 'Cross-section_7'
# 8
profile_A = np.array([[11, 43], [11, 52]])
prof_name = 'Cross-section_8'
# 9
profile_A = np.array([[12, 43], [12, 52]])
prof_name = 'Cross-section_9'
# 10
profile_A = np.array([[13, 43], [13, 52]])
prof_name = 'Cross-section_10'
# 11
profile_A = np.array([[14, 43], [14, 52]])
prof_name = 'Cross-section_11'
# 12
profile_A = np.array([[15, 43], [15, 52]])
prof_name = 'Cross-section_12'
# 13
profile_A = np.array([[16, 43], [16, 52]])
prof_name = 'Cross-section_13'
# 14
profile_A = np.array([[17, 43], [17, 52]])
prof_name = 'Cross-section_14'
# 15
profile_A = np.array([[18, 43], [18, 52]])
prof_name = 'Cross-section_15'
# 16
profile_A = np.array([[19, 43], [19, 52]])
prof_name = 'Cross-section_16'
# 17
profile_A = np.array([[20, 43], [20, 52]])
prof_name = 'Cross-section_17'
# 18
profile_A = np.array([[21, 43], [21, 52]])
prof_name = 'Cross-section_18'
# 19
profile_A = np.array([[22, 43], [22, 52]])
prof_name = 'Cross-section_19'
# 20
profile_A = np.array([[23, 43], [23, 52]])
prof_name = 'Cross-section_20'
# ##### 21 ##############
# profile_A = np.array([[2, 43], [23, 43]])
# prof_name = 'Cross-section_21'
# # # ##### 22 ##############
# profile_A = np.array([[2, 43.7], [23, 43.7]])
# prof_name = 'Cross-section_22'
# # # ##### 23 ##############
# profile_A = np.array([[2, 44.4], [23, 44.4]])
# prof_name = 'Cross-section_23'
# # # ##### 24 ##############
# profile_A = np.array([[2, 45.1], [23, 45.1]])
# prof_name = 'Cross-section_24'
# # ## ##### 25 ##############
# profile_A = np.array([[2, 45.8], [23, 45.8]])
# prof_name = 'Cross-section_25'
# # # ## ##### 26 ##############
# profile_A = np.array([[2, 46.5], [23, 46.5]])
# prof_name = 'Cross-section_26'
# # # ## ##### 27 ##############
# profile_A = np.array([[2, 47.2], [23, 47.2]])
# prof_name = 'Cross-section_27'
# # # ## ##### 28 ##############
# profile_A = np.array([[2, 47.9], [23, 47.9]])
# prof_name = 'Cross-section_28'
# # # ## ##### 29 ##############
# profile_A = np.array([[2, 48.6], [23, 48.6]])
# prof_name = 'Cross-section_29'
# # # ## ##### 30 ##############
# profile_A = np.array([[2, 49.3], [23, 49.3]])
# prof_name = 'Cross-section_30'
# # # ## ##### 31 ##############
# profile_A = np.array([[2, 50], [23, 50]])
# prof_name = 'Cross-section_31'
G2_, sta, xx, zz = plot_migration_sphr.create_2d_profile(mObs_ep, m_params, profile_A, sta, swath=37.5, plot=True)

################
# Smoothing    #
################
mObs = rf_mig.ccp_smooth(G2_, m_params)
# mObs[np.abs(mObs) < np.max(np.abs(mObs)) * 15 / 100] = 0
mObs = rf_mig.ccpFilter(mObs)
# ################
# # Plotting     #
# ################
# plot_migration_sphr.plot_migration_profile(Gp=mObs, xx=xx, zz=zz, migration_param_dict=m_params, sta=sta,
#                                            work_directory=work_dir, filename='iasp91', plot_title='iasp91')
######################################################################################
######################################################################################
# File for creating cross-sections with GMT
for i, x in enumerate(xx):
    for j, z in enumerate(zz):
        print(kilometers2degrees(x), z, mObs[i,j])
        with open(prof_name + '.txt', 'a') as of:
            of.write('{} {} {} \n'.
                     format(kilometers2degrees(x), z, mObs[i, j]))
######################################################################################
######################################################################################