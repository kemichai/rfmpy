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
with open('/home/kmichailos/Desktop/All_EPcrust.npy', 'rb') as f:
    mObs_ep = np.load(f)
with open('/home/kmichailos/Desktop/All_iasp91.npy', 'rb') as f:
    mObs_ia = np.load(f)

# 3D to 2D
profile_A = np.array([[8, 45.5], [15, 50]])
G2_, sta, xx, zz = plot_migration_sphr.create_2d_profile(mObs_ep, m_params, profile_A, sta, swath=25, plot=True)

################
# Smoothing    #
################
mObs = rf_mig.ccp_smooth(G2, m_params)
# mObs[np.abs(mObs) < np.max(np.abs(mObs)) * 15 / 100] = 0
mObs = rf_mig.ccpFilter(mObs)
# ################
# # Plotting     #
# ################
plot_migration_sphr.plot_migration_profile(Gp=G2, xx=xx, zz=zz, migration_param_dict=m_params, sta=sta,
                                           work_directory=work_dir, filename='EPcrust', plot_title='EPcrust')
######################################################################################
######################################################################################
# File for creating cross-sections with GMT
for i, x in enumerate(xx):
    for j, z in enumerate(zz):
        print(kilometers2degrees(x), z, mObs[i,j])
        with open('xyz_smoothed_test.txt', 'a') as of:
            of.write('{} {} {} \n'.
                     format(kilometers2degrees(x), z, mObs[i, j]))
######################################################################################
######################################################################################