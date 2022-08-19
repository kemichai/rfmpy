"""
Pick Moho depths...

Location: Chavannes-pres-renens, CH
Date: Aug 2022
Author: Konstantinos Michailos
"""
import matplotlib
matplotlib.use('TkAgg')
import rfmpy.core.migration_sphr as rf_mig
import rfmpy.utils.migration_plots_spher as plot_migration_sphr
import numpy as np
import platform
import os
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
# TODO: aaa
with open('/home/kmichailos/Desktop/All_EPcrust_new_mantle_vel.npy', 'rb') as f:
    mObs_ep = np.load(f)

with open('/home/kmichailos/Desktop/All_iasp91.npy', 'rb') as f:
    mObs_ia = np.load(f)



# 3D to 2D
# 1
profile_A = np.array([[5, 43], [5, 50]])
prof_name = 'Cross-section 1'
# 2
profile_A = np.array([[7.5, 43], [7.5, 50]])
# 3
profile_A = np.array([[10, 43], [10, 50]])
# # 4
# profile_A = np.array([[12.5, 43], [12.5, 50]])
# # 5
# profile_A = np.array([[15, 43], [15, 50]])
# # 6
# profile_A = np.array([[17.5, 43], [17.5, 50]])
# # A
# profile_A = np.array([[2, 42.5], [20, 42.5]])
# # B
# profile_A = np.array([[2, 45], [10, 45]])
# # C
# profile_A = np.array([[2, 47.5], [20, 47.5]])
# # D
# profile_A = np.array([[2, 50], [20, 50]])


G2_, sta_, xx, zz = plot_migration_sphr.create_2d_profile_4_moho_picker(mObs_ep, m_params, profile_A, sta, swath=37.5, plot=True)

G2 = rf_mig.ccp_smooth(G2_, m_params)
# G2[np.abs(G2) < np.max(np.abs(G2)) * 15 / 100] = 0
G2 = rf_mig.ccpFilter(G2)
# Manually pick moho deps
# IMPORTANT NOTE: only works with cross-sections the have S-N and W-E directions!!!
plot_migration_sphr.moho_picker(Gp=G2, xx=xx, zz=zz, migration_param_dict=m_params,
                                sta=sta_, work_directory=work_dir, profile=profile_A,
                                profile_name=prof_name)
