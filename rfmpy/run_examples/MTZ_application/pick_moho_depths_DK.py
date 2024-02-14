"""
Pick Moho depths...

Location: Chavannes-pres-renens, CH
Date: Aug 2022 -- Modified Jan 2024
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
    codes_root_dir = '/github'
    desktop_dir = '/home/kmichailos/Desktop'
    hard_drive_dir = '/media/kmichailos/SEISMIC_DATA/'
else:
    data_root_dir = '/media/kmichall/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/github'
    desktop_dir = '/home/kmichall/Desktop'
    hard_drive_dir = '/media/kmichall/SEISMIC_DATA/'

# Define paths
work_dir = os.getcwd()
path = desktop_dir + "/RF_test/"

## Define MIGRATION parameters
# Ray-tracing parameters
inc = 2  # km
zmax = 750 # km
# Determine study area (x -> perpendicular to the profile)
minx = -4.0 # degrees 2.optional:2 or -4
maxx = 38.0 # degrees 2.optional:30 or 38
pasx = 0.38 # degrees
miny = 38.0 # degrees 2.optional:41 or 38
maxy = 54.0 # degrees 2.optional:51 or 54
pasy = 0.27 # degrees
minz = -5 # km
# maxz needs to be >= zmax
maxz = 750 # km
pasz = 2 # km
# Pass all the migration parameters in a dictionary to use them in functions called below
m_params = {'minx': minx, 'maxx': maxx,
            'pasx': pasx, 'pasy': pasy, 'miny': miny, 'maxy': maxy,
            'minz': minz, 'maxz': maxz, 'pasz': pasz, 'inc': inc, 'zmax': zmax}
# read stations
sta = rf_mig.read_stations_from_sac(path2rfs=path)


# COPY PASTE FROM HERE |
#                     /|\
with open('/home/kmichailos/Desktop/All_zmodel_m60_vel.npy', 'rb') as f:
    mObs_ep = np.load(f)

# with open('/home/kmichailos/Desktop/All_iasp91.npy', 'rb') as f:
#     mObs_ia = np.load(f)

# # 3D to 2D
# 0
profile_A = np.array([[8, 43], [8, 51]])
prof_name = 'Cross-section_0'

# swath=37.5
G2_, sta_, xx, zz = plot_migration_sphr.create_2d_profile_4_moho_picker(mObs_ep, m_params, profile_A, sta, swath=150, plot=True)
G2 = rf_mig.ccp_smooth(G2_, m_params)
# G2[np.abs(G2) < np.max(np.abs(G2)) * 15 / 100] = 0
G2 = rf_mig.ccpFilter(G2)
# Manually pick moho deps
# IMPORTANT NOTE: only works with cross-sections the have S-N and W-E directions!!!
plot_migration_sphr.moho_picker(Gp=G2, xx=xx, zz=zz, migration_param_dict=m_params,
                                sta=sta_, work_directory=work_dir, profile=profile_A, profile_name=prof_name,
                                path4file= work_dir )#)+ '/rfmpy/visualisation/gmt/maps/files/moho_picks/')
