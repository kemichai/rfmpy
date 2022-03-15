"""
Function for calculating 3D migration of RFs in spherical coordinates for the AlpArray

Note: Based on codes originally written by Matteo Scarponi.

Location: Chavannes-pres-renens, CH
Date: Mar 2022
Author: Konstantinos Michailos
"""
import obspy
import glob
import numpy as np
from rfmpy.utils import migration_utils
import platform
import os

# Set up paths
if platform.node().startswith('kmichailos-laptop'):
    data_root_dir = '/media/kmichailos/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/home/kmichailos/Desktop/codes/github'
    desktop_dir = '/home/kmichailos/Desktop'
    hard_drive_dir = '/media/kmichailos/SEISMIC_DATA/'
else:
    data_root_dir = '/media/kmichall/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/home/kmichall/Desktop/Codes/github'
    desktop_dir = '/home/kmichall/Desktop'
    hard_drive_dir = '/media/kmichall/SEISMIC_DATA/'

# Define paths
work_dir = os.getcwd()
path = work_dir + "/data/RF/"

# TODO: look at what this parameters stands for...
ori_prof = 90
prof_azimuth = 90

sta, dxSta, dySta = migration_utils.read_stations(path2rfs=path, ori_prof=ori_prof)
lon_c = sta["LONSTA"].mean()  # Center of the profile
lat_c = sta["LATSTA"].mean()  # Center of the profile

# Read RF files
stream = migration_utils.Read_Traces(path2rfs=path, sta=sta, ori_prof=ori_prof)

# Define migration parameters
# Ray-tracing parameters
inc = 0.250
zmax = 100
# Determine study area (x -> perpendicular to the profile)
minx = -120 + dxSta
maxx = 120 + dxSta
pasx = 0.50
miny = -100 + dySta
maxy = 100 + dySta
pasy = maxy - miny
minz = -2
# maxz needs to be >= zmax
maxz = 100
pasz = 0.50
# Pass all the migration parameters in a dictionary to use them in functions
m_params = {'minx': minx, 'maxx': maxx, 'pasx': pasx, 'pasy': maxy-miny, 'miny': miny, 'maxy': maxy,
            'minz': minz, 'maxz': maxz, 'pasz': pasz, 'inc': inc, 'zmax': zmax}

# Ray tracing
stream_ray_trace = migration_utils.tracing_1D(tr=stream, ori_prof=ori_prof,
                                              migration_param_dict=m_params,
                                              lon_c=lon_c, lat_c=lat_c, zMoho=50,)

# stream = tools.tracing_2D(stream=stream, ori_prof=ori_prof, path_velocity_model=work_dir,
#                           parameters=migration_params, lon_c=lon_c, lat_c=lat_c, dx=dxSta, dy=dySta,)
# Or with a 1D velocity model

# Migration
mObs = migration_utils.ccpM(stream_ray_trace, m_params, sta, phase="PS", stack=0, dbaz=180, bazmean=180)
# Smoothing
mObs = migration_utils.ccp_smooth(mObs, m_params)
mObs[np.abs(mObs) < np.max(np.abs(mObs)) * 15 / 100] = 0
mObs = migration_utils.ccpFilter(mObs)
# Plotting
migration_utils.Migration(Gp=mObs, migration_param_dict=m_params, sta=sta,
                          work_directory=work_dir,
                          filename=False)
