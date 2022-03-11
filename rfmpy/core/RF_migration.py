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


def Migration_Params(dx, dy):

    # SETTING GEOMETRY PARAMETERS FOR MIGRATION

    minx = -120 + dx
    maxx = 120 + dx
    pasx = 0.50

    miny = -100 + dy
    maxy = 100 + dy
    pasy = maxy - miny

    minz = -2
    maxz = 100
    pasz = 0.50

    inc = 0.250
    zmax = 100  # !! always <= than maxz

    params = (minx, maxx, pasx, miny, maxy, pasy, minz, maxz, pasz, inc, zmax)

    return params

# Define paths
work_dir = os.getcwd()
path = work_dir + "/data/RF/"
# To be removed
rf_list = 'RRF_list.txt'

ori_prof = 90
sta, dxSta, dySta = migration_utils.read_stations(path2rfs=path, ori_prof=ori_prof)

lon_c = sta["LONSTA"].mean()  # Center of the profile
lat_c = sta["LATSTA"].mean()  # Center of the profile

stream = migration_utils.Read_Traces(path2rfs=path, sta=sta, ori_prof=ori_prof)

# MIGRATION PARAMETERS #
migration_params = Migration_Params(dx=dxSta, dy=dySta)
print(migration_params)
minx, maxx, pasx = migration_params[:3]
miny, maxy, pasy = migration_params[3:6]
minz, maxz, pasz = migration_params[6:9]
inc, zmax = migration_params[9:11]
print(minx, maxx)
print(miny, maxy)
print(minz, maxz)
# PERFORM RAY TRACING #
stream = migration_utils.tracing_1D(stream=stream, ori_prof=ori_prof,parameters=migration_params,
                          lon_c=lon_c, lat_c=lat_c, zMoho=50,)
# Migrating #
mObs = migration_utils.ccpM(stream, migration_params, sta, phase="PS", stack=0, dbaz=180, bazmean=180)
# Smoothing #
#############
mObs = migration_utils.ccp_smooth(mObs, migration_params)
mObs[np.abs(mObs) < np.max(np.abs(mObs)) * 15 / 100] = 0
mObs = migration_utils.ccpFilter(mObs)
# Plotting #
############
migration_utils.Migration(Gp=mObs, parameters=migration_params, sta=sta, work_directory=work_dir,
                          filename=False)
