import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Parameters import setup
from Tools import tools
from Tools import plotting_tools as pt
import time
from multiprocessing import Pool
from contextlib import closing

################
# DEFINE PATHS #
################

path = "/Users/mscarpon/Desktop/Projects/Ivrea/"
stations_list = "IA_StCoordinates.txt"  # List of stations
rf_list = "RRF.ZRT.10s.2Hz.list.txt"  # List of RFs
path_velocity_model = "/Users/mscarpon/Desktop/Projects/Ivrea/"  # 2D velocity model

"""Orientation of the profile with respect to the North direction (clockwise)
   ori_prof = 90 -> West-East oriented profile
   ori_prof = 0 -> North-South oriented profile """

ori_prof = 90

####################
# READ AND PROJECT #
####################
# Read stations and project
sta, dxSta, dySta = tools.Read_Seismic_Stations(file=path + "Seismic/IA_StCoordinates.txt", ori_prof=ori_prof)

lon_c = sta["LONSTA"].mean()  # Center of the profile
lat_c = sta["LATSTA"].mean()  # Center of the profile

# Read RFs
stream = tools.Read_Traces(RF_list=path + "Seismic/" + rf_list, sta=sta, ori_prof=ori_prof)


########################
# MIGRATION PARAMETERS #
########################
migration_params = setup.Migration_Params(dx=dxSta, dy=dySta)
print(migration_params)
minx, maxx, pasx = migration_params[:3]
miny, maxy, pasy = migration_params[3:6]
minz, maxz, pasz = migration_params[6:9]
inc, zmax = migration_params[9:11]
print(minx, maxx)
print(miny, maxy)
print(minz, maxz)
#######################
# PERFORM RAY TRACING #
#######################
# With a 2D velocity model
# stream = tools.tracing_2D(stream=stream, ori_prof=ori_prof, path_velocity_model=path_velocity_model,
#                           parameters=migration_params, lon_c=lon_c, lat_c=lat_c, dx=dxSta, dy=dySta,)
# Or with a 1D velocity model
stream = tools.tracing_1D(stream=stream, ori_prof=ori_prof,parameters=migration_params,
                          lon_c=lon_c, lat_c=lat_c, zMoho=50,)
###################################
# PERFORM TIME-TO-DEPTH MIGRATION #
###################################

#############
# Migrating #
#############

""" bazmean and dbaz are two key parameters to decide
    which RFs traces you want to include in the migration
    bazmean is the main backazimuthal direction (e.g. bazmean=90 to include traces from the East)
    and dbaz decides how wide is your selection interval, 
    so that you select everything that falls into bazmean +- dbaz
    bazmean == 180 and dbaz == 180 means you select everything (360 deg interval)"""

mObs = tools.ccpM(stream, migration_params, sta, phase="PS", stack=0, dbaz=180, bazmean=180)

#############
# Smoothing #
#############
mObs = tools.ccp_smooth(mObs, migration_params)
mObs[np.abs(mObs) < np.max(np.abs(mObs)) * 15 / 100] = 0
mObs = tools.ccpFilter(mObs)
############
# Plotting #
############
pt.Migration(Gp=mObs, parameters=migration_params, sta=sta, filename=False)
             # filename="/Users/mscarpon/Desktop/MigrationTest.png",