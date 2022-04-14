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

plt.scatter(sta["LONSTA"], sta["LATSTA"], c='r', marker='v')
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
# plot_migration_sphr.plot_migration_profile(Gp=mObs, migration_param_dict=m_params, sta=sta,
#                                       work_directory=work_dir, filename=False)
#



#############################################
#### CREATE 2D grid profile from the 3d grid???
# OR DO ALL THAT IN SPHERICAL SOMEHOW/???

# Define profile for migration
profile_az = 0
# Define point (lon, lat) to project the station's coordinates with respect
profile_lon = sta["LONSTA"].mean()
profile_lat = sta["LATSTA"].mean()



#### CREATE 2D grid profile from the 3d grid
    # Read migration parameters
    minx = migration_param_dict['minx']
    maxx = migration_param_dict['maxx']
    pasx = migration_param_dict['pasx']
    miny = migration_param_dict['miny']
    maxy = migration_param_dict['maxy']
    pasy = migration_param_dict['pasy']
    minz = migration_param_dict['minz']
    maxz = migration_param_dict['maxz']
    pasz = migration_param_dict['pasz']
    inc = migration_param_dict['inc']
    zmax = migration_param_dict['zmax']

    # Grid preparation
    xx = np.arange(minx, maxx + pasx, pasx)
    yy = np.arange(miny, maxy + pasy, pasy)
    zz = np.arange(minz, maxz + pasz, pasz)
    # G

    pts = np.array([7.4, 46, 11])
    Amp = G_interpolated(pts)


    ref_pnt = np.array([[xx[0], yy[0]], [xx[-1], yy[-1]]])
    profile_width = 100
    num_events = len(xx)
    lon = xx
    lat = yy

    prof_dist, prof_dep = [], []
    cos_lat = np.cos(ref_pnt[0][1] * np.pi / 180)
    vec_ab = ref_pnt[1] - ref_pnt[0]
    vec_ab[0] *= cos_lat
    abs_ab = np.linalg.norm(vec_ab)
    for i in range(num_events):
        loc_c = np.array([lon[i], lat[i]])
        vec_ac = loc_c - ref_pnt[0]
        vec_ac[0] *= cos_lat
        abs_ac = np.linalg.norm(vec_ac)
        cos = vec_ac.dot(vec_ab) / abs_ab / abs_ac
        if abs_ac * (1 - cos ** 2) ** 0.5 > profile_width / 111.: continue
        if cos < 0 or abs_ac * cos > abs_ab: continue
        prof_dist.append(abs_ac * cos * 111)



    XX, ZZ = np.meshgrid(prof_dist[1:], prof_dep[1:])
