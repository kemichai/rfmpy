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
# path = work_dir + "/data/RF/RF/"
# path = desktop_dir + "/RF_test/RF/"
# path='/media/kmichailos/SEISMIC_DATA/RF_calculations/RF/'
path = desktop_dir + "/all_rfs/RF/"
# path='/media/kmichailos/SEISMIC_DATA/RF_calculations/RF_low_quality/'
#################
# Read stations #
#################
# Read station coordinates from the rfs (sac files) in a pandas dataframe
sta = rf_mig.read_stations_from_sac(path2rfs=path)

# plt.scatter(sta["LONSTA"], sta["LATSTA"],
#             c='r', marker='v', edgecolor='k', s=100)
# plt.show()



################
# Read RFs     #
################
# List of stations to exclude during migration due to noisy data
# this list of stations was defined by manually inspecting RF traces vs their back-azimuths
# TODO: finish this list...
# station_to_exclude = ['CH.WOLEN', 'CR.RABC', 'FR.NEEW', 'GU.CARD',
#                       'HU.BSZH', 'IV.SARZ', 'IV.ZCCA', 'Z3.A153A', 'YP.CT37'
#                       'Z3.A251A',
#                       ]
# 'CR.ZAG', 'CZ.PBCC', 'FR.BETS', 'FR.ASEAF', 'OX.BALD',
# 'SI.BOSI', 'SK.KOLS', 'XT.AAE03', 'Z3.A088A',
# Z3.400> the ones with 400> OBS ones
# Maybe remove this ...'CH.WOLEN''FR.NEEW'

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

# 'EPcrust' or 'iasp91'
stream_ray_trace = rf_mig.tracing_3D_sphr(stream=stream, migration_param_dict=m_params,
                                          velocity_model='EPcrust')

#
# piercing_lon = []
# piercing_lat = []
# for i, tr in enumerate(stream_ray_trace):
#     tr.stats.station
#     for j, z in enumerate(tr.Z):
#         if z > 34 and z < 35:
#             # print(tr.Xp[j], tr.Yp[j])
#             piercing_lon.append(tr.Xp[j])
#             piercing_lat.append(tr.Yp[j])
# plt.scatter(piercing_lon, piercing_lat, alpha=.3,
#             c='gray', marker='x', edgecolor='gray', s=50)
# plt.scatter(sta["LONSTA"], sta["LATSTA"],
#             c='r', marker='v', edgecolor='k', s=100)
# plt.show()
# #
# wav_p_lon = []
# wav_p_lat = []
# wav_p_dep = []
# for i, tr in enumerate(stream_ray_trace):
#     tr.stats.station
#     for j, z in enumerate(tr.Z):
#             # print(tr.Xp[j], tr.Yp[j])
#             wav_p_lon.append(tr.Xp[j])
#             wav_p_lat.append(tr.Yp[j])
#             wav_p_dep.append(z)
#
#
#
# plt.scatter(wav_p_lon, wav_p_lat, alpha=0.5,
#             c=wav_p_dep, marker='.', edgecolor=None, s=1)
# plt.scatter(sta["LONSTA"], sta["LATSTA"],
#             c='r', marker='v', edgecolor='k', s=100)
# plt.show()

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
