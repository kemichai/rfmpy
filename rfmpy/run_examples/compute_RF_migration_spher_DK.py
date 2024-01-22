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
import platform
import os
import time
import numpy as np
from obspy.taup import TauPyModel
from scipy.interpolate import RegularGridInterpolator
from scipy import interpolate
import obspy
import glob
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from obspy.geodetics.base import gps2dist_azimuth as gps2dist
from obspy.geodetics import degrees2kilometers, kilometers2degrees
import multiprocessing as mp
from obspy import Stream

def read_traces_sphr(path2rfs, sta):
    """
    Read all SAC files to get the seismic site details.

    :type path2rfs: str
    :param path2rfs: Path to the stored RF SAC files.
    :type sta: Pandas DataFrames.
    :param sta: Seismic site details.

    :return: A stream of the receiver function waveforms.
    """

    # Stations x and y values
    lonsta = sta["LONSTA"].values
    latsta = sta["LATSTA"].values
    altsta = sta["ZSTA"].values

    # Assign a unique number to each single separate station (index)
    enum = enumerate(sta["NAMESTA"].values)
    dictionary = dict((i, j) for j, i in enum)
    # Define model (iasp91)
    model = TauPyModel(model="iasp91")

    stream = obspy.Stream()
    all_rfs = glob.glob(path2rfs + '*.SAC')
    print("|-----------------------------------------------|")
    print("| Reading receiver functions...                 |")
    for i, rf in enumerate(all_rfs):
        trace_ = obspy.read(rf)
        trace = trace_[0]
        print(f"| Reading trace {i} of {len(all_rfs)}")
        # Define station's index number
        station_index = dictionary[trace.stats.station]
        trace.station_index = station_index
        # Station's name
        trace.kstnm = trace.stats.sac.kstnm
        # Delta value
        trace.delta = trace.stats.sac.delta
        trace.lon0 = lonsta[station_index]
        trace.lat0 = latsta[station_index]
        # Station's latitude
        trace.stla = trace.stats.sac.stla
        # Station's longitude
        trace.stlo = trace.stats.sac.stlo
        # Station's altitude in km
        trace.alt = sta["ZSTA"].values[station_index]
        # Back azimuth of station to eq
        trace.baz = trace.stats.sac.baz
        trace.stats.baz = trace.stats.sac.baz
        # Local back-azimuth (ori_prof is zero here)
        trace.lbaz = trace.baz
        trace.gcarc = trace.stats.sac.gcarc
        trace.depth = trace.stats.sac.evdp

        # Should not have distances smaller than 30 degrees...
        if trace.gcarc < 30:
            # Seconds per degrees
            try:
                trace.prai = model.get_travel_times(source_depth_in_km=trace.depth,
                                                distance_in_degree=trace.gcarc,
                                                phase_list=["Pn"], )[0].ray_param_sec_degree
            except Exception as e:
                print(e)
                trace.prai = -10

        # This is the distance range we use..
        elif trace.gcarc >= 30 and trace.gcarc <= 95.1:
            trace.prai = model.get_travel_times(source_depth_in_km=trace.depth,
                                                distance_in_degree=trace.gcarc,
                                                phase_list=["P"], )[0].ray_param_sec_degree
        # Should not have distances greater than 95 degrees...
        elif trace.gcarc > 95.1:
            try:
                trace.prai = model.get_travel_times(source_depth_in_km=trace.depth,
                                                    distance_in_degree=trace.gcarc,
                                                    phase_list=["PKIKP"], )[0].ray_param_sec_degree
            except Exception as e:
                print(e)
                trace.prai = -10

        trace.time = range(len(trace.data)) * trace.delta - 5
        trace.filename = rf
        trace.rms = np.sqrt(np.mean(np.square(trace.data)))

        stream.append(trace)
    print("|-----------------------------------------------|")

    return stream





def write_files_4_piercing_points_and_raypaths(st, sta, piercing_depth=35, plot=True):
    """..."""

    piercing_lon = []
    piercing_lat = []
    for i, tr in enumerate(st):
        # this try statement gets rid of RFs with no ray paths (ones with tr.Z = -1; tr.Xp = -1 etc)
        try:
            depth = tr.Z[0]
        except:
            continue
        for j, z in enumerate(tr.Z):
            if z > piercing_depth and z < piercing_depth + 15:
                piercing_lon.append(tr.Xp[j])
                piercing_lat.append(tr.Yp[j])
                with open('piercing_points.txt', 'a') as of:
                    of.write('{}, {}\n'.format(tr.Xp[j], tr.Yp[j]))

    wav_p_lon = []
    wav_p_lat = []
    wav_p_dep = []
    for i, tr in enumerate(st):
        try:
            depth = tr.Z[0]
        except:
            continue
        for j, z in enumerate(tr.Z):
                wav_p_lon.append(tr.Xp[j])
                wav_p_lat.append(tr.Yp[j])
                wav_p_dep.append(z)
                with open('ray_path.txt', 'a') as of:
                    of.write('{}, {}, {}\n'.
                            format(tr.Xp[j], tr.Yp[j], z))

    if plot:
        # Plot raypaths
        plt.scatter(wav_p_lon, wav_p_lat, alpha=0.5,
                    c=wav_p_dep, marker='.', edgecolor=None, s=5)
        plt.scatter(sta["LONSTA"], sta["LATSTA"],
                    c='r', marker='v', edgecolor='k', s=100)
        plt.show()

        # Plot piercing points
        plt.scatter(piercing_lon, piercing_lat, alpha=.3,
                    c='gray', marker='x', edgecolor='gray', s=50)
        plt.scatter(sta["LONSTA"], sta["LATSTA"],
                    c='r', marker='v', edgecolor='k', s=100)
        plt.show()

    return

t_beg = time.time()
# Set up paths
if platform.node().startswith('kmichailos-laptop'):
    data_root_dir = '/media/kmichailos/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/home/kmichailos/Desktop/codes/github'
    desktop_dir = '/home/kmichailos/Desktop'
    hard_drive_dir = '/media/kmichailos/SEISMIC_DATA/'
else:
    data_root_dir = '/media/konstantinos/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/home/Desktop/codes'
    desktop_dir = '/home/konstantinos/Desktop'
    hard_drive_dir = '/media/konstantinos/SEISMIC_DATA/'

# Define paths
work_dir = os.getcwd()
# Example RFs from a couple of teleseismic events
# path = work_dir + "/data/RF/RF/"
# Path to RFs in the hard drive
# path='/media/kmichailos/SEISMIC_DATA/RF_calculations/RF/'
# Path to RFs in the Desktop

path = desktop_dir + "/RF_test/"

#################
# Read stations #
#################
# Read station coordinates from the rfs (sac files) in a pandas dataframe
sta = rf_mig.read_stations_from_sac(path2rfs=path)

################
# Read RFs     #
################
stream = read_traces_sphr(path2rfs=path, sta=sta)


# Define MIGRATION parameters
# Ray-tracing parameters
inc = 5 # km
zmax = 800 # km
# Determine study area (x -> perpendicular to the profile)
minx = 0.0 # degrees
maxx = 40.0 # degrees
pasx = 0.23 # degrees
miny = 40.0 # degrees
maxy = 60.0 # degrees
pasy = 0.23 # degrees
minz = -5 # km
# maxz needs to be >= zmax
maxz = 800 # km
pasz = 10 # km
# Pass all the migration parameters in a dictionary to use them in functions called below
m_params = {'minx': minx, 'maxx': maxx,
            'pasx': pasx, 'pasy': pasy, 'miny': miny, 'maxy': maxy,
            'minz': minz, 'maxz': maxz, 'pasz': pasz, 'inc': inc, 'zmax': zmax}
################
# Ray tracing  #
################
# Pick one of the two velocity models
# 'EPcrust' or 'iasp91' or 'zmodel_m60'
stream_ray_trace = rf_mig.tracing_3D_sphr(stream=stream, migration_param_dict=m_params,
                                          velocity_model='zmodel_m60')

# Write piercing points in a file
write_files_4_piercing_points_and_raypaths(stream_ray_trace, sta, piercing_depth=40, plot=False)
################
# Migration    #
################
mObs = rf_mig.ccpm_3d(stream_ray_trace, m_params, output_file="/home/kmichailos/Desktop/All_zmodel_m60_vel", phase="PS")


total_time = time.time() - t_beg
print('Ray tracing took ' + str(round(total_time)/60) + ' minutes in total.')
