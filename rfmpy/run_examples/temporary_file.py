"""
Create list of seismic sites used.
Plot seismic network.

Location: Chavannes-pres-renens, CH
Date: Jan 2022
Author: Konstantinos Michailos
"""
import rfmpy.utils.RF_Util as rf_util
import platform

# Set up parameters and paths
if platform.node().startswith('kmichailos-laptop'):
    data_root_dir = '/media/kmichailos/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/home/kmichailos/Desktop/codes/bitbucket'
    desktop_dir = '/home/kmichailos/Desktop'
    hard_drive_dir = '/media/kmichailos/SEISMIC_DATA/'
else:
    data_root_dir = '/media/kmichall/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/home/kmichall/Desktop/Codes/bitbucket'
    desktop_dir = '/home/kmichall/Desktop'
    hard_drive_dir = '/media/kmichall/SEISMIC_DATA/'

# DATA_RFAA_part_1
path_wavs_list_part1 = [hard_drive_dir + 'RF_data/DATA_RFAA_part_1/SWISS/data/',
                        hard_drive_dir + 'RF_data/DATA_RFAA_part_1/EASI/easi_data/',
                        hard_drive_dir + 'RF_data/DATA_RFAA_part_1/FRANCE/data_sort/',
                        hard_drive_dir + 'RF_data/DATA_RFAA_part_1/FRANCE/data/',
                        hard_drive_dir + 'RF_data/DATA_RFAA_part_1/North_Italy/events_fri_ven/']

path_wavs_list_part2 = [hard_drive_dir + 'RF_data/DATA_RFAA_part_2/Austria/data_AAA_corrected/',
                        hard_drive_dir + 'RF_data/DATA_RFAA_part_2/CIFAlps/data_YP2012/',
                        hard_drive_dir + 'RF_data/DATA_RFAA_part_2/data_DINAR/',
                        hard_drive_dir + 'RF_data/DATA_RFAA_part_2/HU_SK/data/']
# DATA_RFAA_part_3
path_wavs_list_part3 = [hard_drive_dir + 'RF_data/DATA_RFAA_part_3/AARF/DATA_MOBST/data/',
                        hard_drive_dir + 'RF_data/DATA_RFAA_part_3/AARF/DATA_PERMST/data/',
                        hard_drive_dir + 'RF_data/DATA_RFAA_part_3/GERMANY/DE_AA_RF/DATA/data/']
path_wavs_list_part4 = [hard_drive_dir + 'RF_data/CIFALPS/data_YP2012/']
# INGV
path_wavs_list_part5 = [hard_drive_dir + 'RF_data/INGV-Permanent-data/',
                        hard_drive_dir + 'RF_data/INGV-Temporary-data/data/']

# New path for FR stations 27th of Jan
path_wavs_list_part_ = [hard_drive_dir + 'RF_data/FR_new/data/']
path_wavs_list_part__ = [hard_drive_dir + 'RF_data/FR_new/data_sort/']

# List of unique seismic sites
sta1 = rf_util.get_station_info(path_wavs_list_part1)
sta2 = rf_util.get_station_info(path_wavs_list_part2)
sta3 = rf_util.get_station_info(path_wavs_list_part3)
sta4 = rf_util.get_station_info(path_wavs_list_part4)
sta5 = rf_util.get_station_info(path_wavs_list_part5)
sta_ = rf_util.get_station_info(path_wavs_list_part_)
sta__ = rf_util.get_station_info(path_wavs_list_part__)


sta = sta1 + sta2 + sta3 + sta4 + sta5

unique_all_sta = []
for s in sta:
    if s not in unique_all_sta:
        unique_all_sta.append(s)
# using this for making the gmt plot

# For plotting see
# wiggle_bins functions in miscellaneous





"""
Download metadata information to doublecheck the horizontal components are oriented to
North and East. Will need to check this before I apply the rotations.

Location: Chavannes-pres-renens, CH
Date: Jan 2022
Author: Konstantinos Michailos
"""

from obspy.clients.fdsn import Client
import platform
from obspy import read_inventory, read_events, read

# Set up parameters and paths
if platform.node().startswith('kmichailos-laptop'):
    data_root_dir = '/media/kmichailos/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/home/kmichailos/Desktop/codes/bitbucket'
    desktop_dir = '/home/kmichailos/Desktop'
    hard_drive_dir = '/media/kmichailos/SEISMIC_DATA/'
else:
    data_root_dir = '/media/kmichall/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/home/kmichall/Desktop/Codes/bitbucket'
    desktop_dir = '/home/kmichall/Desktop'
    hard_drive_dir = '/media/kmichall/SEISMIC_DATA/'



azimuths = []
dips = []
with open('Z3.txt', 'r') as f:
    for line in f:
        if line.startswith('#'):
            print(line)
            continue
        else:
            ln = line.split('|')
            az = ln[8]
            dp = ln[9]
            dips.append(dp)
            azimuths.append(az)
            print(ln[1],ln[3], az, dp)


path_wavs = '/media/kmichall/SEISMIC_DATA/RF_data/DATA_RFAA_part_1/SWISS/data/'

from obspy.signal.rotate import rotate2zne

# for each station we need the azimuth and the dip...

stream = read(path_wavs + 'P_2015.002.08.21.55/' + '*ZUR*')
tr1 = stream[0]
tr2 = stream[1]
tr3 = stream[2]

tr_e = tr1.copy()
tr_n = tr2.copy()
tr_z = tr3.copy()
tr_e.data, tr_n.data, tr_z.data = rotate2zne(tr1.data, 30, -90, tr2.data, 90, 3, tr3.data, 92, 3, inverse=False)


#####################################
"""
Reorient part WIP

"""
import rfmpy.core.RF_Main as RF
import platform
from obspy import read_inventory, read_events, UTCDateTime as UTC
import rfmpy.utils.RF_Util as rf_util
from rfmpy.utils.signal_processing import rotate_trace, remove_response, ConvGauss
from rfmpy.utils.qc import *
from obspy import read_inventory, read_events, UTCDateTime as UTC
import itertools
from pathlib import Path
import glob
import obspy
import numpy as np



path_wavs = '/media/kmichall/SEISMIC_DATA/RF_data/DATA_RFAA_part_1/FRANCE/data_sort/'
path_ev=path_wavs
all_event_dir = glob.glob(path_ev + '*')
all_event_dir.sort()
event_dir = all_event_dir[119]


# todo: Rotate to real N and E
from obspy.signal.rotate import rotate2zne
from obspy import Stream
from obspy import read_inventory

inv = read_inventory('/home/kmichall/Desktop/Codes/github/rfmpy/rfmpy/metadata/*.xml')

for event_dir in all_event_dir:
    vert_comp_traces, north_comp_traces, east_comp_traces = rf_util.get_unique_stations(event_dir)
    print(event_dir)


# def correct_orientations(east, north, vertical):
    v_corr = []
    e_corr = []
    n_corr = []
    for i, trace in enumerate(vert_comp_traces):
        # print(z_trace, north_comp_traces[i])
        trace_n = north_comp_traces[i]
        trace_e = east_comp_traces[i]
        trace_z = vert_comp_traces[i]

        orig_stream = Stream()
        orig_stream.append(trace_e)
        orig_stream.append(trace_n)
        orig_stream.append(trace_z)
        # orig_stream.plot()
        ####################
        def get_trace_name(tr):
            return tr.stats.network + '.' + tr.stats.station + '.' + tr.stats.location + '.' + tr.stats.channel
        e_trace_name = get_trace_name(trace_e)
        n_trace_name = get_trace_name(trace_n)
        z_trace_name = get_trace_name(trace_z)
        # Time of trace to choose the right epoch
        trace_time = trace.stats.starttime
        for net in inv:
            for sta in net:
                for cha in sta:
                    cha_name = net.code + '.' + sta.code + '.' + cha.location_code + '.' + cha.code
                    if e_trace_name == cha_name and trace_time > cha.start_date and trace_time < cha.end_date:
                        e_trace_az = cha.azimuth
                        e_trace_dip = cha.dip
                    if n_trace_name == cha_name and trace_time > cha.start_date and trace_time < cha.end_date:
                        n_trace_az = cha.azimuth
                        n_trace_dip = cha.dip
                    if z_trace_name == cha_name and trace_time > cha.start_date and trace_time < cha.end_date:
                        z_trace_az = cha.azimuth
                        z_trace_dip = cha.dip


        tr_e = trace_e.copy()
        tr_n = trace_n.copy()
        tr_z = trace_z.copy()
        tr_z.data, tr_n.data, tr_e.data = rotate2zne(trace_z.data, z_trace_az, z_trace_dip,
                                                     trace_n.data, n_trace_az, n_trace_dip,
                                                     trace_e.data, e_trace_az, e_trace_dip,
                                                     inverse=False)
        rot_stream = Stream()
        rot_stream.append(tr_e)
        rot_stream.append(tr_n)
        rot_stream.append(tr_z)
        all = Stream()
        all = orig_stream + rot_stream
        if n_trace_az != 0.0 or e_trace_az != 90.0:
            print(n_trace_az, e_trace_az)
            all.plot()
        # Change channel names for borehole sites
        if tr_n.stats.channel[-1] != 'N':
            print(tr_n.stats.channel)
            tr_n.stats.channel = tr_n.stats.channel[0:2] + 'N'
            tr_n.stats.sac.kcmpnm = tr_n.stats.sac.kcmpnm[0:2] + 'N'
        if tr_e.stats.channel[-1] != 'E':
            print(tr_e.stats.channel)
            tr_e.stats.channel = tr_n.stats.channel[0:2] + 'E'
            tr_e.stats.sac.kcmpnm = tr_n.stats.sac.kcmpnm[0:2] + 'E'

        v_corr.append(tr_z)
        e_corr.append(tr_e)
        n_corr.append(tr_n)


    return e_corr, n_corr, v_corr

