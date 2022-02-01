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

all_event_dir = glob.glob(path_ev + '*')
event_dir = all_event_dir[3]

print('Calculating RF for event in: ', event_dir)
station_list = rf_util.get_unique_stations(event_dir)
east_comp_traces = obspy.Stream()
north_comp_traces = obspy.Stream()
vert_comp_traces = obspy.Stream()
for station in station_list:
    single_cha_trace = obspy.read(event_dir + '/*' + station + '*')
    if single_cha_trace[0].stats.channel[-1] == 'Z':
        vert_comp_traces.append(single_cha_trace[0])
    if single_cha_trace[0].stats.channel[-1] == 'E':
        east_comp_traces.append(single_cha_trace[0])
    if single_cha_trace[0].stats.channel[-1] == 'N':
        north_comp_traces.append(single_cha_trace[0])
# todo: Rotate to real N and E
from obspy.signal.rotate import rotate2zne
from obspy import Stream
from obspy import read_inventory

for i, z_trace in enumerate(vert_comp_traces):
    print(z_trace, north_comp_traces[i])
    tr1 = north_comp_traces[i]
    tr2 = east_comp_traces[i]
    tr3 = z_trace
    orig_stream = Stream()
    orig_stream.append(tr1)
    orig_stream.append(tr2)
    orig_stream.append(tr3)
    # orig_stream.plot()
    ####################
    # READING stationxml
    inv = read_inventory('/home/kmichall/Desktop/Codes/github/rfmpy/rfmpy/metadata/RDAlparray.xml')
    print('>>> Read inventory...')
    e_trace_name = tr1.stats.network + '.' + tr1.stats.station + '.' + tr1.stats.channel
    e_trace_time = tr1.stats.starttime
    for net in inv:
        for sta in net:
            for cha in sta:
                cha_name = net.code + '.' + sta.code + '.' + cha.code
                if e_trace_name == cha_name and e_trace_time > cha.start_date and e_trace_time < cha.end_date:
                    print(cha_name, cha.azimuth, cha.dip)

    ####################

    tr_e = tr1.copy()
    tr_n = tr2.copy()
    tr_z = tr3.copy()
    tr_z.data, tr_n.data, tr_e.data = rotate2zne(tr3.data, 0, -90,
                                                 tr2.data, 0, 0,
                                                 tr1.data, 90, 0, inverse=False)
    rot_stream = Stream()
    rot_stream.append(tr_e)
    rot_stream.append(tr_n)
    rot_stream.append(tr_z)
    all = Stream()
    all = orig_stream + rot_stream
    all.plot()

