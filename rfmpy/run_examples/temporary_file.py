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
# List of unique seismic sites
sta1 = rf_util.get_station_info(path_wavs_list_part1)
sta2 = rf_util.get_station_info(path_wavs_list_part2)
sta3 = rf_util.get_station_info(path_wavs_list_part3)
sta4 = rf_util.get_station_info(path_wavs_list_part4)
sta5 = rf_util.get_station_info(path_wavs_list_part5)


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

from obspy import UTCDateTime
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

from obspy.clients.fdsn import RoutingClient
from obspy.clients.fdsn import Client
client=RoutingClient('eida-routing')
client=Client('ETH')

starttime = UTCDateTime("2015-01-01")
endtime = UTCDateTime("2019-01-02")
inventory = client.get_stations(network="Z3*", station="*",
                                level='channel',
                                starttime=starttime,
                                endtime=endtime)
inventory.plot('local')
for cha in inventory[0].stations[0].channels:
    print(cha.azimuth)

inventory.write('Z3_eth.xml',format='STATIONXML')

# GFZ / BGR / LMU nodes archives data collected by German institutions
# ODC node archives data collected by institutions from Austria, Hungary, Czech Republic.


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

stream = read(path_wavs + 'P_2015.002.08.21.55/' + '*ZUR*')

from obspy.signal.rotate import rotate2zne

a, b, c = rotate2zne(stream[2], 30, -90, stream[0], 90, 3, stream[1], 92, 3)
