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
client=Client()

starttime = UTCDateTime("2015-01-01")
endtime = UTCDateTime("2019-01-02")
inventory = client.get_stations(network="BW", station="*",
                                starttime=starttime,
                                endtime=endtime)
inventory.plot('local')

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