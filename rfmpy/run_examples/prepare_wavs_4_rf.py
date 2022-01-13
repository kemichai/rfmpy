"""
Code to create the earthquake database we will later use for the
receiver function calculations. Many of the functions called in this main code,
are imported from RF_Util.py which provides a sequence of function definitions.

NOTES:
    1) This code reads or downloads a catalog of seismic events according to the given parameters.
    2) Goes through the seismic raw database, where data is supposed to be stored in the form of daily files.
    3) Reads the daily file, cuts around the expected earthquake arrival time, plus minus a certain.
       window as indicated in the parameters.
    4) Removes instrument response, filters, resamples and save the Z,E,N components to the EVENTS database folder.

=============================================
Requirements:
    * obspy
    * seismic waveform data stored in a
      specific way (see code below)
    *
=============================================

Original codes by Matteo Scarponi on 30.11.2021
Location: Chavannes-pres-renens, CH
Date: Dec 2021
Author: Konstantinos Michailos
"""

from obspy.clients.fdsn import Client
from rfmpy.core.cut_waveforms import prep_wavs4rf as cut_wavs
import platform
import os
from obspy import read_inventory, read_events, UTCDateTime as UTC
import glob

# Set up parameters and paths
if platform.node().startswith('kmichailos-laptop'):
    data_root_dir = '/media/kmichailos/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/home/kmichailos/Desktop/codes/bitbucket'
    desktop_dir = '/home/kmichailos/Desktop'
else:
    data_root_dir = '/media/kmichall/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/home/kmichall/Desktop/Codes/bitbucket'
    desktop_dir = '/home/kmichall/Desktop'

# Path to raw waveform data
wav_path = '/media/kmichall/SEISMIC_DATA/Data_archive/'
# Path to store traces
path_out = desktop_dir + '/RF_test/EVENTS/'

# Path to store quakeml catalog of teleseismic events
path_cat = desktop_dir + '/RF_test/'
# Path to .xml metadata files for each station
path_metadata = desktop_dir + '/RF_test/stations_metadata/'
path_out_RF = desktop_dir + '/RF_test/RF/'

# Download stationxml inventory file if it doesn't exist
inv_file = path_metadata + 'stationxml.xml'
if not os.path.exists(inv_file):
    client = Client('IRIS')
    inv = client.get_stations(network='YL', channel='BH?', level='channel', minlatitude=26.0, maxlatitude=30.0)
    inv.write(inv_file, 'STATIONXML')
inv = read_inventory(inv_file)
# inventory.plot(label=False)
# fig = inv.plot('local', color_per_network=True)

# Download or read tele-seismic earthquake catalog
coords = inv.get_coordinates('YL.BIRA..BHZ')
lonlat = (coords['longitude'], coords['latitude'])
# Catalog file path
cat_file = path_cat + 'temp_cat.xml'
if not os.path.exists(cat_file):
    client = Client()
    kwargs = {'starttime': UTC('2002-01-01'), 'endtime': UTC('2002-01-15'), 'latitude': lonlat[1],
              'longitude': lonlat[0], 'minradius': 30, 'maxradius': 90, 'minmagnitude': 6.0, 'maxmagnitude': 8.0}
    cat = client.get_events(**kwargs)
    cat.write(cat_file, 'QUAKEML')
cat = read_events(cat_file)
# fig = cat.plot('local')

run_code = cut_wavs(catalog=cat, inventory=inv, wav_directory=wav_path, output_dir=path_out)
