"""
Code to compute receiver functions.

NOTES: Reads data as they were stored using the prepare_wavs_4_rf.py
       script. ../EVENTS/YEAR.JULDAY.HOUR.MINUTE.SECOND/YEAR.JULDAY.HOUR.MINUTE.SECOND.STATION.COMPONENT.SAC

=============================================
Requirements:
    * ObsPy
    * seismic waveform data stored in a
      specific way (see code below)
    *
=============================================

Original code written by Matteo Scarponi on 30.11.2021
Location: Chavannes-pres-renens, CH
Date: Jan 2022
Author: Konstantinos Michailos
"""

import core.main.RF_Main as RF
import platform
from obspy import read_inventory, read_events, UTCDateTime as UTC

# Set up parameters and paths
if platform.node().startswith('kmichailos-laptop'):
    data_root_dir = '/media/kmichailos/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/home/kmichailos/Desktop/codes/bitbucket'
    desktop_dir = '/home/kmichailos/Desktop'
else:
    data_root_dir = '/media/kmichall/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/home/kmichall/Desktop/Codes/bitbucket'
    desktop_dir = '/home/kmichall/Desktop'

# Path in which waveforms are stored
path_wavs = '/media/kmichall/SEISMIC_DATA/RF_data/DATA_RFAA_part_1/SWISS/data/'
# path_wavs = desktop_dir + '/RF_test/EVENTS/'
path_wavs = desktop_dir + '/RF_test/test_data/'

# Path to store RFs
path_out_RF = '/media/kmichall/SEISMIC_DATA/RF_calculations/'
path_out_RF = desktop_dir + '/RF_test/RF_1/'

a = RF.calculate_rf(path_ev=path_wavs, path_out=path_out_RF,
                    iterations=100, c1=10, c2=10, c3=1, c4=1,
                    max_frequency=2, save=True)

