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

import rfmpy.core.RF_Main as RF
import platform
from obspy import read_inventory, read_events, UTCDateTime as UTC

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

# Path in which waveforms are stored
path_wavs = '/media/kmichall/SEISMIC_DATA/RF_data/DATA_RFAA_part_1/SWISS/data/'
# path_wavs = desktop_dir + '/RF_test/EVENTS/'
# path_wavs = desktop_dir + '/RF_test/test_data/'
# DATA_RFAA_part_1
path_wavs_list_part_1 = [hard_drive_dir + 'RF_data/DATA_RFAA_part_1/SWISS/data/',
                         hard_drive_dir + 'RF_data/DATA_RFAA_part_1/EASI/easi_data/',
                         hard_drive_dir + 'RF_data/DATA_RFAA_part_1/FRANCE/data_sort/',
                         hard_drive_dir + 'RF_data/DATA_RFAA_part_1/North_Italy/events_fri_ven/']

path_wavs_list_part_2 = [hard_drive_dir + 'RF_data/DATA_RFAA_part_2/Austria/data_AAA_corrected/',
                         hard_drive_dir + 'RF_data/DATA_RFAA_part_2/CIFAlps/data_YP2012/',
                         hard_drive_dir + 'RF_data/DATA_RFAA_part_2/data_DINAR/',
                         hard_drive_dir + 'RF_data/DATA_RFAA_part_2/HU_SK/data/']
# DATA_RFAA_part_3
path_wavs_list_part_3 = [hard_drive_dir + 'RF_data/DATA_RFAA_part_3/AARF/DATA_MOBST/data/',
                         hard_drive_dir + 'RF_data/DATA_RFAA_part_3/AARF/DATA_PERMST/data/',
                         hard_drive_dir + 'RF_data/DATA_RFAA_part_3/GERMANY/DE_AA_RF/DATA/data/']
path_wavs_list_part_4 = [hard_drive_dir + 'RF_data/CIFALPS/data_YP2012/']
# INGV
path_wavs_list_part_5 = [hard_drive_dir + 'RF_data/INGV-Permanent-data/',
                         hard_drive_dir + 'RF_data/INGV-Temporary-data/data/']

# Path to store RFs
# path_out_RF = '/media/kmichall/SEISMIC_DATA/RF_calculations/'
path_out_RF = desktop_dir + '/RF_test/RF_1/'

a = RF.calculate_rf(path_ev=path_wavs, path_out=path_out_RF,
                    iterations=200, ds=30, c1=10, c2=10, c3=1, c4=1,
                    max_frequency=1, save=True, plot=True)

