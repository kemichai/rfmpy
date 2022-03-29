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
import os
import time

# Set up paths
if platform.node().startswith('kmichailos-laptop'):
    data_root_dir = '/media/kmichailos/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/home/kmichailos/Desktop/codes/github'
    desktop_dir = '/home/kmichailos/Desktop'
    hard_drive_dir = '/media/kmichailos/SEISMIC_DATA/'
else:
    data_root_dir = '/media/kmichall/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/home/kmichall/Desktop/Codes/github'
    desktop_dir = '/home/kmichall/Desktop'
    hard_drive_dir = '/media/kmichall/SEISMIC_DATA/'

# Path in which waveforms are stored
# path_wavs = desktop_dir + '/RF_test/test_data/'

path_wavs = [
    # hard_drive_dir + 'RF_data/DATA_RFAA_part_1/SWISS/data/',
             hard_drive_dir + 'RF_data/DATA_RFAA_part_1/EASI/data/',]
    #          hard_drive_dir + 'RF_data/SLOVENIAN/',
    #          hard_drive_dir + 'RF_data/DATA_RFAA_part_1/FRANCE/data_sort/',
    #          hard_drive_dir + 'RF_data/DATA_RFAA_part_1/FRANCE/data/',
    #          hard_drive_dir + 'RF_data/DATA_RFAA_part_1/North_Italy/events_fri_ven/',
    #          hard_drive_dir + 'RF_data/DATA_RFAA_part_2/Austria/data_AAA_corrected/',
    #          hard_drive_dir + 'RF_data/DATA_RFAA_part_2/CIFAlps/data_YP2012/',
    #          hard_drive_dir + 'RF_data/DATA_RFAA_part_2/data_DINAR/',
    #          hard_drive_dir + 'RF_data/DATA_RFAA_part_2/HU_SK/data/',
    #          hard_drive_dir + 'RF_data/DATA_RFAA_part_3/AARF/DATA_MOBST/data/',
    #          hard_drive_dir + 'RF_data/DATA_RFAA_part_3/AARF/DATA_PERMST/data/',
    #          hard_drive_dir + 'RF_data/DATA_RFAA_part_3/GERMANY/DE_AA_RF/DATA/data/',
    #          hard_drive_dir + 'RF_data/CIFALPS/data_YP2012/',
    #          hard_drive_dir + 'RF_data/INGV-Permanent-data/',
    #          hard_drive_dir + 'RF_data/INGV-Temporary-data/data/']



# Path to store RFs
path_out_RF = '/media/kmichall/SEISMIC_DATA/RF_calculations/'
# path_out_RF = desktop_dir + '/RF_test/RF_Km/'

t_beg = time.time()
# Path for StationXML files
work_dir = os.getcwd()
path_meta = work_dir + '/data/metadata/'
try:
    print('>>> Reading inventory...')
    inv = read_inventory(path_meta + '/*.xml')
    print('>>> Read inventory...')
except Exception as e:
    raise type(e)('>>> Move to the top directory of the repository!')

# =================================================================================================================== #
# Define parameters for calculating receiver functions
# Define sta/lta parameters
sta_lta_qc_parameters = {'sta': 3, 'lta': 50, 'high_cut': 1.0, 'threshold': 2.5}
# Define pre-processing parameters
pre_processing_parameters = {'low_cut': 0.05, 'high_cut': 1.0, 'order': 2, 't_before': 40, 't_after': 60}
for path_wav in path_wavs:
    print(path_wav)
    RF.calculate_rf(path_ev=path_wav, path_out=path_out_RF,
                inventory=inv, iterations=200, ds=30,
                c1=10, c2=10, c3=1, c4=1,
                sta_lta_qc=sta_lta_qc_parameters,
                pre_processing=pre_processing_parameters,
                max_frequency=1, save=True, plot=False)
# =================================================================================================================== #

t_end = time.time()
total_time = t_end - t_beg
print('It took ' + str(round(total_time)) + ' seconds in total.')
