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

Location: Chavannes-pres-renens, CH
Date: Jan 2022
Author: Konstantinos Michailos
"""

import rfmpy.core.RF_Main as RF
import platform
from obspy import read_inventory, read_events
import os
import time

# Set up paths
if platform.node().startswith('kmichailos-laptop'):
    data_root_dir = '/media/kmichailos/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/home/kmichailos/Desktop/codes/github'
    desktop_dir = '/home/kmichailos/Desktop/'
    hard_drive_dir = '/media/kmichailos/SEISMIC_DATA/'
else:
    data_root_dir = '/home/konstantinos/Desktop/ZNE_waveforms/'
    codes_root_dir = '/home/konstantinos/Desktop/codes/'
    desktop_dir = '/home/konstantinos/Desktop'
    hard_drive_dir = '/home/konstantinos/Desktop/ZNE_waveforms/'

# Path in which waveforms are stored
path_wavs = [
             hard_drive_dir + 'SWISS/',
             hard_drive_dir + 'EASI/',
             hard_drive_dir + 'SLOVENIA/',
             hard_drive_dir + 'OBS/',
             hard_drive_dir + 'FRANCE/',
             hard_drive_dir + 'North_ITALY/',
             hard_drive_dir + 'AUSTRIA/',
             hard_drive_dir + 'DINAR/',
             hard_drive_dir + 'HU_SK/',
             hard_drive_dir + 'MOBST/',
             hard_drive_dir + 'PERMST/',
             hard_drive_dir + 'GERMANY/',
             hard_drive_dir + 'CIFALPS/',
             hard_drive_dir + 'INGV/',
             hard_drive_dir + 'AAPA/',]

# Define paths
work_dir = os.getcwd()
# path_wavs = [work_dir + "/data/data_sample/"]

# Path to store RFs
# path_out_RF = '/media/kmichailos/SEISMIC_DATA/RF_calculations/'
# path_out_RF = work_dir + '/data/RF/'
path_out_RF = '/home/konstantinos/Desktop/all_rfs/'
t_beg = time.time()

# Path for StationXML files
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
                c1=10, c2=10,
                sta_lta_qc=sta_lta_qc_parameters,
                pre_processing=pre_processing_parameters,
                max_frequency=1, save=True, plot=False)
# =================================================================================================================== #

total_time = time.time() - t_beg
print('It took ' + str(round(total_time)) + ' seconds in total.')
