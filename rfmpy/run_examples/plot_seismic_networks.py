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
# List of unique seismic sites
sta = rf_util.get_station_info(path_wavs_list_part1)
