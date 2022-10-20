"""
Read RFs to exclude traces from the time-to-depth migration calculations.
Discard any traces with low quality results (noisy).

Code consists of two main parts:
1) move 'bad' RFs to another directory so we excluded them from the time-to-depth migration
2) plot the rms values of all RFs to choose the threshold (in this case we use 0.07)

Location: Chavannes-pres-renens, CH
Date: Jul 2022
Author: Konstantinos Michailos
"""
import platform
import rfmpy.core.migration_sphr as rf_mig
import os
import rfmpy.utils.RF_Util as rf_util
from obspy.taup import TauPyModel
from obspy.geodetics import kilometers2degrees
from rfmpy.visualisation import plotting as plt_rf
import glob
import obspy
import numpy as np
from rfmpy.visualisation import plotting as plt_rf
from rfmpy.visualisation import tools
import matplotlib.pyplot as plt
from obspy.taup import TauPyModel
import platform
from obspy.geodetics import kilometers2degrees
import matplotlib
import rfmpy.utils.RF_Util as rf_util
from obspy import Stream
import shutil

# Set up paths
if platform.node().startswith('kmichailos-laptop'):
    data_root_dir = '/media/kmichailos/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/home/kmichailos/Desktop/codes/github'
    desktop_dir = '/home/kmichailos/Desktop'
    hard_drive_dir = '/media/kmichailos/SEISMIC_DATA/'
else:
    data_root_dir = '/media/kmichall/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/github'
    desktop_dir = '/home/kmichall/Desktop'
    hard_drive_dir = '/media/kmichall/SEISMIC_DATA/'

# Define paths
work_dir = os.getcwd()
# pathRF = work_dir + "/data/RF/RF/"
# pathRF = '/media/kmichailos/SEISMIC_DATA/RF_calculations/RF/'
pathRF = '/home/kmichailos/Desktop/all_rfs/RF/'
pathTRF = '/home/kmichailos/Desktop/all_rfs/TRF/'
# path_badRF = '/media/kmichailos/SEISMIC_DATA/RF_calculations/RF_low_quality/'
path_badRF = '/home/kmichailos/Desktop/all_rfs/RF_low_quality/'
path_badTRF = '/home/kmichailos/Desktop/all_rfs/TRF_low_quality/'
path_TRF_ = '/home/kmichailos/Desktop/all_rfs/TRF_/'

all_files = glob.glob(pathRF + "*")
my_final_list = set(all_files)

all_files_trf = glob.glob(pathTRF + "*")
my_final_list_trf = set(all_files_trf)


if not os.path.exists(path_badRF):
    os.mkdir(path_badRF)
    print("Directory '%s' created" % path_badRF.split('/')[-2])
else:
    print("Directory '%s' exists" % path_badRF.split('/')[-2])

# Read all the available RFs and create a list of all the stations that have RF calculated
path_wavs_list_part = [pathRF]
sta = rf_util.get_station_info(path_wavs_list_part)
unique_all_sta = []
for s in sta:
    if s.split(' ')[0] not in unique_all_sta:
        unique_all_sta.append(s.split(' ')[0])
unique_all_sta.sort()

for station in unique_all_sta:
    files = glob.glob(pathRF + "*" + station + "*")
    files_trf = glob.glob(pathTRF + "*" + station + "*")
    if len(files) == 0:
        continue
    for i, file in enumerate(files):
        # print(file)
        tr = obspy.read(file)
        # Compute rms
        tr.rms = np.sqrt(np.mean(np.square(tr[0].data)))
        if tr.rms >= 0.07:
            print('Discarding trace: ', tr, ' as rms value is larger than the threshold.')
            try:
                shutil.move(file, path_badRF)
                shutil.move(files_trf[i], path_badTRF)
            except Exception as e:
                print(e)




##############################
# Plotting part of the code...
km_to_deg = 111.19
model = TauPyModel(model="iasp91")

# Read all the available RFs and create a list of all the stations
# that have RF calculated
path_wavs_list_part = [pathRF]
sta = rf_util.get_station_info(path_wavs_list_part)
unique_all_sta = []
for s in sta:
    if s.split(' ')[0] not in unique_all_sta:
        unique_all_sta.append(s.split(' ')[0])
unique_all_sta.sort()

# Loop on stations to read RFs
station_number = 0
rms = []
sta_name = []
all_traces_rms = []
for station in unique_all_sta:
    station_number += 1
    files = glob.glob(pathRF + "*" + station + "*")
    if len(files) == 0:
        continue
    stream = Stream()
    for file in files:
        tr = obspy.read(file)
        if len(tr[0].data) == 2000:
            stream.append(tr[0])
    stream_rms = []
    for trace in stream:
        # Compute ray parameter for the single RF in SECONDS PER KM
        trace.prai = (model.get_travel_times(source_depth_in_km=trace.stats.sac.evdp / 1000,
                                             distance_in_degree=trace.stats.sac.dist,
                                             phase_list=["P"], )[0].ray_param_sec_degree / km_to_deg)  # SECONDS/KM
        trace.stats.baz = trace.stats.sac.baz
        trace.rms = np.sqrt(np.mean(np.square(trace.data)))
        stream_rms.append(trace.rms)
        all_traces_rms.append(trace.rms)
    print(np.mean(stream_rms), station)
    stream_rms.sort()
    # plt.hist(stream_rms)
    # plt.title(station)
    # plt.show()

    rms.append(np.mean(stream_rms))
    sta_name.append(station)


# Set figure details
font = {'family': 'normal',
        'weight': 'normal',
        'size': 15}
matplotlib.rc('font', **font)
# Set figure width to 12 and height to 9
fig_size = plt.rcParams["figure.figsize"]
fig_size[1] = 7
fig_size[0] = 12
# All rms values from individual traces
all_traces_rms.sort()
index = []
for i in range(len(all_traces_rms)):
    print(i)
    index.append(i)


# TODO: add this in the SI
bins = np.arange(0.0, len(all_traces_rms), 1)
plt.hist(bins, all_traces_rms, histtype='step', orientation='vertical',
             color='gray',facecolor='gray', alpha=0.7, linewidth=1.5,
             edgecolor='k',fill=True, label='RFs')
# plt.hist(bins, all_traces_rms, label='RFs')
plt.xlabel('Index')
plt.ylabel('RMS')
plt.axhline(y=0.07, color='r', linestyle='-', label='Cut-off limit')
plt.tight_layout()
plt.legend()
plt.savefig('RMS_values.png', format='png', dpi=300)

plt.show()