"""
Read RFs to exclude traces from the time-to-depth migration calculations.
Discard any traces with low quality results (noisy).

Location: Chavannes-pres-renens, CH
Date: JUL 2022
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
pathRF = work_dir + "/data/RF/RF/"
pathRF = '/media/kmichailos/SEISMIC_DATA/RF_calculations/RF/'

# Read all the available RFs and create a list of all the stations
# that have RF calculated
path_wavs_list_part = [pathRF]
sta = rf_util.get_station_info(path_wavs_list_part)
unique_all_sta = []
for s in sta:
    if s.split(' ')[0] not in unique_all_sta:
        unique_all_sta.append(s.split(' ')[0])
unique_all_sta.sort()

# Compute a reference ray parameter (pref) for normal moveout correction
model = TauPyModel(model="iasp91")
km_to_deg = 111.19
pref = kilometers2degrees((model.get_travel_times(source_depth_in_km=0,
                                                  distance_in_degree=65,
                                                  phase_list=["P"])[0].ray_param_sec_degree))
# Plotting and moveout correction parameters
bazstart = 20
bazend = 380
bazstep = 20
amplitude = 2.5

Z, VP, VS = plt_rf.get_iasp91(zmax=200, step=0.25, zmoho=75)

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

list1 = rms
list2 = sta_name
list1, list2 = zip(*sorted(zip(list1, list2)))

# Set figure details
font = {'family': 'normal',
        'weight': 'normal',
        'size': 5}
matplotlib.rc('font', **font)
# Set figure width to 12 and height to 9
fig_size = plt.rcParams["figure.figsize"]
fig_size[1] = 7
fig_size[0] = 10


plt.bar(list2, list1)
plt.xlabel('Seismic stations')
plt.ylabel('Average RMS values')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig('RMS_values.png', format='png', dpi=300)
plt.show()


for i, rms_ in enumerate(list1):
    print(rms_, list2[i])


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

bins = np.arange(0.0, 0.2, 0.001)
# plt.hist(all_traces_rms, bins, histtype='step', orientation='vertical',
#              color='gray',facecolor='gray', alpha=0.7, linewidth=1.5,
#              edgecolor='k',fill=True)
plt.bar(index, all_traces_rms)
plt.title('rms from individual traces')
plt.show()