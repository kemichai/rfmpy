"""
Code for plotting receiver functions.

Original code written by Matteo Scarponi on 30.11.2021
Location: Chavannes-pres-renens, CH
Date: Jan 2022
Author: Konstantinos Michailos
"""
import glob
import obspy
import numpy as np
from rfmpy.visuals import plotting as plt_rf
from rfmpy.visuals import tools
import matplotlib.pyplot as plt
from obspy.taup import TauPyModel
import platform
from obspy.geodetics import kilometers2degrees
import matplotlib
import rfmpy.utils.RF_Util as rf_util


###########################################
# Define list of stations for reading RFs #
###########################################

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

# Set figure details
font = {'family': 'normal',
        'weight': 'normal',
        'size': 18}
matplotlib.rc('font', **font)
# Set figure width to 12 and height to 9
fig_size = plt.rcParams["figure.figsize"]
fig_size[1] = 10.27
fig_size[0] = 13.69
plt.rcParams["figure.figsize"] = fig_size
# Path at which receiver functions are stored
pathRF = '/media/kmichall/SEISMIC_DATA/RF_calculations/RF/'

# Read all the available RFs and create a list of all the stations
# that have RF calculated
path_wavs_list_part = [pathRF]
sta = rf_util.get_station_info(path_wavs_list_part)
unique_all_sta = []
for s in sta:
    if s not in unique_all_sta:
        unique_all_sta.append(s)

# TODO: sort list and read four by four as below...

stations = ["PANIX", "ROTHE","SIMPL", "LAUCH"]  # !!! Based on number of stations check out...
# ... the number of subplots you want to have: LINE 139


# Compute a reference ray parameter (pref) for normal moveout correction
model = TauPyModel(model="iasp91")
km_to_deg = 111.19
pref = kilometers2degrees((model.get_travel_times(source_depth_in_km=0,
                                                  distance_in_degree=65,
                                                  phase_list=["P"])[0].ray_param_sec_degree))
# Plotting and moveout correction parameters
bazstart = 10
bazend = 370
bazstep = 10
amplitude = 5

Z, VP, VS = plt_rf.get_iasp91(zmax=200, step=0.25, zmoho=75)


f = plt.figure(1)
# Loop on stations to read RFs
station_number = 0
for station in stations:
    station_number += 1
    files = glob.glob(pathRF + "*" + station + "*")
    if len(files) == 0:
        continue
    stream = obspy.read(pathRF + "*" + station + "*")
    for trace in stream:
        # Compute ray parameter for the single RF in SECONDS PER KM
        trace.prai = (model.get_travel_times(source_depth_in_km=trace.stats.sac.evdp / 1000,
                                             distance_in_degree=trace.stats.sac.dist,
                                             phase_list=["P"],)[0].ray_param_sec_degree/ km_to_deg)  # SECONDS/KM
        trace.stats.baz = trace.stats.sac.baz
        # trace.stats.sac.a = 5

        # Apply normal moveout correction
        # This step is very useful before stacking RFs in time
        # as it homogenizes the RF delay times, accounting for different earthquake sources
        # which happened at different depths and distances with respect to the seismic station
        # You can compare the effect of with/without this correction by commenting line 99
        trace.data = tools.rfmops(trace, pref, Z, VP, VS)
    stream.sort(keys=["baz"])

    # Prepare to plot RFs, STACKED by back-azimuthal direction for each station
    allbaz = np.arange(bazstart, bazend, bazstep)
    stack = np.zeros((len(allbaz), len(stream[0].data)))
    count = np.zeros(len(allbaz))

    ibaz = []
    index = 0
    for baz in range(len(allbaz)):
        counter = 0
        indexes = []
        for trace in stream:
            if abs(trace.stats.baz - allbaz[baz]) <= bazstep / 2.0:
                stack[baz, :] += trace.data[:]
                counter += 1
                count[baz] += 1
                indexes.append(int(index))
            index += 1
        if counter != 0:
            stack[baz, :] = stack[baz, :] / counter
        else:
            stack[baz, :] = np.nan
        ibaz.append(np.squeeze(indexes))

    # Prepare subplots - > number of subplots has to be adjusted based on number of stations!!! #
    # Watch out for paramter trace.stats.sac.a:
    # it stores the delay of the direct-P arrival with respect to the
    # beginning of the trace.
    # This parameter is determined in the first place in the beginning of the RF computation
    # where the Z trace is cropped of "a" seconds to avoid the RF starting from zero.

    ax = plt.subplot(1, len(stations), station_number)
    tt = range(len(stream[0].data)) * trace.stats.sac.delta - trace.stats.sac.a
    # Loop on the back-azimuths and plot the stacked traces for the given interval
    for i in range(len(allbaz)):
        ax.plot(tt, stack[i, :] * amplitude + i, "k", lw=0.5, zorder=-i)
        ax.fill_between(tt, y1=stack[i, :] * amplitude + i, y2=i,
                        where=[stack[i, j] * 10 + i > i for j in range(len(stream[0].data))],
                        color="dodgerblue", zorder=-i,)
        ax.fill_between(tt, y1=stack[i, :] * amplitude + i, y2=i,
                        where=[stack[i, j] * 10 + i < i for j in range(len(stream[0].data))],
                        color="red", zorder=-i,)
    ax.axvline(x=0, ymin=0, ymax=1, color="g", alpha=0.5, linestyle="--", lw=1, zorder=-20)
    # ax.axvline(x=1, ymin=0, ymax=1, color="r", alpha=0.25, linestyle="--", lw=1, zorder=-20)
    # ax.axvline(x=2, ymin=0, ymax=1, color="r", alpha=0.25, linestyle="--", lw=1, zorder=-20)
    print(-i)


    ax.set_axisbelow(True)
    ax.set_yticks(np.arange(len(allbaz)))
    major_ticks = np.arange(0, 25, 5)
    minor_ticks = np.arange(-5, 25, 1)
    ax.set_xticks(major_ticks)
    ax.set_xticks(minor_ticks, minor=True)
    # Or if you want different settings for the grids:
    ax.grid(which='minor', alpha=0.7, linestyle="--", zorder=0)
    ax.grid(which='major', alpha=0.9, zorder=0)
    ax.set_xlim([-5, 25])
    ax.set_ylim([-1, len(allbaz) + 1])
    ax.set_yticklabels(allbaz)
    ax.text(0.975, 0.975, station, transform=ax.transAxes, va="center", ha="right")
    if station_number == 1:
        ax.set_ylabel("Back Azimuth (degrees)")
    else:
        ax.set_yticklabels([])
    ax.set_xlabel("Time (s)")

    # WRITING NUMBER OF TRACES ON THE RIGHT y-AXIS #
    ax = ax.twinx()  # instantiate a second axes that shares the same x-axis
    color = "tab:blue"
    ax.set_ylim([-1, len(allbaz) + 1])
    ax.set_yticks(np.arange(len(allbaz)))
    ax.set_yticklabels(np.array(count, dtype="int"))
    ax.tick_params(axis="y", labelcolor=color)
    # ax.set_ylabel("Number of traces", rotation=90)

    if station_number == 4 or station_number == 10:
        ax.set_ylabel("Number of traces", rotation=90)
    ax.yaxis.label.set_color("tab:blue")

plt.tight_layout()
plt.savefig('/home/kmichall/Desktop/RF_plots/' + station, format='png', dpi=300)
plt.show()
plt.close()
