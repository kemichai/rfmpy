import glob
import obspy
import numpy as np
from Tools import tools
import matplotlib.pyplot as plt
from obspy.taup import TauPyModel

###########################################
# Define list of stations for reading RFs #
###########################################

stations = ["IA02A", "IA06A"]  # !!! Based on number of stations check out...
# ... the number of subplots you want to have: LINE 139
#############################
# Path where RFs are stored #
#############################

pathRF = "/Users/mscarpon/Desktop/Projects/Ivrea/Seismic/SAC/01/RF/RRF/RRF_10_10_1_1/"

##########################################################################
# Compute a reference ray parameter (pref) for normal moveout correction #
##########################################################################

model = TauPyModel(model="iasp91")

km_to_deg = 111.19

pref = (
    model.get_travel_times(
        source_depth_in_km=0, distance_in_degree=65, phase_list=["P"]
    )[0].ray_param_sec_degree
    / km_to_deg
)

##############################################
# Plotting and moveout correction parameters #
##############################################

bazstart = 10
bazend = 370
bazstep = 20
amplitude = 4.5

inc = 0.25
zmax = 200
zMoho = 75

Z = np.arange(0, zmax + inc, inc)
VP = tools.pvelIASP(Z, zMoho)
VS = tools.svelIASP(Z, zMoho)

##########
##########
# FIGURE #
##########
##########

f = plt.figure(1, figsize=[8, 6])

################################
# Loop on stations to read RFs #
################################

station_number = 0
for station in stations:
    station_number += 1
    files = glob.glob(pathRF + "*" + station + "*")
    if len(files) == 0:
        continue

    stream = obspy.read(pathRF + "*" + station + "*")

    for trace in stream:

        #############################################################
        # Compute ray parameter for the single RF in SECONDS PER KM #
        #############################################################

        trace.prai = (
            model.get_travel_times(
                source_depth_in_km=trace.stats.sac.evdp / 1000,
                distance_in_degree=trace.stats.sac.dist,
                phase_list=["P"],
            )[0].ray_param_sec_degree
            / km_to_deg
        )  # SECONDS PER KM
        trace.stats.baz = trace.stats.sac.baz

        ###################################
        # Apply normal moveout correction #
        ###################################

        # This step is very useful before stacking RFs in time
        # as it homogenizes the RF delay times, accounting for different earthquake sources
        # which happened at different depths and distances with respect to the seismic station

        # You can compare the effect of with/without this correction by commenting line 99

        trace.data = tools.rfmops(trace, pref, Z, VP, VS)

    stream.sort(keys=["baz"])

    #############################################################################
    # Prepare to plot RFs, STACKED by back-azimuthal direction for each station #
    #############################################################################

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

    #############################################################################################
    # Prepare subplots - > number of subplots has to be adjusted based on number of stations!!! #
    #############################################################################################

    # Watch out for paramter trace.stats.sac.a:
    # it stores the delay of the direct-P arrival with respect to the
    # beginning of the trace.
    # This parameter is determined in the first place in the beginning of the RF computation
    # where the Z trace is cropped of "a" seconds to avoid the RF starting from zero.

    ax = plt.subplot(1, 2, station_number)

    tt = range(len(stream[0].data)) * trace.stats.sac.delta - trace.stats.sac.a

    # Loop on the back-azimuths and plot the stacked traces for the given interval

    for i in range(len(allbaz)):
        ax.plot(tt, stack[i, :] * amplitude + i, "k", lw=0.5, zorder=-i)
        ax.fill_between(
            tt,
            y1=stack[i, :] * amplitude + i,
            y2=i,
            where=[stack[i, j] * 10 + i > i for j in range(len(stream[0].data))],
            color="r",
            zorder=-i,
        )
        ax.fill_between(
            tt,
            y1=stack[i, :] * amplitude + i,
            y2=i,
            where=[stack[i, j] * 10 + i < i for j in range(len(stream[0].data))],
            color="b",
            zorder=-i,
        )

    ax.axvline(
        x=0, ymin=0, ymax=1, color="g", alpha=0.5, linestyle="--", lw=1, zorder=-20
    )
    ax.axvline(
        x=1, ymin=0, ymax=1, color="r", alpha=0.25, linestyle="--", lw=1, zorder=-20
    )
    ax.axvline(
        x=2, ymin=0, ymax=1, color="r", alpha=0.25, linestyle="--", lw=1, zorder=-20
    )

    ax.grid(True, alpha=0.5, lw=0.4, linestyle="--")
    ax.set_xticks([0, 5, 10, 15, 20, 30])
    ax.set_yticks(np.arange(len(allbaz)))
    ax.set_xlim([-2, 20])
    ax.set_ylim([-1, len(allbaz) + 1])
    ax.set_yticklabels(allbaz)
    ax.text(0.975, 0.975, station, transform=ax.transAxes, va="center", ha="right")

    if station_number == 1:
        ax.set_ylabel("Baz (deg)")
    else:
        ax.set_yticklabels([])
    ax.set_xlabel("Time (s)")

    ################################################
    # WRITING NUMBER OF TRACES ON THE RIGHT y-AXIS #
    ################################################

    ax = ax.twinx()  # instantiate a second axes that shares the same x-axis

    color = "tab:blue"
    ax.set_ylim([-1, len(allbaz) + 1])
    ax.set_yticks(np.arange(len(allbaz)))
    ax.set_yticklabels(np.array(count, dtype="int"))
    ax.tick_params(axis="y", labelcolor=color)

    if station_number == 2 or station_number == 10:
        ax.set_ylabel("Number of traces", rotation=90)
    ax.yaxis.label.set_color("tab:blue")

plt.tight_layout()
plt.show()
plt.close()