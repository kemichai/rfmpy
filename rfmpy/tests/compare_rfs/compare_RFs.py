"""
Compare RFs with other codes for the same event.

NOTE: code only works running from the same directory as this python file is stored!

Location: Chavannes-pres-renens, CH
Date: Feb 2022
Author: Konstantinos Michailos
"""

import matplotlib.pyplot as plt
from obspy import read
import numpy as np
import matplotlib
import os
font = {'family': 'normal',
        'weight': 'normal',
        'size': 18}
matplotlib.rc('font', **font)
# Set figure width to 12 and height to 9
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 13
fig_size[1] = 8
plt.rcParams["figure.figsize"] = fig_size
subplot_rect = {'left': 0.08, 'right': 0.92, 'bottom': 0.08, 'top': 0.95, 'wspace': 0.1, 'hspace': 0.1}

# Define working directory and read all the data
work_dir = os.getcwd()
try:
    print('>>> Reading data...')
    rf_km_path = work_dir + '/wavs/RF_Km/RF/'
    st_km = read(rf_km_path + '*RRF.SAC')
    rf_gh_path = work_dir + '/wavs/RF_gh/'
    st_gh = read(rf_gh_path + '*.SAC')
    rf_js_path = work_dir + '/wavs/RF_js/'
    st_js = read(rf_js_path + '*.SAC')
except Exception as e:
    raise type(e)('>>> CD to the same directory as the python file! TYPE: cd rfmpy/tests/compare_rfs')


for tr_km in st_km:
    sta_name_km = tr_km.stats.station
    for tr_gh in st_gh:
        sta_name_gh = tr_gh.stats.station
        for tr_js in st_js:
            # Match the same station to plot
            if tr_js.stats.station == sta_name_gh and tr_js.stats.station == sta_name_km:
                tt = np.arange(len(tr_km.data))/20
                tt_js = np.arange(len(tr_js.data))/10 + 25
                tt_gh = np.arange(len(tr_gh.data))/20
                fig = plt.figure()
                ax = fig.add_subplot(1, 1, 1)
                ax.set_title(sta_name_km, fontsize=18)
                ax.axvline(30, color="gray",linewidth=1., zorder=3, linestyle=":", alpha=0.5, label='P-wave arrival')
                ax.plot(tt, tr_km.data, color='black', linestyle='-', label='Python - rfmpy',lw=1.8,zorder=4, alpha=0.7)
                ax.plot(tt_gh, tr_gh.data, color='steelblue', linestyle="--", label='Matlab - GH',lw=1.8,zorder=5)
                ax.plot(tt_js, tr_js.data, color='orangered', linestyle="-.", label='Matlab - JS',lw=1.8,zorder=6)
                ax.set_xlim(25, 50)
                ax.legend(loc='best', fontsize=14)
                plt.savefig(work_dir + '/plots/' + sta_name_km + '.png', bbox_inches="tight", format='png', dpi=300)
                plt.show()
