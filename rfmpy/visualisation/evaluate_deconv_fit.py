"""
Sort script to read RFs and make plots for evaluating how the fit
of the deconvolution function looks like in different places and
vs magnitudes, depths and distances from the teleseismic event.

Location: Chavannes-pres-renens, CH
Date: Feb 2022
Author: Konstantinos Michailos
"""

from obspy import read
import platform
import glob
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from scipy.stats import gaussian_kde

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

# Path to store RFs
path_to_RFs = '/media/kmichailos/SEISMIC_DATA/RF_calculations/RF/'
path_to_RFs = desktop_dir + "/all_rfs/RRF/RF/"

all_rf_files = glob.glob(path_to_RFs + '*SAC')

# cross correlation value between R component and convolved one
cc = []
st_lat = []
st_lon = []
ev_dist = []
ev_dep = []
ev_mag = []
sta_name = []
for rf_file in all_rf_files:
    rf = read(rf_file)
    cc.append(rf[0].stats.sac.f)
    st_lat.append(rf[0].stats.sac.stla)
    st_lon.append(rf[0].stats.sac.stlo)
    ev_dist.append(rf[0].stats.sac.dist)
    ev_dep.append(rf[0].stats.sac.evdp)
    ev_mag.append(rf[0].stats.sac.mag)
    sta_name.append(rf[0].stats.station)

# Plotting part of the script

font = {'family': 'normal',
        'weight': 'normal',
        'size': 18}
matplotlib.rc('font', **font)
# Set figure width to 12 and height to 9
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 15
fig_size[1] = 8
plt.rcParams["figure.figsize"] = fig_size
subplot_rect = {'left': 0.08, 'right': 0.92, 'bottom': 0.08, 'top': 0.95, 'wspace': 0.1, 'hspace': 0.1}



#############################
# Calculate the point density
bins = np.arange(27.5, 102.5, 5)
# ................................................
ax1 = plt.subplot2grid((2, 4), (1, 0), colspan=1)
ax1.scatter(ev_dist, cc, facecolor='darkgrey',
            alpha=0.6, edgecolor='black',
            linewidth=0.6)
# ax1.scatter(ev_dist, cc, c='k', alpha=0.1)
ax1.set_ylabel('Cross-correlation fit', fontsize=20)
ax1.set_xlabel('Event-receiver distance (deg)', fontsize=20)
ax1.tick_params(bottom=True, top=True, left=True, right=True)
plt.xticks(np.arange(30, 110, 10))
plt.yticks(np.arange(0.1, 1.01, 0.1))
ax1.grid('True', linestyle="-", color='gray', linewidth=0.1, alpha=0.5)
ax1.set_xlim([25, 100])

ax2 = plt.subplot2grid((2, 4), (0, 0), colspan=1)
ax2.hist(ev_dist, bins, histtype='step', orientation='vertical',
         color='black',facecolor='black', alpha=0.8, linewidth=2, linestyle='-',
         edgecolor='black',fill=False)
ax2.grid('True', linestyle="-", color='gray', linewidth=0.1, alpha=0.5)
plt.xticks(np.arange(30, 110, 10))
ax2.set_ylabel('Number of data', fontsize=20)
ax2.tick_params(bottom=True, top=True, left=True, right=True)
ax2.set_xticklabels([])
ax2.set_xlim([25, 100])
ax2.set_ylim([0, 25000])


bins = np.arange(5.45, 9.05, .1)
# ................................................
ax3 = plt.subplot2grid((2, 4), (1, 1), colspan=1)
ax3.scatter(ev_mag, cc, facecolor='darkgrey',
            alpha=0.6, edgecolor='black',
            linewidth=0.6)
ax3.set_xlabel('Magnitude', fontsize=20)
plt.yticks(np.arange(0.1, 1.01, 0.1))
plt.xticks(np.arange(5.5, 9, 0.5))
ax3.set_yticklabels([])
ax3.tick_params(bottom=True, top=True, left=True, right=True)
ax3.grid('True', linestyle="-", color='gray', linewidth=0.1, alpha=0.5)
ax3.tick_params(bottom=True, top=True, left=True, right=True)
ax3.set_xlim([5.3, 8.5])

ax4 = plt.subplot2grid((2, 4), (0, 1), colspan=1)
ax4.hist(ev_mag, bins, histtype='step', orientation='vertical',
         color='black',facecolor='black', alpha=0.8, linewidth=2, linestyle='-',
         edgecolor='black',fill=False)
ax4.grid('True', linestyle="-", color='gray', linewidth=0.1, alpha=0.5)
plt.xticks(np.arange(5.5, 9, 0.5))
ax4.set_xticklabels([])
ax4.set_xlim([5.3, 8.5])
ax4.tick_params(bottom=True, top=True, left=True, right=True)
ax4.set_ylim([0, 25000])
ax4.set_yticklabels([])


bins = np.arange(-5, 705.0, 10)
# ................................................
ax5 = plt.subplot2grid((2, 4), (1, 2), colspan=1)
ax5.scatter(ev_dep, cc, facecolor='darkgrey',
            alpha=0.6, edgecolor='black',
            linewidth=0.6)
ax5.set_xlabel('Depth (km)', fontsize=20)
ax5.set_yticklabels([])
plt.yticks(np.arange(0.1, 1.01, 0.1))
plt.xticks(np.arange(0.0, 700, 100))
ax5.tick_params(bottom=True, top=True, left=True, right=True)
ax5.grid('True', linestyle="-", color='gray', linewidth=0.1, alpha=0.5)
ax5.tick_params(bottom=True, top=True, left=True, right=True)
ax5.set_xlim([-20, 710])

ax6 = plt.subplot2grid((2, 4), (0, 2), colspan=1)
ax6.hist(ev_dep, bins, histtype='step', orientation='vertical',
         color='black',facecolor='black', alpha=0.8, linewidth=2, linestyle='-',
         edgecolor='black',fill=False)
ax6.grid('True', linestyle="-", color='gray', linewidth=0.1, alpha=0.5)
plt.xticks(np.arange(0, 700, 100))
ax6.set_xticklabels([])
ax6.set_xlim([-20, 710])
ax6.tick_params(bottom=True, top=True, left=True, right=True)
ax6.set_ylim([0, 25000])
ax6.set_yticklabels([])


bins = np.arange(0.1005, 1.005, 0.01)
ax2 = plt.subplot2grid((2, 4), (1, 3), colspan=1, rowspan=1)
ax2.hist(cc, bins, histtype='step', orientation='horizontal',
         color='black',facecolor='black', alpha=0.8, linewidth=2, linestyle='-',
         edgecolor='black',fill=False)
ax2.grid('True', linestyle="-", color='gray', linewidth=0.1, alpha=0.5)
# plt.xticks(np.arange(30, 110, 10))
ax2.set_xlabel('Number of data', fontsize=20)
ax2.tick_params(bottom=True, top=True, left=True, right=True)
ax2.set_yticklabels([])
ax2.set_ylim([0.1, 1.035])
# ax2.yaxis.tick_right()
# plt.yticks(np.arange(0.1, 1.01, 0.1))
ax2.set_yticks(ax5.get_yticks()[1:])
plt.yticks(np.arange(0.1, 1.01, 0.1))



plt.tight_layout()
plt.savefig('cc_vs_dep_mag_dist.png', bbox_inches="tight", format='png', dpi=300)
plt.show()


# Map...
# plot map with colored stations and cc

# calculate average cc for each station
unique_stations = []
lat = []
lon = []
for i, station in enumerate(sta_name):
    if station not in unique_stations:
        unique_stations.append(station)
        lat.append(st_lat[i])
        lon.append(st_lon[i])


ccs = [[] for i in range(len(unique_stations))]
for i, sta in enumerate(unique_stations):
    for j, sta_ in enumerate(sta_name):
        print(j)
        if sta == sta_:
            print('same station')
            ccs[i].append(cc[j])

cc_ave = []
for i, sta in enumerate(unique_stations):
    print(i, np.average(ccs[i]))
    cc_ave.append(np.average(ccs[i]))

# Write file
for i, ln in enumerate(lon):
    with open('cc_values.txt', 'a') as of:
        of.write('{}, {}, {}\n'.format(lon[i], lat[i], cc_ave[i]))



from mpl_toolkits.basemap import Basemap


MIN_LAT, MAX_LAT, MIN_LON, MAX_LON = (40, 52, 0, 20)

font = {'family': 'normal',
        'weight': 'normal',
        'size': 18}
matplotlib.rc('font', **font)
# Set figure width to 12 and height to 9
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 9
fig_size[1] = 7
plt.rcParams["figure.figsize"] = fig_size
subplot_rect = {'left': 0.08, 'right': 0.92, 'bottom': 0.08, 'top': 0.95, 'wspace': 0.1, 'hspace': 0.1}


fig = plt.figure()
# Map view
ax1 = plt.subplot2grid((3, 3), (0, 0), colspan=3, rowspan=3)
# fig.sca(ax1)
# Define Basemap
# Note: You can define the resolution of the map you just created. Higher
# resolutions take longer to create.
#    'c' - crude
#    'l' - low
#    'i' - intermediate
#    'h' - high
#    'f' - full
bmap = Basemap(llcrnrlon=MIN_LON, llcrnrlat=MIN_LAT, urcrnrlon=MAX_LON,
               urcrnrlat=MAX_LAT, resolution='l', projection='merc',
               lat_0=MIN_LAT, lon_0=MIN_LON, ax=ax1, fix_aspect=False)
# Draw some map elements on the map
bmap.drawcoastlines()
bmap.drawmapboundary(fill_color='white')
bmap.drawcountries()
bmap.fillcontinents(color='white',lake_color='white')
bmap.drawparallels(np.arange(MIN_LAT, MAX_LAT, 2), labels=[1, 0, 0, 0],
                   linewidth=0.1, dashes=[1, 5])
bmap.drawmeridians(np.arange(MIN_LON, MAX_LON, 2), labels=[0, 0, 0, 1],
                   linewidth=0.1, dashes=[1, 5])
xx, yy = bmap(lon, lat)
conf = bmap.scatter(xx, yy, edgecolor="k", alpha=0.95, s=180, marker='v',
             label='Tibet 1D', c=cc_ave, cmap='viridis', ax=ax1, zorder=2)
bmap.drawmapscale(
    2, 46.5, 2, 46.5,
    200,
    units='km', fontsize=10,
    yoffset=None,
    barstyle='simple', labelstyle='simple',
    fillcolor1='w', fillcolor2='#000000',
    fontcolor='#000000',
    zorder=101)
cax = fig.add_axes([0.15, 0.7, 0.2, 0.02])
cb = fig.colorbar(conf, cax=cax, orientation="horizontal",
                  extend='both',
                  ticks=[0.7, 0.8, 0.9, 1.0])
cb.set_label("Average CC values (fit)", fontsize=16)
# plt.subplots_adjust(**subplot_rect)
plt.tight_layout()
plt.savefig('cc_values_map.png', bbox_inches="tight", format='png', dpi=300)
plt.show()

