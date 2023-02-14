"""
Calculate residuals between our moho estimates and from previous studies
such as Spada et al. or Grad and Tira.

KM = Michailos et al (our study)
SP = Spada et al
GD = Grad and Tira

Location: Chavannes-pres-renens, CH
Date: Feb 2023
Author: Konstantinos Michailos
"""

import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Qt5Agg")


def dist_calc(loc1, loc2):
    """..."""
    r = 6371.009  # Radius of the Earth in km
    dlat = np.radians(abs(loc1[0] - loc2[0]))
    dlong = np.radians(abs(loc1[1] - loc2[1]))
    ddepth = abs(loc1[2] - loc2[2])
    mean_lat = np.radians((loc1[0] + loc2[0]) / 2)
    dist_ = r * np.sqrt(dlat ** 2 + (np.cos(mean_lat) * dlong) ** 2)
    dist = np.sqrt(dist_ ** 2 + ddepth ** 2)
    return dist

def map_absolute_differences(longitudes, latitudes, longitudes_, latitudes_, depth_differences):
    """
    Plots a map of absolute differences in depth between two sets of locations.

    :param longitudes: A 1D numpy array of longitudes for the first set of locations.
    :type longitudes: numpy.ndarray
    :param latitudes: A 1D numpy array of latitudes for the first set of locations.
    :type latitudes: numpy.ndarray
    :param longitudes_: A 1D numpy array of longitudes for the second set of locations.
    :type longitudes_: numpy.ndarray
    :param latitudes_: A 1D numpy array of latitudes for the second set of locations.
    :type latitudes_: numpy.ndarray
    :param depth_differences: A 1D numpy array of depth differences between the two sets of locations.
    :type depth_differences: numpy.ndarray

    :return: A Matplotlib figure object with a map of the absolute depth differences.
    :rtype: matplotlib.figure.Figure
    """

    x_obs = np.asarray(longitudes)
    y_obs = np.asarray(latitudes)
    x = np.asarray(longitudes_)
    y = np.asarray(latitudes_)
    v = np.asarray(depth_differences)

    font = {'family': 'normal',
            'weight': 'normal',
            'size': 14}
    matplotlib.rc('font', **font)
    # Set figure width to 12 and height to 9
    fig_size = plt.rcParams["figure.figsize"]
    fig_size[0] = 12
    fig_size[1] = 15
    plt.rcParams["figure.figsize"] = fig_size

    min_lat, max_lat, min_lon, max_lon = (41, 52, 0, 25)

    fig = plt.figure()
    ax1 = plt.subplot2grid((2,2), (0, 0), colspan=2, rowspan=2)
    bmap = Basemap(llcrnrlon=min_lon, llcrnrlat=min_lat, urcrnrlon=max_lon,
                   urcrnrlat=max_lat, resolution='l', projection='merc',
                   lat_0=min_lat, lon_0=min_lon, ax=ax1, fix_aspect=False)
    bmap.drawcoastlines()
    bmap.drawmapboundary(fill_color='white')
    bmap.fillcontinents(color='whitesmoke', lake_color='white')
    bmap.drawparallels(np.arange(min_lat, max_lat, 2), labels=[1, 0, 0, 0],
                       linewidth=0.5, dashes=[1, 10])
    bmap.drawmeridians(np.arange(min_lon, max_lon, 2), labels=[0, 0, 0, 1],
                       linewidth=0.5, dashes=[1, 10])
    xi, yi = bmap(x, y)
    x_obs_, yxobs_ = bmap(x_obs, y_obs)

    conf = bmap.scatter(xi, yi, c=v, s=150, cmap='coolwarm', ax=ax1,
                        vmin=-20, vmax=20, zorder=101, label=r'Moho - Moho$_{GT}$')
    bmap.scatter(x_obs_, yxobs_, c='black', s=10, zorder=101, label="Location of Moho depths (this study)")
    cax = fig.add_axes([0.75, 0.25, 0.2, 0.02])
    cb = fig.colorbar(conf, cax=cax, orientation="horizontal",
                      extend='both',
                      ticks=[-20, -10, 0, 10, 20])
    cb.set_label("Moho depth residuals (km)", fontsize=16)
    fig.legend(loc='upper right')

    plt.show()


def plot_histo_differences(dep_diff_grad_, dep_diff_spada_):
    """"""
    # Set figure width to 12 and height to 9
    fig_size = plt.rcParams["figure.figsize"]
    fig_size[0] = 12
    fig_size[1] = 9
    plt.rcParams["figure.figsize"] = fig_size
    # Plot histogram of the absolute differences in Moho depths
    bins = np.arange(-25.5, 25.5, 1)
    ax2 = plt.subplot2grid((1, 1), (0, 0), colspan=1, rowspan=1)
    plt.axvspan(0, 30, facecolor='gray',alpha=0.09)
    # plt.axvline(x=0, c='black', linewidth=0.5)
    ax2.hist(dep_diff_grad_, bins, histtype='step', orientation='vertical',
             color='black',facecolor='black', alpha=0.8, linewidth=2.,
             linestyle='-', edgecolor='black',fill=False, label=r'Moho - Moho$_{GT}$')
    ax2.hist(dep_diff_spada_, bins, histtype='step', orientation='vertical',
             color='black',facecolor='black', alpha=0.8,
             linewidth=2., linestyle='--', edgecolor='black',
             fill=False, label=r'Moho - Moho$_{SP}$')
    # plt.axvline(np.mean(dep_diff_spada_), color='None',
    #             linewidth=2, label='Median difference with Spada et al. =' +
    #             str(round(np.mean(dep_diff_spada_), 1)))

    ax2.grid('True', linestyle=":", color='gray', linewidth=0.1, alpha=0.4)
    ax2.set_ylabel('Number of observations', fontsize=20)
    ax2.set_xlabel('Moho depth residuals (km)', fontsize=20)
    ax2.tick_params(bottom=True, top=True, left=True, right=True)
    ax2.legend(loc='best')
    ax2.set_xlim([-25, 25])
    plt.text(.75, 0.98, 'Moho estimates from this study are relatively deeper',
         horizontalalignment='center',style='italic',
         verticalalignment='center',
         transform = ax2.transAxes)
    plt.tight_layout()
    plt.show()

# Read files from KM, SP, GD
lons_KM = []
lats_KM = []
deps_KM = []
with open("moho_depths_all.dat", "r") as f:
    for line in f:
        ln = line.split(',')
        lon = float(ln[0])
        lat = float(ln[1])
        mh = float(ln[2])
        lons_KM.append(lon)
        lats_KM.append(lat)
        deps_KM.append(mh)

# Spada et al
lons_SP = []
lats_SP = []
deps_SP = []
with open("Spada_moho.dat", "r") as f:
    for line in f:
        ln = line.split(' ')
        lon = float(ln[0])
        lat = float(ln[1])
        mh = float(ln[2])
        lons_SP.append(lon)
        lats_SP.append(lat)
        deps_SP.append(mh)

# Spada et al
lons_GD = []
lats_GD = []
deps_GD = []
with open("grad.csv", "r") as f:
    for line in f:
        ln = line.split(',')
        lon = float(ln[0])
        lat = float(ln[1])
        mh = float(ln[2])
        if lon > 0 and lon < 25 and 55 > lat > 40:
            lons_GD.append(lon)
            lats_GD.append(lat)
            deps_GD.append(mh)


loc_array_SP = []
for i, j in enumerate(lons_SP):
    loc_array_SP.append([lons_SP[i], lats_SP[i], deps_SP[i]])
loc_array_GD = []
for i, j in enumerate(lons_GD):
    loc_array_GD.append([lons_GD[i], lats_GD[i], deps_GD[i]])

plt.scatter(lons_SP, lats_SP, c=deps_SP, cmap='viridis',
            s=50, edgecolors="None",
            marker="v", label="Spada et al.", )
plt.scatter(lons_KM, lats_KM, c=deps_KM, cmap='viridis',
            s=50, facecolors="grey", edgecolors="k",
            marker="o", label="Our study", )
plt.colorbar()
plt.legend()
plt.show()

plt.scatter(lons_GD, lats_GD, c=deps_GD, cmap='viridis',
            s=50, edgecolors="None",
            marker="v", label="Grad and Tira", )
plt.scatter(lons_KM, lats_KM, c=deps_KM, cmap='viridis',
            s=50, facecolors="grey", edgecolors="k",
            marker="o", label="Our study", )
plt.colorbar()
plt.legend()
plt.show()


lon_diff_SP = []
lat_diff_SP = []
dep_diff_SP = []
radius = 10
for i, ln in enumerate(lons_KM):
    # set depth to zero to find the distance from the surface
    location = [ln, lats_KM[i], 0]
    for j, sp_ln in enumerate(lons_SP):
        SP_location = [sp_ln, lats_SP[j], 0]
        dist = dist_calc(location, SP_location)
        if dist < radius:
            print(f"Nearby observations:  {deps_KM[i]}, {deps_SP[j]} !")
            depth_difference = deps_KM[i] - deps_SP[j]
            lon_diff_SP.append(sp_ln)
            lat_diff_SP.append(lats_SP[j])
            dep_diff_SP.append(depth_difference)

lon_diff_GD = []
lat_diff_GD = []
dep_diff_GD = []
for i, ln in enumerate(lons_KM):
    # set depth to zero to find the distance from the surface
    location = [ln, lats_KM[i], 0]
    for j, gd_ln in enumerate(lons_GD):
        GD_location = [gd_ln, lats_GD[j], 0]
        dist = dist_calc(location, GD_location)
        if dist < radius:
            print(f"Nearby observations:  {deps_KM[i]}, {deps_GD[j]} !")
            depth_difference = deps_KM[i] - deps_GD[j]
            lon_diff_GD.append(gd_ln)
            lat_diff_GD.append(lats_GD[j])
            dep_diff_GD.append(depth_difference)

# for i, ln in enumerate(lon_diff_SP):
#     with open('diff2spada.txt', 'a') as of:
#         of.write('{} {} {}\n'.format(lon_diff_SP[i], lat_diff_SP[i], dep_diff_SP[i]))
#
# for i, ln in enumerate(lon_diff_GD):
#     with open('diff2Grad.txt', 'a') as of:
#         of.write('{} {} {}\n'.format(lon_diff_GD[i], lat_diff_GD[i], dep_diff_GD[i]))


spada = map_absolute_differences(lons_KM, lats_KM, lon_diff_SP, lat_diff_SP, dep_diff_SP)
grad = map_absolute_differences(lons_KM, lats_KM, lon_diff_GD, lat_diff_GD, dep_diff_GD)
plot_histo_differences(dep_diff_GD, dep_diff_SP)

print(f"Spada (Mean: {round(np.mean(dep_diff_SP), 1)} ± {round(np.std(dep_diff_SP), 1)} km, Median: {round(np.median(dep_diff_SP), 1)} km)")
print(f"Grad (Mean: {round(np.mean(dep_diff_GD), 1)} ± {round(np.std(dep_diff_GD), 1)} km, Median: {round(np.median(dep_diff_GD), 1)} km)")

# Spada (Mean: 1.2 ± 8.4 km, Median: 0.5 km)
# Grad (Mean: 0.1 ± 6.9 km, Median: -1.1 km)
