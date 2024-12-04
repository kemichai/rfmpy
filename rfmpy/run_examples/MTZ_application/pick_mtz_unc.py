"""
Pick Moho depths...

Location: Chavannes-pres-renens, CH
Date: Aug 2022
Author: Konstantinos Michailos
"""
import matplotlib
matplotlib.use('TkAgg')
import migration_sphr as rf_mig
import migration_plots_spher as plot_migration_sphr
import numpy as np
import platform
import os
import time
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from collections import OrderedDict
import matplotlib.patches as patches
from scipy.interpolate import RegularGridInterpolator
from pyproj import Geod
import pyproj
from obspy.geodetics import degrees2kilometers, kilometers2degrees


def add_legend(ax, fontsize=12):
    fontsize = 10
    h, l = ax.get_legend_handles_labels()
    by_l = OrderedDict(zip(l, h))
    legend = ax.legend(by_l.values(), by_l.keys(), loc="best", fontsize=fontsize)
    frame = legend.get_frame()
    frame.set_facecolor("white")
    frame.set_edgecolor("black")
    frame.set_linewidth(0.5)
    frame.set_alpha(1)

    return legend


def add_colorbar(ax, m, title=False, ticks=False, ticks_Flag=False, visible=True):
    # colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.05)
    if ticks_Flag:
        cbar = plt.colorbar(m, cax=cax, ticks=ticks)
    else:
        cbar = plt.colorbar(m, cax=cax)
    cbar.ax.set_visible(visible)
    if title:
        cbar.ax.set_title(title, fontsize=fontsize)
    cbar.ax.tick_params(axis="both", which="major", labelsize=fontsize)
    return cbar


fontsize = 12
markersize = 100

def get_end_point(lat1, lon1, baz, d):
    """
    Calculates the end point in lon, lat given we know:
     1) the initial point,
     2) the distance and
     3) the back azimuth
    """
    # Radius of the Earth
    R = 6371
    # Convert degrees to radians
    baz_ = np.radians(baz)
    # Current lat point converted to radians
    lat1 = np.radians(lat1)
    # Current long point converted to radians
    lon1 = np.radians(lon1)

    lat2 = np.arcsin(np.sin(lat1) * np.cos(d / R) + np.cos(lat1) * np.sin(d / R) * np.cos(baz_))
    lon2 = lon1 + np.arctan2(np.sin(baz_) * np.sin(d / R) * np.cos(lat1),
                             np.cos(d / R) - np.sin(lat1) * np.sin(lat2))
    # Convert to degrees
    lat2 = np.degrees(lat2)
    lon2 = np.degrees(lon2)
    return lat2, lon2


def az2baz(angle):
    """"""
    if angle >= 180.:
        baz = angle - 180.
    else:
        baz = angle + 180.
    return baz


def get_perpendicular_azimuth(az):
    if az > 270:
        perp_az1 = (az + 90) - 360
        perp_az2 = az - 90
    elif az < 90:
        perp_az1 = az + 90
        perp_az2 = az - 90 + 360
    else:
        perp_az1 = az + 90
        perp_az2 = az - 90
    return perp_az1, perp_az2


def project(station_lats, station_lons, point_lat, point_lon, angle):
    """
    Projects stations coordinates to a given point (lon, lat) in respect to an angle to the north.

    NOTE: Takes station coordinates and projects them with respect to the
          center of the profile and the angle of the profile with respect to the North direction.
          Output is in [km] for x,y coordinates with respect to lono and lato

    :param station_lats: Seismic station's latitudes in degrees.
    :param station_lons: Seismic station's longitudes in degrees.
    :param point_lat: Given point's latitude in degrees.
    :param point_lon: Given point's longitude in degrees.
    :param angle: Azimuthal angle in degrees.

    :returns: Distance in km parallel and perpendicular to the given line.
    """

    ylat = (station_lats - point_lat) * 111.19
    xlon = (station_lons - point_lon) * 111.19 * np.cos(np.radians(station_lats))

    M = np.array([[np.cos(np.radians(angle)), np.sin(np.radians(angle))],
                  [-np.sin(np.radians(angle)), np.cos(np.radians(angle))],])
    R = np.dot(np.column_stack((xlon, ylat)), M)

    distx = R[:, 1]
    disty = R[:, 0]

    return distx, disty


def project_stations(sta, ori_prof, point_lat, point_lon):
    """
    Projects stations to a given line.

    :type sta: Pandas DataFrames.
    :param sta: Station details.
    :param ori_prof: Azimuthal angle in degrees.
    :param point_lat: Given point's latitude in degrees.
    :param point_lon: Given point's longitude in degrees.

    :return: Pandas DataFrame with station details and distance along profile and elevation
    """

    xsta, ysta = project(sta["LATSTA"].values, sta["LONSTA"].values, point_lat, point_lon, ori_prof)

    dx, dy = 0, 0
    sta["XSTA"] = xsta + dx
    sta["YSTA"] = ysta + dy
    # Set elevation with negative numbers in km
    sta["ZSTA"] = (-1) * sta["ALTSTA"].values / 1000

    return sta, dx, dy


def plot_moho_picker1(Gp, xx, zz, migration_param_dict, sta, work_directory, profile, profile_name, path4file):
    """

    :param Gp:
    :param xx:
    :param zz:
    :param migration_param_dict:
    :param work_directory:
    :return:
    """
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    import matplotlib.gridspec as gridspec
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    # for picking moho deps
    import pandas as pd
    from matplotlib.colors import LinearSegmentedColormap
    import matplotlib
    matplotlib.use('TkAgg')
    import numpy as np
    import matplotlib.pyplot as plt
    from obspy.geodetics import degrees2kilometers, kilometers2degrees
    from math import radians, degrees, sin, cos, asin, acos, sqrt

    def great_circle(lon1, lat1, lon2, lat2):
        lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
        return 6371 * (acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2)))


    fontsize = 12
    markersize = 100

    XX, ZZ = np.meshgrid(xx, zz)

    pal_col = work_directory + "/data/colormaps/vik.txt"
    pal_col = pd.read_csv(pal_col, header=None, index_col=False, sep="\s+", names=["R", "G", "B"])
    cm = LinearSegmentedColormap.from_list("blue2red", pal_col.values, len(pal_col))
    c = np.min([np.max(Gp), 0.1])
    c = 0.15
    CL = 2

    plt.close('all')
    # PLOT
    f = plt.figure(1, figsize=[15,5])
    gs0 = gridspec.GridSpec(nrows=1, ncols=1, figure=f,
                            hspace=0.08, right=0.91, left=0.09, bottom=0.08, top=0.96,)

    ax = f.add_subplot(gs0[0])  # Ray tracing
    # bwr, seismic, coolwarm, RdBu
    m = ax.pcolormesh(XX, ZZ, Gp.T, cmap=cm, vmin=-c / CL, vmax=c / CL,
                      zorder=1, shading="auto")
    add_colorbar(ax, m)
    #ax.scatter(sta["XSTA"].values, sta["ZSTA"].values,
     #          markersize, facecolors="grey", edgecolors="k",
     #          marker="v", lw=0.95, zorder=3, clip_on=False,
     #          label="Seismic stations", )

    ax.set_aspect("equal")
    ax.set_xlabel("x [km]", fontsize=fontsize)
    ax.set_ylabel("z [km]", fontsize=fontsize)

    majorLocator = MultipleLocator(50)
    minorLocator = MultipleLocator(10)
    ax.xaxis.set_major_locator(majorLocator)
    ax.xaxis.set_minor_locator(minorLocator)
    ax.set_xticks(np.arange(0, 3000, 100))
    ax.set_xlim([xx[0] - 10 , xx[-1] + 40])

    majorLocator = MultipleLocator(50)
    minorLocator = MultipleLocator(5)
    ax.yaxis.set_major_locator(majorLocator)
    ax.yaxis.set_minor_locator(minorLocator)
    ax.set_yticks(np.arange(50, zz[-5], 50))
    ax.set_ylim([800, 0])

    ax.tick_params(axis="both", which="major", labelsize=fontsize)
    ax.tick_params(axis="both", which="minor", labelsize=fontsize)

    print("|-----------------------------------------------|\n"
          "|              MTZ picker manual               |\n"
          "|-----------------------------------------------|\n"
          "| Make your picks using mouse left button and:  |\n"
          "| - the button c for ~ 410f depth,      |\n"
          "| - the button m for ~ 410a depth.         |\n"
          "| - the button o for ~ 660f depth.         |\n"
          "| - the button b for ~ 660a depth.         |\n"
          "|-----------------------------------------------|")
    # Is it a N-S or a E-W cross section?
    if profile[0][1] == profile[1][1]:
        orientation = 'W-E'
    elif profile[0][0] == profile[1][0]:
        orientation = 'S-N'
    else:
        assert False, (
                'OH NO! Cross-section is neither E-W or S-W!'
                'MTZ picker only supports East to West or South to North orientations for the time being...')

    def onkey(event):
        # print(event.key)
        if event.key == 'c':
            if event.xdata is not None and event.ydata is not None:
                print('Dist:', event.xdata, '410f:', event.ydata)
                if orientation == 'S-N':
                    lat = profile[0][1] + kilometers2degrees(event.xdata)
                    print('Lon: ', profile[0][0], 'Lat: ', lat, '410f:', event.ydata)
                    lon = profile[0][0]
                else:
                    lon = profile[0][0] + kilometers2degrees(event.xdata)
                    print('Lon: ', lon, 'Lat: ', profile[0][1], '410f:', event.ydata)
                    lat = profile[0][1]
                # write 410f depths
                # this line is for the cross-sections great circle distance
                gc_dist = kilometers2degrees(great_circle(profile[0][0],profile[0][1],lon,lat))
                with open(path4file + '/410f_depths_' + profile_name + '.txt', 'a') as of:
                    of.write('{}, {}, {}, {}\n'.
                             format(lon, lat, event.ydata, gc_dist))
                ax.plot(event.xdata, event.ydata, markeredgecolor='cyan', marker='D',
                        markerfacecolor='cyan', linestyle='',markersize=7,
                        linewidth=2, alpha=1, label='410f depth')
                f.canvas.draw()
        elif event.key == 'm':
            if event.xdata is not None and event.ydata is not None:
                print('Dist:', event.xdata, '410a:', event.ydata)
                if orientation == 'S-N':
                    lat = profile[0][1] + kilometers2degrees(event.xdata)
                    print('Lon: ', profile[0][0], 'Lat: ', lat, '410a:', event.ydata)
                    lon = profile[0][0]
                else:
                    lon = profile[0][0] + kilometers2degrees(event.xdata)
                    print('Lon: ', lon, 'Lat: ', profile[0][1], '410a:', event.ydata)
                    lat = profile[0][1]
                # Write 410a depths
                # this line is for the cross-sections great circle distance
                gc_dist = kilometers2degrees(great_circle(profile[0][0],profile[0][1],lon,lat))
                with open(path4file + '/410a_depths_' + profile_name + '.txt', 'a') as of:
                    of.write('{}, {}, {}, {}\n'.
                             format(lon, lat, event.ydata, gc_dist))
                ax.plot(event.xdata, event.ydata, markeredgecolor='magenta', marker='D',
                        markerfacecolor='magenta', linestyle='',markersize=7,
                        linewidth=2, alpha=1, label='410a depth')
                f.canvas.draw()
        elif event.key == 'o':
            if event.xdata is not None and event.ydata is not None:
                print('Dist:', event.xdata, '660f:', event.ydata)
                if orientation == 'S-N':
                    lat = profile[0][1] + kilometers2degrees(event.xdata)
                    print('Lon: ', profile[0][0], 'Lat: ', lat, '660f:', event.ydata)
                    lon = profile[0][0]
                else:
                    lon = profile[0][0] + kilometers2degrees(event.xdata)
                    print('Lon: ', lon, 'Lat: ', profile[0][1], '660f:', event.ydata)
                    lat = profile[0][1]
                # Write 660f depths
                # this line is for the cross-sections great circle distance
                gc_dist = kilometers2degrees(great_circle(profile[0][0],profile[0][1],lon,lat))
                with open(path4file + '/660f_depths_' + profile_name + '.txt', 'a') as of:
                    of.write('{}, {}, {}, {}\n'.
                             format(lon, lat, event.ydata, gc_dist))
                ax.plot(event.xdata, event.ydata, markeredgecolor='orange', marker='D',
                        markerfacecolor='orange', linestyle='',markersize=7,
                        linewidth=2, alpha=1, label='660f depth')
                f.canvas.draw()
        elif event.key == 'b':
            if event.xdata is not None and event.ydata is not None:
                print('Dist:', event.xdata, '660a:', event.ydata)
                if orientation == 'S-N':
                    lat = profile[0][1] + kilometers2degrees(event.xdata)
                    print('Lon: ', profile[0][0], 'Lat: ', lat, '660a:', event.ydata)
                    lon = profile[0][0]
                else:
                    lon = profile[0][0] + kilometers2degrees(event.xdata)
                    print('Lon: ', lon, 'Lat: ', profile[0][1], '660a:', event.ydata)
                    lat = profile[0][1]
                # Write 660a depths
                # this line is for the cross-sections great circle distance
                gc_dist = kilometers2degrees(great_circle(profile[0][0],profile[0][1],lon,lat))
                with open(path4file + '/660a_depths_' + profile_name + '.txt', 'a') as of:
                    of.write('{}, {}, {}, {}\n'.
                             format(lon, lat, event.ydata, gc_dist))
                ax.plot(event.xdata, event.ydata, markeredgecolor='black', marker='D',
                        markerfacecolor='black', linestyle='',markersize=7,
                        linewidth=2, alpha=1, label='660a depth')

                f.canvas.draw()

    if orientation == 'S-N':
        ax.text(xx[0], 13, 'S', fontsize=16, color='black')
        ax.text(xx[-1] - 20 , 13, 'N', fontsize=16)
    elif orientation == 'W-E':
        ax.text(xx[0], 13, 'W', fontsize=16, color='black')
        ax.text(xx[-1] - 20 , 13, 'E', fontsize=16)
    # ax1.text(0.04, 0.17, 'Swath = $\pm$%skm'%(prof_width) , transform=ax1.transAxes)
    ax.set_title(profile_name, pad=20)

    # Plotting a single point outside the window we are plotting so
    # the markers are plotted in the legend
    ax.plot(-20, 10, label='410f', color='cyan', marker='D',
            markerfacecolor='cyan', linestyle='',
            markersize=7, linewidth=2, alpha=1)
    ax.plot(-20, 10, label='410a', color='magenta', marker='D',
            markerfacecolor='magenta', linestyle='',
            markersize=7, linewidth=2, alpha=1)
    ax.plot(-20, 10, label='660f', color='orange', marker='D',
            markerfacecolor='orange', linestyle='',
            markersize=7, linewidth=2, alpha=1)
    ax.plot(-20, 10, markeredgecolor='black', marker='D',
            markerfacecolor='black', linestyle='', markersize=7,
            linewidth=2, alpha=1, label='660a')
    ax.legend(loc="lower right")


    # f.canvas.mpl_connect('pick_event', onkey)
    cid2 = f.canvas.mpl_connect('key_press_event', onkey)

    plt.show()

    return

# Set up paths
if platform.node().startswith('kalmard-laptop'):
    data_root_dir = '/media/kalmard/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/home/kalmard/Desktop/codes/github'
    mig_file_dir = '/home/kalmard/Post_doc/Prf_mtz/Python_3Dmigration'
    hard_drive_dir = '/media/kalmard/SEISMIC_DATA/'
else:
    data_root_dir = '/media/kalmard/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/github'
    mig_file_dir = '/home/kalmard/Post_doc/Prf_mtz/Python_3Dmigration'
    hard_drive_dir = '/media/kalmard/SEISMIC_DATA/'

# Define paths
work_dir = os.getcwd()
# path = work_dir + "/data/RF/RF/"
# path = desktop_dir + "/RF_test/RF/"
# path='/media/kmichailos/SEISMIC_DATA/RF_calculations/RF/'
path = '/home/kalmard/Post_doc/Prf_mtz/matlab_deconvolution/test2/'

## Define MIGRATION parameters
# Ray-tracing parameters
inc = 2  # km
zmax = 800 # km
# Determine study area (x -> perpendicular to the profile)
minx = -13.0 # degrees 2.optional:2 or -4
maxx = 46.0 # degrees 2.optional:30 or 38
pasx = 0.26 # degrees oldest 0.38
miny = 30.0 # degrees 2.optional:41 or 38
maxy = 64.0 # degrees 2.optional:51 or 54
pasy = 0.18 # degrees oldest 0.27
minz = -5 # km
# maxz needs to be >= zmax
maxz = 800 # km
pasz = 2 # km
# Pass all the migration parameters in a dictionary to use them in functions called below
m_params = {'minx': minx, 'maxx': maxx,
            'pasx': pasx, 'pasy': pasy, 'miny': miny, 'maxy': maxy,
            'minz': minz, 'maxz': maxz, 'pasz': pasz, 'inc': inc, 'zmax': zmax}
# read stations
sta = rf_mig.read_stations_from_sac(path2rfs=path)


# COPY PASTE FROM HERE |
#                     /|\
with open('/home/kalmard/Post_doc/Prf_mtz/Python_3Dmigration/All_zhu_model_800.npy', 'rb') as f:
    mObs_ep = np.load(f)

# with open('/home/kmichailos/Desktop/All_iasp91.npy', 'rb') as f:
#     mObs_ia = np.load(f)

##Paralell cross-sectinos
##EU01_UNC
#profile_A = np.array([[2, 41.5],   [2, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU01'
##EU02_UNC
#profile_A = np.array([[2.5, 41.5],   [2.5, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU02'
##EU03_UNC
#profile_A = np.array([[3, 41.5],   [3, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU03'
##EU04_UNC
#profile_A = np.array([[3.5, 41.5],   [3.5, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU04'
##EU05_UNC
#profile_A = np.array([[4, 41.5],   [4, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU5'
##EU06_UNC
#profile_A = np.array([[4.5, 41.5],   [4.5, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU6'
##EU07_UNC
#profile_A = np.array([[5, 41.5],   [5, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU7'
##EU08_UNC
#profile_A = np.array([[5.5, 41.5],   [5.5, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU8'
##EU09_UNC
#profile_A = np.array([[6, 41.5],   [6, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU9'
##EU10_UNC
#profile_A = np.array([[6.5, 41.5],   [6.5, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU10'
##EU11_UNC
#profile_A = np.array([[7, 41.5],   [7, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU11'
##EU12_UNC
#profile_A = np.array([[7.5, 41.5],   [7.5, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU12'
##EU13_UNC
#profile_A = np.array([[8, 41.5],   [8, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU13'
##EU14_UNC
#profile_A = np.array([[8.5, 41.5],   [8.5, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU14'
##EU15_UNC
#profile_A = np.array([[9, 41.5],   [9, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU15'
##EU16_UNC
#profile_A = np.array([[9.5, 41.5],   [9.5, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU16'
##EU17_UNC
#profile_A = np.array([[10, 41.5],   [10, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU17'
##EU18_UNC
#profile_A = np.array([[10.5, 41.5],   [10.5, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU18'
##EU19_UNC
#profile_A = np.array([[11, 41.5],   [11, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU19'
##EU20_UNC
#profile_A = np.array([[11.5, 41.5],   [11.5, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU20'
##EU21_UNC
#profile_A = np.array([[12, 41.5],   [12, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU21'
##EU22_UNC
#profile_A = np.array([[12.5, 41.5],   [12.5, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU22'
##EU23_UNC
#profile_A = np.array([[13, 41.5],   [13, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU23'
##EU24_UNC
#profile_A = np.array([[13.5, 41.5],   [13.5, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU24'
##EU25_UNC
#profile_A = np.array([[14, 41.5],   [14, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU25'
##EU26_UNC
#profile_A = np.array([[14.5, 41.5],   [14.5, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU26'
##EU27_UNC
#profile_A = np.array([[15, 41.5],   [15, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU27'
##EU28_UNC
#profile_A = np.array([[15.5, 41.5],   [15.5, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU28'
##EU29_UNC
#profile_A = np.array([[16, 41.5],   [16, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU29'
##EU30_UNC
#profile_A = np.array([[16.5, 41.5],   [16.5, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU30'
##EU31_UNC
#profile_A = np.array([[17, 41.5],   [17, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU31'
##EU32_UNC
#profile_A = np.array([[17.5, 41.5],   [17.5, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU32'
##EU33_UNC
#profile_A = np.array([[18, 41.5],   [18, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU33'
##EU34_UNC
#profile_A = np.array([[18.5, 41.5],   [18.5, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU34'
##EU35_UNC
#profile_A = np.array([[19, 41.5],   [19, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU35'
##EU36_UNC
#profile_A = np.array([[19.5, 41.5],   [19.5, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU36'
##EU37_UNC
#profile_A = np.array([[20, 41.5],   [20, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU37'
##EU38_UNC
#profile_A = np.array([[20.5, 41.5],   [20.5, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU38'
##EU39_UNC
#profile_A = np.array([[21, 41.5],   [21, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU39'
##EU40_UNC
#profile_A = np.array([[21.5, 41.5],   [21.5, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU40'
##EU41_UNC
#profile_A = np.array([[22, 41.5],   [22, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU41'
##EU42_UNC
#profile_A = np.array([[22.5, 41.5],   [22.5, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU42'
##EU43_UNC
#profile_A = np.array([[23, 41.5],   [23, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU43'
##EU44_UNC
#profile_A = np.array([[23.5, 41.5],   [23.5, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU44'
##EU45_UNC
#profile_A = np.array([[24, 41.5],   [24, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU45'
##EU46_UNC
#profile_A = np.array([[24.5, 41.5],   [24.5, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU46'
##EU47_UNC
#profile_A = np.array([[25, 41.5],   [25, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU47'
##EU48_UNC
#profile_A = np.array([[25.5, 41.5],   [25.5, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU48'
##EU49_UNC
#profile_A = np.array([[26, 41.5],   [26, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU49'
##EU50_UNC
#profile_A = np.array([[26.5, 41.5],   [26.5, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU50'
###EU51_UNC
#profile_A = np.array([[27, 41.5],   [27, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU51'
##EU52_UNC
#profile_A = np.array([[27.5, 41.5],   [27.5, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU52'
##EU53_UNC
#profile_A = np.array([[28, 41.5],   [28, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU53'
##EU54_UNC
#profile_A = np.array([[28.5, 41.5],   [28.5, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU54'
##EU55_UNC
#profile_A = np.array([[29, 41.5],   [29, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU55'
##EU56_UNC
#profile_A = np.array([[29.5, 41.5],   [29.5, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU56'
##EU57_UNC
#profile_A = np.array([[30, 41.5],   [30, 51]]) 
#prof_name = 'Cross_section_Zhu_unc_EU57'
##EU58_UNC
#profile_A = np.array([[2, 51],   [30, 51]])
#prof_name = 'Cross_section_Zhu_unc_EU58'
##EU59_UNC
#profile_A = np.array([[2, 50.5], [30, 50.5]])
#prof_name = 'Cross_section_Zhu_unc_EU59'
##EU60_UNC
#profile_A = np.array([[2, 50], [30, 50]])
#prof_name = 'Cross_section_Zhu_unc_EU60'
##EU61_UNC
#profile_A = np.array([[2, 49.5], [30, 49.5]])
#prof_name = 'Cross_section_Zhu_unc_EU61'
##EU62_UNC
#profile_A = np.array([[2, 49], [30, 49]])
#prof_name = 'Cross_section_Zhu_unc_EU62'
##EU63_UNC
#profile_A = np.array([[2, 48.5], [30, 48.5]])
#prof_name = 'Cross_section_Zhu_unc_EU63'
##EU64_UNC
#profile_A = np.array([[2, 48], [30, 48]])
#prof_name = 'Cross_section_Zhu_unc_EU64'
##EU65_UNC
#profile_A = np.array([[2, 47.5], [30, 47.5]])
#prof_name = 'Cross_section_Zhu_unc_EU65'
##EU66_UNC
#profile_A = np.array([[2, 47], [30, 47]])
#prof_name = 'Cross_section_Zhu_unc_EU66'
###EU67_UNC
#profile_A = np.array([[2, 46.5], [30, 46.5]])
#prof_name = 'Cross_section_Zhu_unc_EU67'
##EU68_UNC
#profile_A = np.array([[2, 46], [30, 46]])
#prof_name = 'Cross_section_Zhu_unc_EU68'
##EU69_UNC
#profile_A = np.array([[2, 45.5], [30, 45.5]])
#prof_name = 'Cross_section_Zhu_unc_EU69'
##EU70_UNC
#profile_A = np.array([[2, 45], [30, 45]])
#prof_name = 'Cross_section_Zhu_unc_EU70'
##EU71_UNC
#profile_A = np.array([[2, 44.5], [30, 44.5]])
#prof_name = 'Cross_section_Zhu_unc_EU71'
##EU72_UNC
#profile_A = np.array([[2, 44], [30, 44]])
#prof_name = 'Cross_section_Zhu_unc_EU72'
##EU73_UNC
#profile_A = np.array([[2, 43.5], [30, 43.5]])
#prof_name = 'Cross_section_Zhu_unc_EU73'
##EU74_UNC
#profile_A = np.array([[2, 43], [30, 43]])
#prof_name = 'Cross_section_Zhu_unc_EU74'
##EU75_UNC
#profile_A = np.array([[2, 42.5], [30, 42.5]])
#prof_name = 'Cross_section_Zhu_unc_EU75'
##EU76_UNC
#profile_A = np.array([[2, 42], [30, 42]])
#prof_name = 'Cross_section_Zhu_unc_EU76'
##EU77_UNC
profile_A = np.array([[2, 41.5], [30, 41.5]])
prof_name = 'Cross_section_Zhu_unc_EU77'

# swath=37.5
G2_, sta_, xx, zz = plot_migration_sphr.create_2d_profile_4_moho_picker(mObs_ep, m_params, profile_A, sta, swath=300, plot=True)
G2 = rf_mig.ccp_smooth(G2_, m_params)
# G2[np.abs(G2) < np.max(np.abs(G2)) * 15 / 100] = 0
G2 = rf_mig.ccpFilter(G2)
# Manually pick moho deps
# IMPORTANT NOTE: only works with cross-sections the have S-N and W-E directions!!!
plot_moho_picker1(Gp=G2, xx=xx, zz=zz, migration_param_dict=m_params,
                                sta=sta_, work_directory=work_dir, profile=profile_A, profile_name=prof_name,
                                path4file= work_dir )#)+ '/rfmpy/visualisation/gmt/maps/files/moho_picks/')
