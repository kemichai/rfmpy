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
          "| - the button m for ~ 410 depth,      |\n"
          "| - the button u for ~ 660 depth.         |\n"
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
        if event.key == 'm':
            if event.xdata is not None and event.ydata is not None:
                print('Dist:', event.xdata, '410:', event.ydata)
                if orientation == 'S-N':
                    lat = profile[0][1] + kilometers2degrees(event.xdata)
                    print('Lon: ', profile[0][0], 'Lat: ', lat, '410:', event.ydata)
                    lon = profile[0][0]
                else:
                    lon = profile[0][0] + kilometers2degrees(event.xdata)
                    print('Lon: ', lon, 'Lat: ', profile[0][1], '410:', event.ydata)
                    lat = profile[0][1]
                # write 410 depths
                # this line is for the cross-sections great circle distance
                gc_dist = kilometers2degrees(great_circle(profile[0][0],profile[0][1],lon,lat))
                with open(path4file + '/410_depths_' + profile_name + '.txt', 'a') as of:
                    of.write('{}, {}, {}, {}\n'.
                             format(lon, lat, event.ydata, gc_dist))
                ax.plot(event.xdata, event.ydata, label='~410 depth',
                        color='black', marker='D',markerfacecolor='white',linestyle='',
                        markersize=7, linewidth=2, alpha=1)
                f.canvas.draw()
        elif event.key == 'u':
            if event.xdata is not None and event.ydata is not None:
                print('Dist:', event.xdata, '660:', event.ydata)
                if orientation == 'S-N':
                    lat = profile[0][1] + kilometers2degrees(event.xdata)
                    print('Lon: ', profile[0][0], 'Lat: ', lat, '660:', event.ydata)
                    lon = profile[0][0]
                else:
                    lon = profile[0][0] + kilometers2degrees(event.xdata)
                    print('Lon: ', lon, 'Lat: ', profile[0][1], '660:', event.ydata)
                    lat = profile[0][1]
                # Write 660 depths
                # this line is for the cross-sections great circle distance
                gc_dist = kilometers2degrees(great_circle(profile[0][0],profile[0][1],lon,lat))
                with open(path4file + '/660_depths_' + profile_name + '.txt', 'a') as of:
                    of.write('{}, {}, {}, {}\n'.
                             format(lon, lat, event.ydata, gc_dist))
                ax.plot(event.xdata, event.ydata, markeredgecolor='black', marker='D',
                        markerfacecolor='gray', linestyle='',markersize=7,
                        linewidth=2, alpha=1, label='660')

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
    ax.plot(-20, 10, label='410', color='black', marker='D',
            markerfacecolor='white', linestyle='',
            markersize=7, linewidth=2, alpha=1)
    ax.plot(-20, 10, markeredgecolor='black', marker='D',
            markerfacecolor='gray', linestyle='', markersize=7,
            linewidth=2, alpha=1, label='660')
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
#EU01
#profile_A = np.array([[2, 41.5],   [2, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU01'
##EU02
#profile_A = np.array([[2.2, 41.5],   [2.2, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU02'
##EU03
#profile_A = np.array([[2.4, 41.5],   [2.4, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU03'
##EU04
#profile_A = np.array([[2.6, 41.5],   [2.6, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU04'
##EU05
#profile_A = np.array([[2.8, 41.5],   [2.8, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU05'
##EU06
#profile_A = np.array([[3, 41.5],   [3, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU06'
##EU07
#profile_A = np.array([[3.2, 41.5],   [3.2, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU07'
##EU08
#profile_A = np.array([[3.4, 41.5],   [3.4, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU08'
##EU09
#profile_A = np.array([[3.6, 41.5],   [3.6, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU09'
##EU10
#profile_A = np.array([[3.8, 41.5],   [3.8, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU10'
##EU11
#profile_A = np.array([[4, 41.5],   [4, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU11'
##EU12
#profile_A = np.array([[4.2, 41.5],   [4.2, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU12'
##EU13
#profile_A = np.array([[4.4, 41.5],   [4.4, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU13'
##EU14
#profile_A = np.array([[4.6, 41.5],   [4.6, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU14'
##EU15
#profile_A = np.array([[4.8, 41.5],   [4.8, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU15'
##EU16
#profile_A = np.array([[5, 41.5],   [5, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU16'
##EU17
#profile_A = np.array([[5.2, 41.5],   [5.2, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU17'
##EU18
#profile_A = np.array([[5.4, 41.5],   [5.4, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU18'
##EU19
#profile_A = np.array([[5.6, 41.5],   [5.6, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU19'
##EU20
#profile_A = np.array([[5.8, 41.5],   [5.8, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU20'
##EU21
#profile_A = np.array([[6, 41.5],   [6, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU21'
##EU22
#profile_A = np.array([[6.2, 41.5],   [6.2, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU22'
##EU23
#profile_A = np.array([[6.4, 41.5],   [6.4, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU23'
##EU24
#profile_A = np.array([[6.6, 41.5],   [6.6, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU24'
##EU25
#profile_A = np.array([[6.8, 41.5],   [6.8, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU25'
##EU26
#profile_A = np.array([[7, 41.5],   [7, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU26'
##EU27
#profile_A = np.array([[7.2, 41.5],   [7.2, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU27'
##EU28
#profile_A = np.array([[7.4, 41.5],   [7.4, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU28'
##EU29
#profile_A = np.array([[7.6, 41.5],   [7.6, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU29'
#EU30
#profile_A = np.array([[7.8, 41.5],   [7.8, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU30'
##EU31
#profile_A = np.array([[8, 41.5],   [8, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU31'
##EU32
#profile_A = np.array([[8.2, 41.5],   [8.2, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU32'
##EU33
#profile_A = np.array([[8.4, 41.5],   [8.4, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU33'
##EU34
#profile_A = np.array([[8.6, 41.5],   [8.6, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU34'
##EU35
#profile_A = np.array([[8.8, 41.5],   [8.8, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU35'
##EU36
#profile_A = np.array([[9, 41.5],   [9, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU36'
##EU37
#profile_A = np.array([[9.2, 41.5],   [9.2, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU37'
##EU38
#profile_A = np.array([[9.4, 41.5],   [9.4, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU38'
##EU39
#profile_A = np.array([[9.6, 41.5],   [9.6, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU39'
##EU40
#profile_A = np.array([[9.8, 41.5],   [9.8, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU40'
##EU41
#profile_A = np.array([[10, 41.5],   [10, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU41'
##EU42
#profile_A = np.array([[10.2, 41.5],   [10.2, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU42'
##EU43
#profile_A = np.array([[10.4, 41.5],   [10.4, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU43'
##EU44
#profile_A = np.array([[10.6, 41.5],   [10.6, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU44'
##EU45
#profile_A = np.array([[10.8, 41.5],   [10.8, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU45'
##EU46
#profile_A = np.array([[11, 41.5],   [11, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU46'
##EU47
#profile_A = np.array([[11.2, 41.5],   [11.2, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU47'
##EU48
#profile_A = np.array([[11.4, 41.5],   [11.4, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU48'
##EU49
#profile_A = np.array([[11.6, 41.5],   [11.6, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU49'
##EU50
#profile_A = np.array([[11.8, 41.5],   [11.8, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU50'
##EU51
#profile_A = np.array([[12, 41.5],   [12, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU51'
##EU52
#profile_A = np.array([[12.2, 41.5],   [12.2, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU52'
##EU53
#profile_A = np.array([[12.4, 41.5],   [12.4, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU53'
##EU54
#profile_A = np.array([[12.6, 41.5],   [12.6, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU54'
##EU55
#profile_A = np.array([[12.8, 41.5],   [12.8, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU55'
##EU56
#profile_A = np.array([[13, 41.5],   [13, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU56'
##EU57
#profile_A = np.array([[13.2, 41.5],   [13.2, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU57'
##EU58
#profile_A = np.array([[13.4, 41.5],   [13.4, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU58'
#EU59
#profile_A = np.array([[13.6, 41.5],   [13.6, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU59'
##EU60
#profile_A = np.array([[13.8, 41.5],   [13.8, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU60'
##EU61
#profile_A = np.array([[14, 41.5],   [14, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU61'
##EU62
#profile_A = np.array([[14.2, 41.5],   [14.2, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU62'
##EU63
#profile_A = np.array([[14.4, 41.5],   [14.4, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU63'
##EU64
#profile_A = np.array([[14.6, 41.5],   [14.6, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU64'
##EU65
#profile_A = np.array([[14.8, 41.5],   [14.8, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU65'
##EU66
#profile_A = np.array([[15, 41.5],   [15, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU66'
##EU67
#profile_A = np.array([[15.2, 41.5],   [15.2, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU67'
##EU68
#profile_A = np.array([[15.4, 41.5],   [15.4, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU68'
##EU69
#profile_A = np.array([[15.6, 41.5],   [15.6, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU69'
##EU70
#profile_A = np.array([[15.8, 41.5],   [15.8, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU70'
##EU71
#profile_A = np.array([[16, 41.5],   [16, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU71'
##EU72
#profile_A = np.array([[16.2, 41.5],   [16.2, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU72'
##EU73
#profile_A = np.array([[16.4, 41.5],   [16.4, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU73'
##EU74
#profile_A = np.array([[16.6, 41.5],   [16.6, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU74'
##EU75
#profile_A = np.array([[16.8, 41.5],   [16.8, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU75'
##EU76
#profile_A = np.array([[17, 41.5],   [17, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU76'
##EU77
#profile_A = np.array([[17.2, 41.5],   [17.2, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU77'
##EU78
#profile_A = np.array([[17.4, 41.5],   [17.4, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU78'
##EU79
#profile_A = np.array([[17.6, 41.5],   [17.6, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU79'
##EU80
#profile_A = np.array([[17.8, 41.5],   [17.8, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU80'
##EU81
#profile_A = np.array([[18, 41.5],   [18, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU81'
##EU82
#profile_A = np.array([[18.2, 41.5],   [18.2, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU82'
##EU83
#profile_A = np.array([[18.4, 41.5],   [18.4, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU83'
##EU84
#profile_A = np.array([[18.6, 41.5],   [18.6, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU84'
##EU85
#profile_A = np.array([[18.8, 41.5],   [18.8, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU85'
##EU86
#profile_A = np.array([[19, 41.5],   [19, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU86'
##EU87
#profile_A = np.array([[19.2, 41.5],   [19.2, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU87'
#EU88
#profile_A = np.array([[19.4, 41.5],   [19.4, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU88'
##EU89
#profile_A = np.array([[19.6, 41.5],   [19.6, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU89'
##EU90
#profile_A = np.array([[19.8, 41.5],   [19.8, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU90'
##EU91
#profile_A = np.array([[20, 41.5],   [20, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU91'
##EU92
#profile_A = np.array([[20.2, 41.5],   [20.2, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU92'
##EU93
#profile_A = np.array([[20.4, 41.5],   [20.4, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU93'
##EU94
#profile_A = np.array([[20.6, 41.5],   [20.6, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU94'
##EU95
#profile_A = np.array([[20.8, 41.5],   [20.8, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU95'
##EU96
#profile_A = np.array([[21, 41.5],   [21, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU96'
##EU97
#profile_A = np.array([[21.2, 41.5],   [21.2, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU97'
##EU98
#profile_A = np.array([[21.4, 41.5],   [21.4, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU98'
##EU99
#profile_A = np.array([[21.6, 41.5],   [21.6, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU99'
##EU100
#profile_A = np.array([[21.8, 41.5],   [21.8, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU100'
##EU101
#profile_A = np.array([[22, 41.5],   [22, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU101'
##EU102
#profile_A = np.array([[22.2, 41.5],   [22.2, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU102'
##EU103
#profile_A = np.array([[22.4, 41.5],   [22.4, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU103'
##EU104
#profile_A = np.array([[22.6, 41.5],   [22.6, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU104'
##EU105
#profile_A = np.array([[22.8, 41.5],   [22.8, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU105'
##EU106
#profile_A = np.array([[23, 41.5],   [23, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU106'
##EU107
#profile_A = np.array([[23.2, 41.5],   [23.2, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU107'
##EU108
#profile_A = np.array([[23.4, 41.5],   [23.4, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU108'
##EU109
#profile_A = np.array([[23.6, 41.5],   [23.6, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU109'
##EU110
#profile_A = np.array([[23.8, 41.5],   [23.8, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU110'
##EU111
#profile_A = np.array([[24, 41.5],   [24, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU111'
##EU112
#profile_A = np.array([[24.2, 41.5],   [24.2, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU112'
##EU113
#profile_A = np.array([[24.4, 41.5],   [24.4, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU113'
##EU114
#profile_A = np.array([[24.6, 41.5],   [24.6, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU114'
##EU115
#profile_A = np.array([[24.8, 41.5],   [24.8, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU115'
##EU116
#profile_A = np.array([[25, 41.5],   [25, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU116'
#EU117
#profile_A = np.array([[25.2, 41.5],   [25.2, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU117'
##EU118
#profile_A = np.array([[25.4, 41.5],   [25.4, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU118'
##EU119
#profile_A = np.array([[25.6, 41.5],   [25.6, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU119'
##EU120
#profile_A = np.array([[25.8, 41.5],   [25.8, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU120'
##EU121
#profile_A = np.array([[26, 41.5],   [26, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU121'
##EU122
#profile_A = np.array([[26.2, 41.5],   [26.2, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU122'
##EU123
#profile_A = np.array([[26.4, 41.5],   [26.4, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU123'
##EU124
#profile_A = np.array([[26.6, 41.5],   [26.6, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU124'
##EU125
#profile_A = np.array([[26.8, 41.5],   [26.8, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU125'
##EU126
#profile_A = np.array([[27, 41.5],   [27, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU126'
##EU127
#profile_A = np.array([[27.2, 41.5],   [27.2, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU127'
##EU128
#profile_A = np.array([[27.4, 41.5],   [27.4, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU128'
##EU129
#profile_A = np.array([[27.6, 41.5],   [27.6, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU129'
##EU130
#profile_A = np.array([[27.8, 41.5],   [27.8, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU130'
##EU131
#profile_A = np.array([[28, 41.5],   [28, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU131'
##EU132
#profile_A = np.array([[28.2, 41.5],   [28.2, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU132'
##EU133
#profile_A = np.array([[28.4, 41.5],   [28.4, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU133'
##EU134
#profile_A = np.array([[28.6, 41.5],   [28.6, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU134'
#EU135
#profile_A = np.array([[28.8, 41.5],   [28.8, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU135'
##EU136
#profile_A = np.array([[29, 41.5],   [29, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU136'
##EU137
#profile_A = np.array([[29.2, 41.5],   [29.2, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU137'
##EU138
#profile_A = np.array([[29.4, 41.5],   [29.4, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU138'
##EU139
#profile_A = np.array([[29.6, 41.5],   [29.6, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU139'
##EU140
#profile_A = np.array([[29.8, 41.5],   [29.8, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU140'
##EU141
#profile_A = np.array([[30, 41.5],   [30, 51]]) 
#prof_name = 'Cross_section_Zhu_m_EU141'
##EU142
#profile_A = np.array([[2, 51],   [30, 51]])
#prof_name = 'Cross_section_Zhu_m_EU142'
##EU143
#profile_A = np.array([[2, 50.8], [30, 50.8]])
#prof_name = 'Cross_section_Zhu_m_EU143'
##EU144
#profile_A = np.array([[2, 50.6], [30, 50.6]])
#prof_name = 'Cross_section_Zhu_m_EU144'
##EU145
#profile_A = np.array([[2, 50.4], [30, 50.4]])
#prof_name = 'Cross_section_Zhu_m_EU145'
##EU146
#profile_A = np.array([[2, 50.2], [30, 50.2]])
#prof_name = 'Cross_section_Zhu_m_EU146'
##EU147
#profile_A = np.array([[2, 50], [30, 50]])
#prof_name = 'Cross_section_Zhu_m_EU147'
##EU148
#profile_A = np.array([[2, 49.8], [30, 49.8]])
#prof_name = 'Cross_section_Zhu_m_EU148'
##EU149
#profile_A = np.array([[2, 49.6], [30, 49.6]])
#prof_name = 'Cross_section_Zhu_m_EU149'
##EU150
#profile_A = np.array([[2, 49.4], [30, 49.4]])
#prof_name = 'Cross_section_Zhu_m_EU150'
##EU151
#profile_A = np.array([[2, 49.2], [30, 49.2]])
#prof_name = 'Cross_section_Zhu_m_EU151'
##EU152
#profile_A = np.array([[2, 49], [30, 49]])
#prof_name = 'Cross_section_Zhu_m_EU152'
##EU153
#profile_A = np.array([[2, 48.8], [30, 48.8]])
#prof_name = 'Cross_section_Zhu_m_EU153'
##EU154
#profile_A = np.array([[2, 48.6], [30, 48.6]])
#prof_name = 'Cross_section_Zhu_m_EU154'
##EU155
#profile_A = np.array([[2, 48.4], [30, 48.4]])
#prof_name = 'Cross_section_Zhu_m_EU155'
##EU156
#profile_A = np.array([[2, 48.2], [30, 48.2]])
#prof_name = 'Cross_section_Zhu_m_EU156'
##EU157
#profile_A = np.array([[2, 48], [30, 48]])
#prof_name = 'Cross_section_Zhu_m_EU157'
##EU158
#profile_A = np.array([[2, 47.8], [30, 47.8]])
#prof_name = 'Cross_section_Zhu_m_EU158'
##EU159
#profile_A = np.array([[2, 47.6], [30, 47.6]])
#prof_name = 'Cross_section_Zhu_m_EU159'
##EU160
#profile_A = np.array([[2, 47.4], [30, 47.4]])
#prof_name = 'Cross_section_Zhu_m_EU160'
##EU161
#profile_A = np.array([[2, 47.2], [30, 47.2]])
#prof_name = 'Cross_section_Zhu_m_EU161'
##EU162
#profile_A = np.array([[2, 47], [30, 47]])
#prof_name = 'Cross_section_Zhu_m_EU162'
##EU163
#profile_A = np.array([[2, 46.8], [30, 46.8]])
#prof_name = 'Cross_section_Zhu_m_EU163'
##EU164
#profile_A = np.array([[2, 46.6], [30, 46.6]])
#prof_name = 'Cross_section_Zhu_m_EU164'
##EU165
#profile_A = np.array([[2, 46.4], [30, 46.4]])
#prof_name = 'Cross_section_Zhu_m_EU165'
##EU166
#profile_A = np.array([[2, 46.2], [30, 46.2]])
#prof_name = 'Cross_section_Zhu_m_EU166'
##EU167
#profile_A = np.array([[2, 46], [30, 46]])
#prof_name = 'Cross_section_Zhu_m_EU167'
##EU168
#profile_A = np.array([[2, 45.8], [30, 45.8]])
#prof_name = 'Cross_section_Zhu_m_EU168'
##EU169
#profile_A = np.array([[2, 45.6], [30, 45.6]])
#prof_name = 'Cross_section_Zhu_m_EU169'
##EU170
#profile_A = np.array([[2, 45.4], [30, 45.4]])
#prof_name = 'Cross_section_Zhu_m_EU170'
##EU171
#profile_A = np.array([[2, 45.2], [30, 45.2]])
#prof_name = 'Cross_section_Zhu_m_EU171'
##EU172
#profile_A = np.array([[2, 45], [30, 45]])
#prof_name = 'Cross_section_Zhu_m_EU172'
##EU173
#profile_A = np.array([[2, 44.8], [30, 44.8]])
#prof_name = 'Cross_section_Zhu_m_EU173'
##EU174
#profile_A = np.array([[2, 44.6], [30, 44.6]])
#prof_name = 'Cross_section_Zhu_m_EU174'
##EU175
#profile_A = np.array([[2, 44.4], [30, 44.4]])
#prof_name = 'Cross_section_Zhu_m_EU175'
##EU176
#profile_A = np.array([[2, 44.2], [30, 44.2]])
#prof_name = 'Cross_section_Zhu_m_EU176'
##EU177
#profile_A = np.array([[2, 44], [30, 44]])
#prof_name = 'Cross_section_Zhu_m_EU177'
##EU178
#profile_A = np.array([[2, 43.8], [30, 43.8]])
#prof_name = 'Cross_section_Zhu_m_EU178'
##EU179
#profile_A = np.array([[2, 43.6], [30, 43.6]])
#prof_name = 'Cross_section_Zhu_m_EU179'
##EU180
#profile_A = np.array([[2, 43.4], [30, 43.4]])
#prof_name = 'Cross_section_Zhu_m_EU180'
##EU181
#profile_A = np.array([[2, 43.2], [30, 43.2]])
#prof_name = 'Cross_section_Zhu_m_EU181'
##EU182
#profile_A = np.array([[2, 43], [30, 43]])
#prof_name = 'Cross_section_Zhu_m_EU182'
##EU183
#profile_A = np.array([[2, 42.8], [30, 42.8]])
#prof_name = 'Cross_section_Zhu_m_EU183'
##EU184
#profile_A = np.array([[2, 42.6], [30, 42.6]])
#prof_name = 'Cross_section_Zhu_m_EU184'
##EU185
#profile_A = np.array([[2, 42.4], [30, 42.4]])
#prof_name = 'Cross_section_Zhu_m_EU185'
##EU186
#profile_A = np.array([[2, 42.2], [30, 42.2]])
#prof_name = 'Cross_section_Zhu_m_EU186'
##EU187
#profile_A = np.array([[2, 42], [30, 42]])
#prof_name = 'Cross_section_Zhu_m_EU187'
##EU188
#profile_A = np.array([[2, 41.8], [30, 41.8]])
#prof_name = 'Cross_section_Zhu_m_EU188'
##EU189
#profile_A = np.array([[2, 41.6], [30, 41.6]])
#prof_name = 'Cross_section_Zhu_m_EU189'
##EU190
profile_A = np.array([[2, 41.4], [30, 41.4]])
prof_name = 'Cross_section_Zhu_m_EU190'

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
