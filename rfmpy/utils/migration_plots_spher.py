"""
Functions for calculating 3D migration of RFs in cartesian coordinates.

Note: Based on codes originally written by Matteo Scarponi.

Location: Chavannes-pres-renens, CH
Date: Mar 2022
Author: Konstantinos Michailos
"""

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


def plot_migration_profile_old(Gp, migration_param_dict, sta, work_directory, filename=False):


    # Read migration parameters
    minx = migration_param_dict['minx']
    maxx = migration_param_dict['maxx']
    pasx = migration_param_dict['pasx']
    miny = migration_param_dict['miny']
    maxy = migration_param_dict['maxy']
    pasy = migration_param_dict['pasy']
    minz = migration_param_dict['minz']
    maxz = migration_param_dict['maxz']
    pasz = migration_param_dict['pasz']
    inc = migration_param_dict['inc']
    zmax = migration_param_dict['zmax']

    # Grid preparation
    xx = np.arange(minx, maxx + pasx, pasx)
    yy = np.arange(miny, maxy + pasy, pasy)
    zz = np.arange(minz, maxz + pasz, pasz)

    XX, ZZ = np.meshgrid(xx, zz)


    # COLOR PALETTE AND COLORMAP
    # cork, bam, broc, vik
    pal_col = work_directory + "/data/colormaps/vik.txt"
    pal_col = pd.read_csv(pal_col, header=None, index_col=False, sep="\s+", names=["R", "G", "B"])
    cm = LinearSegmentedColormap.from_list("blue2red", pal_col.values, len(pal_col))
    c = np.min([np.max(Gp), 0.1])
    c = 0.06
    CL = 2

    # PLOT
    f = plt.figure(1, figsize=[10, 8])
    gs0 = gridspec.GridSpec(nrows=1, ncols=1, figure=f,
                            hspace=0.08, right=0.91, left=0.09, bottom=0.08, top=0.96,)

    ax = f.add_subplot(gs0[0])  # Ray tracing
    # bwr, seismic, coolwarm, RdBu
    m = ax.pcolormesh(XX, ZZ, Gp.T, cmap=cm, vmin=-c / CL, vmax=c / CL,
                      zorder=1, shading="auto")
    add_colorbar(ax, m)
    ax.scatter(sta["LONSTA"].values, sta["ZSTA"].values,
               markersize, facecolors="grey", edgecolors="k",
               marker="v", lw=0.95, zorder=3, clip_on=False,
               label="Seismic stations", )

    ax.set_aspect("equal")
    ax.set_xlabel("x [km]", fontsize=fontsize)
    ax.set_ylabel("z [km]", fontsize=fontsize)

    majorLocator = MultipleLocator(10)
    minorLocator = MultipleLocator(2.5)
    ax.xaxis.set_major_locator(majorLocator)
    ax.xaxis.set_minor_locator(minorLocator)
    # ax.set_xticks(np.arange(0, 140, 10))

    majorLocator = MultipleLocator(10)
    minorLocator = MultipleLocator(2.5)
    ax.yaxis.set_major_locator(majorLocator)
    ax.yaxis.set_minor_locator(minorLocator)
    ax.set_yticks(np.arange(10, 540, 100))
    ax.set_ylim([500, 0])

    ax.tick_params(axis="both", which="major", labelsize=fontsize)
    ax.tick_params(axis="both", which="minor", labelsize=fontsize)
    if filename:
        f.savefig(filename, dpi=200)

    plt.show()
    plt.close()

    return


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


def create_2d_profile(G3, migration_param_dict, profile_points, sta, swath=200, plot=True):
    """

    :param G3:
    :type migration_param_dict: dict
    :param migration_param_dict: Dictionary of grid points for the migration.
    :param profile_points:
    :param sta:
    :param swath: Swath of profile on both sides in km.
    :param plot:
    :return:
    """

    # Read migration parameters
    minx = migration_param_dict['minx']
    maxx = migration_param_dict['maxx']
    pasx = migration_param_dict['pasx']
    miny = migration_param_dict['miny']
    maxy = migration_param_dict['maxy']
    pasy = migration_param_dict['pasy']
    minz = migration_param_dict['minz']
    maxz = migration_param_dict['maxz']
    pasz = migration_param_dict['pasz']

    # Define grid
    grid_3d_x = np.arange(minx, maxx + pasx, pasx)
    grid_3d_y = np.arange(miny, maxy + pasy, pasy)
    grid_3d_z = np.arange(minz, maxz + pasz, pasz)

    # Interpolate in 3D the RF amplitudes
    G_interpolated = RegularGridInterpolator((grid_3d_x, grid_3d_y, grid_3d_z), G3)

    # Profile start and end
    lon0, lat0 = profile_points[0][0], profile_points[0][1]
    lon1, lat1 = profile_points[1][0], profile_points[1][1]

    profile_swath = swath
    # Profile azimuth
    geoid = pyproj.Geod(ellps='WGS84')
    profile_az, back_azimuth, profile_len_ = geoid.inv(lon0, lat0, lon1, lat1)
    # Profile length (km)
    profile_len = profile_len_ / 1000
    # Two perpendicular azimuths
    az1, az2 = get_perpendicular_azimuth(profile_az)

    # By doing this we have roughly the same spacing as the grid
    num_of_points = int(round(profile_len/degrees2kilometers(pasx))) - 1

    # Coordinates of the points along the profile knowing start and end of profile
    # TODO: when I define a finer grid I won't need the * here!!!!!!
    n_extra_points = num_of_points  # number of these points
    print("Number of points along the profile: ", n_extra_points, " Length of profile: ", profile_len)

    geoid = Geod(ellps="WGS84")
    extra_points = np.array(geoid.npts(lon0, lat0, lon1, lat1, n_extra_points))
    # Create new lists of lon, lat, dep and amps (interpolated)
    lon_points_along_prof = extra_points[:, 0]
    lat_points_along_prof = extra_points[:, 1]

    # For each point of the profile find the two end point perpendicular
    # to them in a distance equal to the swath of the profile
    amps = []
    for i, lon in enumerate(lon_points_along_prof):
        # Two points perpendicular to the azimuth of the profile at each point of the profile
        lat_1, lon_1 = get_end_point(lat_points_along_prof[i], lon_points_along_prof[i], az1, profile_swath)
        lat_2, lon_2 = get_end_point(lat_points_along_prof[i], lon_points_along_prof[i], az2, profile_swath)
        # TODO: when I define a finer grid I won't need the * here!!!!!!
        n_extra_points_ =  (int(round(swath/degrees2kilometers(pasx))) - 1 ) # number of these points
        points_perpendicular_2_prof = np.array(geoid.npts(lon_1, lat_1, lon_2, lat_2, n_extra_points_))

        temp_lon = points_perpendicular_2_prof[:, 0]
        temp_lat = points_perpendicular_2_prof[:, 1]
        # interpolate for each of the 5 points perpendicular to the profile
        amps_matrix_temp = np.zeros((len(grid_3d_z)))
        nG = np.zeros((len(grid_3d_z)))

        for j, lon_ in enumerate(temp_lon):
            # print(j)
            amps_temp = np.zeros((len(grid_3d_z)))
            for k, z in enumerate(grid_3d_z):
                # print(z)
                point = np.array([lon_, temp_lat[j], z])
                VPinterp = G_interpolated(point)
                amps_temp[k] = VPinterp[0]
                # print(VPinterp[0])
            amps_matrix_temp = amps_matrix_temp + amps_temp
            nG = nG + 1
        # add and divide by the number of stacks
        G = np.divide(amps_matrix_temp, nG)
        amps.append(G.tolist())
        # Just add (stack)
        # amps.append(amps_matrix_temp.tolist())
    G2 = np.array(amps)  # 2 dimensions
    print("Number of points perpendicular to the profile: ", n_extra_points_, " Swath: ", swath)

    sta, dxSta, dySta = project_stations(sta=sta, ori_prof=profile_az,
                                         point_lat=lat0, point_lon=lon0)
    if plot:
        # Plot stations and profile
        lons = [lon0, lon1]
        lats = [lat0, lat1]
        plt.plot(lons, lats, c='dodgerblue')
        plt.plot(temp_lon, temp_lat, c='gray', linestyle=':', label='Swath (km)')
        plt.scatter(lon0, lat0, c='dodgerblue', marker='s', edgecolor='k', s=50, label='Start')
        plt.scatter(lon1, lat1, c='dodgerblue', marker='o', edgecolor='k', s=50, label='End')
        plt.scatter(sta["LONSTA"], sta["LATSTA"], c='r', marker='v', edgecolor='k', s=100)
        plt.legend()
        plt.ylim(miny, maxy)
        plt.xlim(minx, maxx)
        plt.show()

    # Grid preparation
    xx = np.arange(0, profile_len, profile_len / n_extra_points)
    zz = np.arange(minz, maxz + pasz, pasz)

    return G2, sta, xx, zz


def plot_migration_profile(Gp, xx, zz, migration_param_dict, sta, work_directory, filename=False):
    """

    :param Gp:
    :param xx:
    :param zz:
    :param migration_param_dict:
    :param sta:
    :param work_directory:
    :param filename:
    :return:
    """

    XX, ZZ = np.meshgrid(xx, zz)

    # zz_corr = []
    # R = 6371.009
    # for i, dist in enumerate(prof_dist):
    #     correction = np.sqrt(dist**2 + R)
    #     print(correction)
    #     zz_corr = zz + correction



    # COLOR PALETTE AND COLORMAP
    # cork, bam, broc, vik
    pal_col = work_directory + "/data/colormaps/vik.txt"
    pal_col = pd.read_csv(pal_col, header=None, index_col=False, sep="\s+", names=["R", "G", "B"])
    cm = LinearSegmentedColormap.from_list("blue2red", pal_col.values, len(pal_col))
    c = np.min([np.max(Gp), 0.1])
    c = 0.06
    CL = 2

    # PLOT
    f = plt.figure(1, figsize=[10, 8])
    gs0 = gridspec.GridSpec(nrows=1, ncols=1, figure=f,
                            hspace=0.08, right=0.91, left=0.09, bottom=0.08, top=0.96,)

    ax = f.add_subplot(gs0[0])  # Ray tracing
    # bwr, seismic, coolwarm, RdBu
    m = ax.pcolormesh(XX, ZZ, Gp.T, cmap=cm, vmin=-c / CL, vmax=c / CL,
                      zorder=1, shading="auto")
    add_colorbar(ax, m)
    ax.scatter(sta["XSTA"].values, sta["ZSTA"].values,
               markersize, facecolors="grey", edgecolors="k",
               marker="v", lw=0.95, zorder=3, clip_on=False,
               label="Seismic stations", )

    ax.set_aspect("equal")
    ax.set_xlabel("x [km]", fontsize=fontsize)
    ax.set_ylabel("z [km]", fontsize=fontsize)

    majorLocator = MultipleLocator(10)
    minorLocator = MultipleLocator(2.5)
    ax.xaxis.set_major_locator(majorLocator)
    ax.xaxis.set_minor_locator(minorLocator)
    # ax.set_xticks(np.arange(0, 140, 10))

    majorLocator = MultipleLocator(10)
    minorLocator = MultipleLocator(2.5)
    ax.yaxis.set_major_locator(majorLocator)
    ax.yaxis.set_minor_locator(minorLocator)
    ax.set_yticks(np.arange(10, zz[-1], 10))
    ax.set_ylim([100, 0])

    ax.tick_params(axis="both", which="major", labelsize=fontsize)
    ax.tick_params(axis="both", which="minor", labelsize=fontsize)
    if filename:
        f.savefig(filename, dpi=200)

    plt.show()
    plt.close()

    return



def plot_ray_tracing(st):
    """..."""
    from mpl_toolkits import mplot3d
    import matplotlib.pyplot as plt
    print("|-----------------------------------------------|")
    print("| Plotting ray traces...                        |")

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    cl = ['r', 'b', 'gray', 'dodgerblue']
    for i, tr in enumerate(st):
        ax.plot3D(tr.Xp, tr.Yp, tr.Z, color='dodgerblue', linestyle='dashed',
                  linewidth=1.5,)
        ax.scatter3D(tr.Xp[0], tr.Yp[0], tr.Z[0],
                     c='red', marker='v', edgecolor='k', s=100)
    ax.invert_zaxis()
    plt.show()
    print("|-----------------------------------------------|")
    return


def moho_picker(Gp, xx, zz, migration_param_dict, sta, work_directory, profile):
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


    fontsize = 12
    markersize = 100

    XX, ZZ = np.meshgrid(xx, zz)

    pal_col = work_directory + "/data/colormaps/vik.txt"
    pal_col = pd.read_csv(pal_col, header=None, index_col=False, sep="\s+", names=["R", "G", "B"])
    cm = LinearSegmentedColormap.from_list("blue2red", pal_col.values, len(pal_col))
    c = np.min([np.max(Gp), 0.1])
    c = 0.1
    CL = 1

    plt.close('all')
    # PLOT
    f = plt.figure(1, figsize=[15, 8])
    gs0 = gridspec.GridSpec(nrows=1, ncols=1, figure=f,
                            hspace=0.08, right=0.91, left=0.09, bottom=0.08, top=0.96, )
    ax = f.add_subplot(gs0[0])  # Ray tracing
    ax.scatter(XX, ZZ, c=Gp.T, cmap=cm, s=50, vmin=-c / CL, vmax=c / CL, alpha=.5,
               zorder=1, picker=True, edgecolors=None, marker='8')

    ax.scatter(sta["XSTA"].values, sta["ZSTA"].values,
               markersize, facecolors="grey", edgecolors="k",
               marker="v", lw=0.95, zorder=3, clip_on=False,
               label="Seismic stations", )

    ax.set_aspect("equal")
    ax.set_xlabel("x [km]", fontsize=fontsize)
    ax.set_ylabel("z [km]", fontsize=fontsize)

    majorLocator = MultipleLocator(10)
    minorLocator = MultipleLocator(2.5)
    ax.xaxis.set_major_locator(majorLocator)
    ax.xaxis.set_minor_locator(minorLocator)
    # ax.set_xticks(np.arange(0, 140, 10))

    majorLocator = MultipleLocator(10)
    minorLocator = MultipleLocator(2.5)
    ax.yaxis.set_major_locator(majorLocator)
    ax.yaxis.set_minor_locator(minorLocator)
    ax.set_yticks(np.arange(10, zz[-1], 10))
    ax.set_ylim([80, 0])
    ax.set_xlim([xx[0] - 10 , xx[-1] + 10])

    ax.tick_params(axis="both", which="major", labelsize=fontsize)
    ax.tick_params(axis="both", which="minor", labelsize=fontsize)

    # def onpick(event):
    #     index = event.ind
    #     index = index[0]
    #     xy = event.artist.get_offsets()
    #     print('Dist:', xy[index][0], 'Moho:', xy[index][1])
    #     lat = profile[0][1] + kilometers2degrees(xy[index][0])
    #     print('Lon: ', profile[0][0], 'Lat: ', lat, 'Moho:', xy[index][1], )

    def onkey(event):
        # print(event.key)
        if event.key == 'm':
            if event.xdata is not None and event.ydata is not None:
                print('Dist:', event.xdata, 'Moho:', event.ydata)
                lat = profile[0][1] + kilometers2degrees(event.xdata)
                print('Lon: ', profile[0][0], 'Lat: ', lat, 'Moho:', event.ydata)
                # write moho depths
                with open('moho_depths.txt', 'a') as of:
                    of.write('{}, {}, {}\n'.
                             format(profile[0][0],
                                    lat,event.ydata))
                ax.plot(event.xdata, event.ydata, 'bo-', label='Moho depth')
                f.canvas.draw()
        elif event.key == 'u':
            if event.xdata is not None and event.ydata is not None:
                print('Dist:', event.xdata, 'Uncertain Moho:', event.ydata)
                lat = profile[0][1] + kilometers2degrees(event.xdata)
                print('Lon: ', profile[0][0], 'Lat: ', lat, 'Uncertain Moho:', event.ydata)
                # write moho depths
                with open('unc_moho_depths.txt', 'a') as of:
                    of.write('{}, {}, {}\n'.
                             format(profile[0][0],
                                    lat,event.ydata))
                ax.plot(event.xdata, event.ydata, 'r^-', label='Uncertain Moho')
                f.canvas.draw()
        ax.legend()


    # f.canvas.mpl_connect('pick_event', onkey)
    cid2 = f.canvas.mpl_connect('key_press_event', onkey)

    plt.show()

    return




