"""
Functions for calculating 3D migration of RFs in spherical coordinates.

Location: Chavannes-pres-renens, CH
Date: Mar 2022
Author: Konstantinos Michailos
"""

import numpy as np
from obspy.taup import TauPyModel
from scipy.interpolate import RegularGridInterpolator
from scipy import interpolate
import obspy
import glob
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from scipy import signal
from scipy.interpolate import interp1d
from obspy.geodetics.base import gps2dist_azimuth as gps2dist
from obspy.geodetics import degrees2kilometers, kilometers2degrees


def project(station_lats, station_lons, point_lat, point_lon, angle):
    """
    Projects stations coordinates to a given point (lon, lat) in respect to an angle to the north.

    NOTE: Takes station coordinates and projects them with respect to the center of the profile and the angle
          of the profile with respect to the North direction.
          Output is in [km] for x,y coordinates with respect to lono and lato

    :param station_lats: Seismic station's latitudes in degrees.
    :param station_lons: Seismic station's longitudes in degrees.
    :param point_lat: Given point's latitude in degrees.
    :param point_lon: Given point's longitude in degrees.
    :param angle: Azimuthal angle in degrees.

    :returns: Distance in km parallel and perpendicular to the line.
    """

    ylat = (station_lats - point_lat) * 111.19
    xlon = (station_lons - point_lon) * 111.19 * np.cos(np.radians(station_lats))

    m = np.array([[np.cos(np.radians(angle)), np.sin(np.radians(angle))],
                  [-np.sin(np.radians(angle)), np.cos(np.radians(angle))],])
    r = np.dot(np.column_stack((xlon, ylat)), m)

    distx = r[:, 1]
    disty = r[:, 0]

    return distx, disty


def read_stations_from_sac(path2rfs):
    """
    Read all SAC files to get the seismic site details.

    :type path2rfs: str
    :param path2rfs: Path to the stored RF SAC files.

    :return: Pandas DataFrame of the seismic site details.
    """

    import pandas as pd

    all_rfs = glob.glob(path2rfs + '*.SAC')
    sta_names = []
    sta_lats = []
    sta_lons = []
    sta_eles = []
    for rf in all_rfs:
        trace = obspy.read(rf)
        # This way we only append items for unique stations
        if trace[0].stats.sac.kstnm not in sta_names:
            sta_names.append(trace[0].stats.sac.kstnm)
            sta_lats.append(trace[0].stats.sac.stla)
            sta_lons.append(trace[0].stats.sac.stlo)
            sta_eles.append(trace[0].stats.sac.stel)

    d_ = {}
    d_['NAMESTA'] = sta_names
    d_['LATSTA'] = sta_lats
    d_['LONSTA'] = sta_lons
    d_['ALTSTA'] = sta_eles

    sta = pd.DataFrame(d_)
    # Add elevations
    sta["ZSTA"] = (-1) * sta["ALTSTA"].values / 1000
    print(sta)

    return sta


def read_traces_sphr(path2rfs, sta):
    """
    Read all SAC files to get the seismic site details.

    :type path2rfs: str
    :param path2rfs: Path to the stored RF SAC files.
    :type sta: Pandas DataFrames.
    :param sta: Seismic site details.

    :return: A stream of the receiver function waveforms.
    """

    # Stations x and y values
    lonsta = sta["LONSTA"].values
    latsta = sta["LATSTA"].values
    altsta = sta["ZSTA"].values

    # Assign a unique number to each single separate station (index)
    enum = enumerate(sta["NAMESTA"].values)
    dictionary = dict((i, j) for j, i in enum)
    # Define model (iasp91)
    model = TauPyModel(model="iasp91")

    stream = obspy.Stream()
    all_rfs = glob.glob(path2rfs + '*.SAC')
    print("|-----------------------------------------------|")
    print("| Reading receiver functions...                 |")
    for i, rf in enumerate(all_rfs):
        trace_ = obspy.read(rf)
        trace = trace_[0]
        print(f"| Reading trace {i} of {len(all_rfs)}")
        # Define station's index number
        station_index = dictionary[trace.stats.station]
        trace.station_index = station_index
        # Station's name
        trace.kstnm = trace.stats.sac.kstnm
        # Delta value
        trace.delta = trace.stats.sac.delta
        trace.lon0 = lonsta[station_index]
        trace.lat0 = latsta[station_index]
        # Station's latitude
        trace.stla = trace.stats.sac.stla
        # Station's longitude
        trace.stlo = trace.stats.sac.stlo
        # Station's altitude in km
        trace.alt = sta["ZSTA"].values[station_index]
        # Back azimuth of station to eq
        trace.baz = trace.stats.sac.baz
        trace.stats.baz = trace.stats.sac.baz
        # Local back-azimuth (ori_prof is zero here)
        trace.lbaz = trace.baz
        trace.gcarc = trace.stats.sac.dist
        trace.depth = trace.stats.sac.evdp

        # Should not have distances smaller than 30 degrees...
        if trace.gcarc < 30:
            # Seconds per degrees
            trace.prai = model.get_travel_times(source_depth_in_km=trace.depth,
                                                distance_in_degree=trace.gcarc,
                                                phase_list=["Pn"], )[0].ray_param_sec_degree
        # This is the distance range we use..
        elif trace.gcarc >= 30 and trace.gcarc <= 95.1:
            trace.prai = model.get_travel_times(source_depth_in_km=trace.depth,
                                                distance_in_degree=trace.gcarc,
                                                phase_list=["P"], )[0].ray_param_sec_degree
        # Should not have distances greater than 95 degrees...
        elif trace.gcarc > 95.1:
            try:
                trace.prai = model.get_travel_times(source_depth_in_km=trace.depth,
                                                    distance_in_degree=trace.gcarc,
                                                    phase_list=["PKIKP"], )[0].ray_param_sec_degree
            except Exception as e:
                print(e)
                trace.prai = -10

        trace.time = range(len(trace.data)) * trace.delta - trace.stats.sac.a
        trace.filename = rf
        trace.rms = np.sqrt(np.mean(np.square(trace.data)))

        stream.append(trace)
    print("|-----------------------------------------------|")

    return stream


def get_iasp91(x_, y, z, zmoho):
    """
    Retrieves P-wave, S-wave velocities and depths
    from IASPEI91 global velocity model.

    :type x_: numpy.array
    :param x_: Numpy array of x values of the grid points.
    :type y_: numpy.array
    :param y_: Numpy array of y values of the grid points.
    :type z: numpy.array
    :param z: Numpy array of z values of the grid points.
    :type zmoho: int
    :param zmoho: Moho depth in km.

    :rtype: numpy.ndarrays
    :returns: Array of P-wave, S-wave velocities and their depths.
    """

    R = 6371  # Earth's radius
    x = (R - z) / R
    VP = np.zeros((x_.size, y.size, z.size))
    VS = np.zeros((x_.size, y.size, z.size))
    for i in range(z.size):
        if z[i] < 20:
            VP[:, :, i] = 5.8
            VS[:, :, i] = 3.36
        elif z[i] < zmoho:
            VP[:, :, i] = 6.5
            VS[:, :, i] = 3.75
        elif z[i] < 120.0:
            VP[:, :, i] = 8.78541 - 0.74953 * x[i]
            VS[:, :, i] = 6.706231 - 2.248585 * x[i]
        elif z[i] < 210.0:
            VS[:, :, i] = 5.75020 - 1.27420 * x[i]
            VP[:, :, i] = 25.41389 - 17.69722 * x[i]
        elif z[i] < 410.0:
            VS[:, :, i] = 15.24213 - 11.08552 * x[i]
            VP[:, :, i] = 30.78765 - 23.25415 * x[i]
        elif z[i] < 660.0:
            VP[:, :, i] = 29.38896 - 21.40656 * x[i]
            VS[:, :, i] = 17.70732 - 13.50652 * x[i]
        elif z[i] < 760.0:
            VP[:, :, i] = 25.96984 - 16.93412 * x[i]
            VS[:, :, i] = 20.76890 - 16.53147 * x[i]
        elif z[i] < 2740.0:
            VP[:, :, i] = (25.1486 - 41.1538 * x[i] + 51.9932 * x[i] * x[i] - 26.6083 * x[i] * x[i] * x[i])
            VS[:, :, i] = (12.9303 - 21.2590 * x[i] + 27.8988 * x[i] * x[i] - 14.1080 * x[i] * x[i] * x[i])
        elif z[i] < 2889.0:
            VP[:, :, i] = 14.49470 - 1.47089 * x[i]
            VS[:, :, i] = 8.16616 - 1.58206 * x[i]
        elif z[i] < 5153.9:
            VP[:, :, i] = 10.03904 + 3.75665 * x[i] - 13.67046 * x[i] * x[i]
            VS[:, :, i] = 1.0e-20
        elif z[i] < 6371.0:
            VP[:, :, i] = 11.24094 - 4.09689 * x[i] * x[i]
            VS[:, :, i] = 3.56454 - 3.45241 * x[i] * x[i]
        else:
            VP[:, :, i] = -1
            VS[:, :, i] = -1.0
    return VP, VS


def get_epcrust(min_lon=0, max_lon=25, min_lat=40, max_lat=55):
    """
    Retrieves P-wave, S-wave velocities and depths
    from EPcrust velocity model.

    EPcrust link: http://eurorem.bo.ingv.it/EPcrust_solar/

    :type : numpy.array
    :param : Numpy array of x values of the grid points.

    :returns: 3D interpolation of P-wave, S-wave velocities.
    """

    from scipy.interpolate import LinearNDInterpolator
    import os

    # Read x, y, z, etc .txt file of EPcrust velocity model
    longitudes = []
    latitudes = []
    topos = []
    thick_sediments = []
    thick_upper = []
    thick_lower = []
    vp_sediments = []
    vp_upper = []
    vp_lower = []
    vs_sediments = []
    vs_upper = []
    vs_lower = []

    work_dir = os.getcwd()
    # Path to EPcrust file
    path_epcrust = work_dir + '/data/EPcrust/'
    try:
        with open(path_epcrust + 'EPcrust_0_5.txt', 'r') as f:
            for line in f:
                if line.startswith('#'):
                    print('|Reading EPcrust velocity model...              |')
                    continue
                else:
                    ln = line.split()
                    lon_ = float(ln[0])
                    lat_ = float(ln[1])
                    topo = float(ln[2])
                    thick_sed = float(ln[3])
                    thick_upp = float(ln[4])
                    thick_low = float(ln[5])
                    vp_sed = float(ln[6])
                    vp_upp = float(ln[7])
                    vp_low = float(ln[8])
                    vs_sed = float(ln[9])
                    vs_upp = float(ln[10])
                    vs_low = float(ln[11])
                    if lon_ < max_lon and lon_ > min_lon and lat_ > min_lat and lat_ < max_lat:
                        longitudes.append(lon_)
                        latitudes.append(lat_)
                        topos.append(topo)
                        thick_sediments.append(thick_sed)
                        thick_upper.append(thick_upp)
                        thick_lower.append(thick_low)
                        vp_sediments.append(vp_sed)
                        vp_upper.append(vp_upp)
                        vp_lower.append(vp_low)
                        vs_sediments.append(vs_sed)
                        vs_upper.append(vs_upp)
                        vs_lower.append(vs_low)
    except Exception as e:
        raise type(e)('>>> Move to the top directory of the repository!')

    lon = np.array(longitudes)
    lat = np.array(latitudes)
    ele = np.array(topos)

    # P wave velocities
    vp_sediments = np.array(vp_sediments)
    vp_upper = np.array(vp_upper)
    vp_lower = np.array(vp_lower)
    # S wave velocities
    vs_sediments = np.array(vs_sediments)
    vs_upper = np.array(vs_upper)
    vs_lower = np.array(vs_lower)
    # Thickness of three layers
    thick_sediments = np.array(thick_sediments)
    thick_upper = np.array(thick_upper)
    thick_lower = np.array(thick_lower)

    # Define depth profiles for each EPcrust's grid points
    points = []
    p_velocities = []
    s_velocities = []
    for i, _ in enumerate(lon):

        # Extend velocity model on at least 5 km above sea level
        z_0_ = -5.0
        point0_ = [_, lat[i], z_0_]
        points.append(point0_)
        # Checked that there is always vp and vp sedimentary.
        p_velocities.append(vp_sediments[i])
        s_velocities.append(vs_sediments[i])

        # First point at Earth's surface (including the topography)
        # If we do not include the topography the Moho depth value will be affected.
        elevation = ele[i]
        # If elevation is above sea surface
        if elevation > 0.0:
            # We define the depth profile starting from sea level (z = 0) with positive values as we move deeper.
            # And negative values for the topography above sea level.
            # So when elevation has a positive value we change it to negative.
            z_0 = (-1) * elevation
        elif elevation < 0.0:
            # If elevation is below sea level it will be negative. So similarly to above we multiply with -1
            # to change the value to positive and get the point where the velocity depth profile begins.
            z_0 = (-1) * elevation
        else:
            z_0 = 0.0
        point0 = [_, lat[i], z_0]
        points.append(point0)
        p_velocities.append(vp_sediments[i])
        s_velocities.append(vs_sediments[i])
        # Second point at the lower limit of the sediments.
        z_1 = thick_sediments[i] + z_0
        point1 = [_, lat[i], z_1]
        points.append(point1)
        p_velocities.append(vp_sediments[i])
        s_velocities.append(vs_sediments[i])
        # Third point at the lower limit of the sediments with the velocity below.
        z_2 = thick_sediments[i] + 0.01 + z_0
        point2 = [_, lat[i], z_2]
        points.append(point2)
        p_velocities.append(vp_upper[i])
        s_velocities.append(vs_upper[i])
        # Fourth point at the lower part of the upper crust.
        z_3 = thick_sediments[i] + thick_upper[i] + z_0
        point3 = [_, lat[i], z_3]
        points.append(point3)
        p_velocities.append(vp_upper[i])
        s_velocities.append(vs_upper[i])
        # Fifth point at the lower part of the upper crust...
        z_4 = thick_sediments[i] + thick_upper[i] + 0.01 + z_0
        point4 = [_, lat[i], z_4]
        points.append(point4)
        p_velocities.append(vp_lower[i])
        s_velocities.append(vs_lower[i])
        # Sixth point at the bottom of the crust...
        z_5 = thick_sediments[i] + thick_upper[i] + thick_lower[i] + z_0
        point5 = [_, lat[i], z_5]
        points.append(point5)
        p_velocities.append(vp_lower[i])
        s_velocities.append(vs_lower[i])
        # Seventh point at the bottom of the crust with mantle's velocity
        z_6 = thick_sediments[i] + thick_upper[i] + thick_lower[i] + 0.01 + z_0
        point6 = [_, lat[i], z_6]
        points.append(point6)
        p_velocities.append(8.1)
        s_velocities.append(6.7)
        # Eighth point at the mantle...
        z_7 = 120
        point7 = [_, lat[i], z_7]
        points.append(point7)
        p_velocities.append(8.1)
        s_velocities.append(4.5)

    points = np.array(points)
    values_p = np.array(p_velocities)
    values_s = np.array(s_velocities)
    # rescale here is important for making the steps sharp (look at the following link:
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.LinearNDInterpolator.html
    liner_interpolation_of_velocities_p = LinearNDInterpolator(points, values_p, rescale=True)
    liner_interpolation_of_velocities_s = LinearNDInterpolator(points, values_s, rescale=True)

    return liner_interpolation_of_velocities_p, liner_interpolation_of_velocities_s


def get_end_point(lat1, lon1, baz, d):
    """
    Calculates the end point in lon, lat given we know:
     1) the initial point,
     2) the distance and
     3) the back azimuth

    :param lat1: Starting point's latitude.
    :param lon1: Starting point's longitude.
    :param baz: Back azimuth value in degrees.
    :param d: Length of the line in km.
    :return: Latitude and longitude of the end of the profile.
    """
    # Radius of the Earth
    R = 6371
    # Convert degrees to radians
    back_azimuth_ = np.radians(baz)
    # Current lat point converted to radians
    lat1 = np.radians(lat1)
    # Current long point converted to radians
    lon1 = np.radians(lon1)

    lat2 = np.arcsin(np.sin(lat1) * np.cos(d / R) + np.cos(lat1) * np.sin(d / R) * np.cos(back_azimuth_))
    lon2 = lon1 + np.arctan2(np.sin(back_azimuth_) * np.sin(d / R) * np.cos(lat1),
                             np.cos(d / R) - np.sin(lat1) * np.sin(lat2))
    # Convert to degrees
    lat_2 = np.degrees(lat2)
    lon_2 = np.degrees(lon2)

    return lat_2, lon_2


def tracing_3D_sphr(stream, migration_param_dict, velocity_model='EPcrust'):
    """
    Function to calculate the theoretical ray paths of the receiver functions in spherical coordinates
    in three dimensions.

    :type stream: obspy.core.stream.Stream
    :param stream: Stream of traces.
    :type migration_param_dict: dict
    :param migration_param_dict: Dictionary of grid points for the migration.
    :type zmoho: int
    :param zmoho: Moho depth in km (only used for iasp91 model; currently using EPcrust)

    :returns: Stream of traces that contain the calculated theoretical ray paths.
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
    inc = migration_param_dict['inc']
    zmax = migration_param_dict['zmax']

    # Sanity checks for our grid...
    if maxz < zmax:
        print("Problem: maxz < zmax !!")
        quit()
    if pasz < inc:
        print("Problem: pasz < inc !!")
        quit()

    # Velocity model
    x = np.arange(minx, maxx, pasx)
    y = np.arange(miny, maxy, pasy)
    # Update naming here for ...
    z = np.arange(minz, zmax + (2*minz) + inc, inc)
    # Define the velocity values on each point of the grid
    # EPcrust
    if velocity_model == 'EPcrust':
        P_vel, S_vel = get_epcrust()
    if velocity_model == 'iasp91':
        zmoho = 35
        z_ = np.arange(minz, zmax + inc, inc)
        VP, VS = get_iasp91(x, y, z_, zmoho)
        # Interpolate
        P_vel_3D_grid = RegularGridInterpolator((x, y, z_), VP)
        S_vel_3D_grid = RegularGridInterpolator((x, y, z_), VS)
    if velocity_model != 'EPcrust' and velocity_model != 'iasp91':
        raise IOError('Velocity model should either be EPcrust or iasp91!')

    # Ray tracing
    st = stream.copy()
    st_len = len(st)
    print("|-----------------------------------------------|")
    print("| 3D ray tracing...                             |")
    for i, tr in enumerate(st):
        if tr.prai > -1:
            print('| ' + str(i + 1) + ' of ' + str(st_len))
            # Ray parameter
            p = tr.prai / 111.19
            # Interpolated velocities
            VPinterp = np.zeros(len(z))
            VSinterp = np.zeros(len(z))
            # S-ray parameter at surface longitude
            Xs = np.zeros(len(z))
            Xs[0] = tr.lon0
            # S-ray parameter at surface latitude
            Ys = np.zeros(len(z))
            Ys[0] = tr.lat0
            # P-ray parameter at surface longitude
            Xp = np.zeros(len(z))
            Xp[0] = tr.lon0
            # P-ray parameter at surface latitude
            Yp = np.zeros(len(z))
            Yp[0] = tr.lat0
            # Back azimuth
            baz_p = np.zeros(len(z))
            _, _, baz_p[0] = gps2dist(tr.stats.sac.evla, tr.stats.sac.evlo, tr.lat0, tr.lon0)
            baz_s = np.zeros(len(z))
            _, _, baz_s[0] = gps2dist(tr.stats.sac.evla, tr.stats.sac.evlo, tr.lat0, tr.lon0)
            # Time it takes for waves to travel through each layer
            Tp = np.zeros(len(z))
            Ts = np.zeros(len(z))
            # -------------------------------
            # Migrate with 3-D velocity model
            # -------------------------------
            for iz in range(len(z) - 1):
                # Apply correction -1 * minz to move to sea level and from there add tr.alt to begin
                # from the stations elevation and not from the top of the grid (zmin).
                z_sta = z[iz] + (-1) * minz + tr.alt

                pts = np.array([Xp[iz], Yp[iz], z_sta])
                # IASP91
                if velocity_model == 'iasp91':
                    VPinterp[iz] = P_vel_3D_grid(pts)
                    # print(z[iz], VPinterp[iz])
                # EPcrust
                if velocity_model == 'EPcrust':
                    VPinterp[iz] = P_vel(pts)[0]
                    # print(z[iz], VPinterp[iz])
                r_earth = 6371
                # Calculate departing incidence angle of the ray (p = r_earth * sin(incidence_angle) / V)
                id_p = np.arcsin(p * VPinterp[iz])
                id_degrees_p = np.rad2deg(id_p)
                # Calculate great - circle distance travelled delta_i - 1 (delta)
                ia_i_p = np.arcsin((np.sin(id_p)) / (r_earth - (z[iz+1] + (-1) * minz + tr.alt)) * (r_earth - (z[iz] + (-1) * minz + tr.alt)))
                # 180 - this
                ia_i_degrees_p = 180 - np.rad2deg(ia_i_p)
                # Angle
                delta_p = 180 - id_degrees_p - ia_i_degrees_p
                # Distance from A to B in km
                gc_dist_p = np.radians(delta_p) * (r_earth - (z[iz] + (-1) * minz + tr.alt))
                # Location of B
                lat_2, lon_2 = get_end_point(Yp[iz], Xp[iz], baz_p[iz], gc_dist_p )
                Yp[iz + 1] = lat_2
                Xp[iz + 1] = lon_2
                _, _, baz_p[iz + 1] = gps2dist(tr.stats.sac.evla, tr.stats.sac.evlo, Yp[iz + 1], Xp[iz + 1], )
                Tp[iz + 1] = Tp[iz] + (inc / np.cos(id_p)) / VPinterp[iz]

                # Same as above for S wave
                # IASP91
                if velocity_model == 'iasp91':
                    VSinterp[iz] = S_vel_3D_grid(pts)
                # EPcrust
                if velocity_model == 'EPcrust':
                    VSinterp[iz] = S_vel(pts)[0]
                # Calculate departing incidence angle of the ray (p = r_earth * sin(incidence_angle) / V)
                ################################33
                # with open('/home/kmichailos/Desktop/' + str(i) +velocity_model +'.txt', 'a') as of:
                #     of.write('{}, {}\n'.
                #              format(z_sta, VSinterp[iz]))
                ################################33
                id_s = np.arcsin(p * VSinterp[iz])
                id_degrees_s = np.rad2deg(id_s)
                # Calculate great - circle distance travelled delta_i - 1 (delta)
                ia_i_s = np.arcsin(
                    (np.sin(id_s)) / (r_earth - (z[iz + 1] + (-1) * minz + tr.alt)) * (r_earth - (z[iz] + (-1) * minz + tr.alt)))
                # 180 - this
                ia_i_degrees_s = 180 - np.rad2deg(ia_i_s)
                # Angle
                delta_s = 180 - id_degrees_s - ia_i_degrees_s
                # Distance from A to B in km
                gc_dist_s = np.radians(delta_s) * (r_earth - (z[iz] + (-1) * minz + tr.alt))

                lat_2, lon_2 = get_end_point(Ys[iz], Xs[iz], baz_s[iz], gc_dist_s)
                Ys[iz + 1] = lat_2
                Xs[iz + 1] = lon_2
                _, _, baz_s[iz + 1] = gps2dist(tr.stats.sac.evla, tr.stats.sac.evlo, Ys[iz + 1], Xs[iz + 1], )
                Ts[iz + 1] = Ts[iz] + (inc / np.cos(id_s)) / VSinterp[iz]

            # ____________________end of 3D migration_______
            D = np.sqrt(np.square(Xp - Xs) + np.square(Yp - Ys))
            E = np.sqrt(np.square(Xp - Xp[0]) + np.square(Yp - Yp[0]))

            Td = D * p
            Te = 2 * E * p

            tr.Z = z + (-1) * minz + tr.alt
            tr.Xp = Xp
            tr.Yp = Yp
            tr.Xs = Xs
            tr.Ys = Ys
            tr.Ts = Ts
            tr.Tp = Tp

            interp = interp1d(tr.time, tr.data, bounds_error=False, fill_value=np.nan)
            tps = -Tp + Ts + Td
            tpps = Tp + Ts + Td - Te
            tpss = 2 * Ts + 2 * Td - Te
            tr.amp_ps = interp(tps)
            tr.amp_pps = interp(tpps)
            tr.amp_pss = interp(tpss)

            # Theoretical traces
            interp = interpolate.interp1d(tpps, tr.amp_ps, bounds_error=False, fill_value=np.nan)
            tr.amp_pps_theo = interp(tps)
            interp = interpolate.interp1d(tpss, tr.amp_ps, bounds_error=False, fill_value=np.nan)
            tr.amp_pss_theo = interp(tps)

            tr.tps = tps
            tr.tpps = tpps
            tr.tpss = tpss

        else:
            print("prai: ", tr.prai)
            tr.Xp = -1
            tr.Yp = -1
            tr.Xs = -1
            tr.Ys = -1
            tr.Z = -1
            tr.amp_ps = -1
            tr.amp_pps = -1
            tr.amp_pss = -1
    print("| End of 3D ray tracing...                      |")
    print("|-----------------------------------------------|")

    return st


# todo: docstring
def ccpm_3d(stream, migration_param_dict, output_file, phase="PS"):
    """
    ...

    :type stream: obspy.core.stream.Stream
    :param stream: Stream of traces.
    :type migration_param_dict: dict
    :param migration_param_dict: Dictionary of grid points for the migration.
    :type output_file:
    :param output_file:
    :type phase:
    :param phase:

    :returns:
    """

    # Time to depth Migration
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

    # Grid preparation
    x = np.arange(minx, maxx + pasx, pasx)
    y = np.arange(miny, maxy + pasy, pasy)
    z = np.arange(minz, maxz + pasz, pasz)

    # rms selection
    rms = np.zeros(len(stream), dtype="float")
    for k in range(len(stream)):
        rms[k] = stream[k].rms
    rms = np.sort(rms)
    i_rms = int(np.floor(len(stream) * 0.98) - 1)
    rms_max = rms[i_rms]

    # Amplitude matrix
    G = np.zeros((len(x), len(y), len(z)))
    # Number of amplitudes in each cell of the matrix
    nG = np.zeros((len(x), len(y), len(z))) + 1e-8
    # Longitudes and latitudes of the matrix

    print("|-----------------------------------------------|")
    print("| Start of common conversion point stacking...  |")

    for i, tr in enumerate(stream):
        if tr.prai >-1 and tr.rms <= rms_max:
            print('| ' + str(i + 1) + ' of ' + str(len(stream)))
            # Find box that the trace is in
            ix = np.floor((tr.Xs - minx) / pasx)
            iy = np.floor((tr.Ys - miny) / pasy)
            iz = np.floor((tr.Z - minz) / pasz)
            ix = np.array(ix, dtype="int")
            iy = np.array(iy, dtype="int")
            iz = np.array(iz, dtype="int")
            for iz_ in iz:
                if iz_ > len(z)-1:
                    print(iz_)
            if phase == "PS":
                # Adding up all the elements
                G[ix, iy, iz] = G[ix, iy, iz] + tr.amp_ps
            # Need to figure out where this i1z came from?
            # Since we only use PS phases we won't need them.
            elif phase == "PPS":
                G[ix, iy, iz] = G[ix, iy, iz] + stream[i].amp_pps[:i1z]
            elif phase == "PSS":
                G[ix, iy, iz] = G[ix, iy, iz] - stream[i].amp_pss[:i1z]
            # Number of observations in each cell
            nG[ix, iy, iz] = nG[ix, iy, iz] + 1
        else:
            print(f'Removing trace, {tr}, because of high rms-value.')
    # G2 = np.squeeze((np.sum(G, axis=1)))
    # nG2 = np.squeeze((np.sum(nG, axis=1)))
    # G2 = G2 / nG2
    # Get the average number of the amplitudes
    G = np.divide(G, nG)

    print("| End of common conversion point stacking...    |")
    np.save(output_file, G, allow_pickle=False)
    print("|-----------------------------------------------|")

    return G


def ccp_smooth(G2, migration_param_dict):
    """
    Apply smoothing to the migrated RF image.

    :type G2:
    :param G2:
    :type migration_param_dict: dict
    :param migration_param_dict: Dictionary of grid points for the migration.

    :returns: Smoothed amplitudes in 2D grid.
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
    inc = migration_param_dict['inc']
    zmax = migration_param_dict['zmax']

    zz = np.arange(minz, maxz + pasz, pasz)
    # NOTE: l0, dl are parameters to decide from
    #       which depth the smoothing becomes important and how much

    zbegin_lisse = -2
    # pasx is in degrees so we modify the line below
    l0 = 1
    l0 = 1./111.11
    dl = 1000000
    # dl = 100

    with np.errstate(divide="warn"):
        G3 = G2
        for iz in range(G2.shape[1]):
            if zz[iz] < zbegin_lisse:
                G3[:, iz] = G2[:, iz]
            else:
                # sigmal = (zz[iz] / dl + l0) / (pasx *111)
                sigmal = (l0) / pasx
                print(sigmal)

                nbml = G2.shape[0]

                mm = int(np.round(nbml / 2))
                C = (
                    1
                    / (sigmal * np.sqrt(2 * np.pi))
                    * np.exp(-0.5 * np.square(np.arange(nbml) - mm) / (sigmal * sigmal))
                )
                C = C / np.sum(C)
                temp = np.convolve(G2[:, iz], C)
                G3[:, iz] = temp[mm : mm + G2.shape[0]]

    return G3


def ccpFilter(G2):
    """
    Convolution with a Gaussian bell for local smooth.

    :type G2:
    :param G2:

    :returns:
    """

    nbm = 7
    b, a = 3, 1.5
    sigma = 2
    C = np.zeros((nbm, nbm))
    mm = np.floor(nbm / 2)
    for i in range(nbm):
        for j in range(nbm):
            r = np.sqrt(((i - mm) ** 2) / (a ** 2) + ((j - mm) ** 2 / (b ** 2)))
            if r < mm:
                C[i, j] = np.exp(-0.5 * (r / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))
    C = C / np.sum(C)
    miniG2 = signal.convolve2d(G2, C, "same")
    return miniG2
