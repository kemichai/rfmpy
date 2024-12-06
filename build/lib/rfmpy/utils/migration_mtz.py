import matplotlib
matplotlib.use('TkAgg')
import platform
import os
import time
import numpy as np
from obspy.taup import TauPyModel
from scipy.interpolate import RegularGridInterpolator
from scipy import interpolate
import obspy
import glob
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from obspy.geodetics.base import gps2dist_azimuth as gps2dist
from obspy.geodetics import degrees2kilometers, kilometers2degrees
from obspy import Stream
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
        trace.gcarc = trace.stats.sac.gcarc
        trace.depth = trace.stats.sac.evdp

        # Should not have distances smaller than 30 degrees...
        if trace.gcarc < 28:
            # Seconds per degrees
            try:
                trace.prai = model.get_travel_times(source_depth_in_km=trace.depth,
                                                distance_in_degree=trace.gcarc,
                                                phase_list=["Pn"], )[0].ray_param_sec_degree
            except Exception as e:
                print(e)
                trace.prai = -10

        # This is the distance range we use..
        elif trace.gcarc >= 28 and trace.gcarc <= 95.1:
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

        trace.time = range(len(trace.data)) * trace.delta - 5
        trace.filename = rf
        trace.rms = np.sqrt(np.mean(np.square(trace.data)))

        stream.append(trace)
    print("|-----------------------------------------------|")

    return stream


def write_files_4_piercing_points_and_raypaths(st, sta, piercing_depth=35, plot=True):
    """..."""

    piercing_lon = []
    piercing_lat = []
    for i, tr in enumerate(st):
        # this try statement gets rid of RFs with no ray paths (ones with tr.Z = -1; tr.Xp = -1 etc)
        try:
            depth = tr.Z[0]
        except:
            continue
        for j, z in enumerate(tr.Z):
            if z > piercing_depth and z < piercing_depth + 15:
                piercing_lon.append(tr.Xp[j])
                piercing_lat.append(tr.Yp[j])
                with open('../piercing_points.txt', 'a') as of:
                    of.write('{}, {}\n'.format(tr.Xp[j], tr.Yp[j]))

    wav_p_lon = []
    wav_p_lat = []
    wav_p_dep = []
    for i, tr in enumerate(st):
        try:
            depth = tr.Z[0]
        except:
            continue
        for j, z in enumerate(tr.Z):
                wav_p_lon.append(tr.Xp[j])
                wav_p_lat.append(tr.Yp[j])
                wav_p_dep.append(z)
                with open('../ray_path.txt', 'a') as of:
                    of.write('{}, {}, {}\n'.
                            format(tr.Xp[j], tr.Yp[j], z))

    if plot:
        # Plot raypaths
        plt.scatter(wav_p_lon, wav_p_lat, alpha=0.5,
                    c=wav_p_dep, marker='.', edgecolor=None, s=5)
        plt.scatter(sta["LONSTA"], sta["LATSTA"],
                    c='r', marker='v', edgecolor='k', s=100)
        plt.show()

        # Plot piercing points
        plt.scatter(piercing_lon, piercing_lat, alpha=.3,
                    c='gray', marker='x', edgecolor='gray', s=50)
        plt.scatter(sta["LONSTA"], sta["LATSTA"],
                    c='r', marker='v', edgecolor='k', s=100)
        plt.show()

    return


def read_vel_model(migration_param_dict, velocity_model='zmodel_m60'):
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

    # EPcrust
    if velocity_model == 'EPcrust':
        P_vel, S_vel = get_epcrust()
    elif velocity_model == 'zmodel_m60':
        P_vel, S_vel = get_zmodel_m60()
    if velocity_model == 'iasp91':
        zmoho = 35
        z_ = np.arange(minz, zmax + inc, inc)
        VP, VS = get_iasp91(x, y, z_, zmoho)
        # Interpolate
        P_vel = RegularGridInterpolator((x, y, z_), VP)
        S_vel = RegularGridInterpolator((x, y, z_), VS)
    if velocity_model != 'EPcrust' and velocity_model != 'iasp91' and velocity_model != 'zmodel_m60':
        raise IOError('Velocity model should either be EPcrust, iasp91 or zmodel_m60!')
    return P_vel, S_vel




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
            try:
                sta_names.append(trace[0].stats.sac.kstnm)
                sta_lats.append(trace[0].stats.sac.stla)
                sta_lons.append(trace[0].stats.sac.stlo)
                sta_eles.append(trace[0].stats.sac.stel)
            except Exception as e:
                print(e)
                print(f"Station elevation missing from trace {rf}")
                continue

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


#def read_traces_sphr(path2rfs, sta):
    """
    Read all SAC files to get the seismic site details.

    :type path2rfs: str
    :param path2rfs: Path to the stored RF SAC files.
    :type sta: Pandas DataFrames.
    :param sta: Seismic site details.

    :return: A stream of the receiver function waveforms.
    """

    # Stations x and y values
#    lonsta = sta["LONSTA"].values
#    latsta = sta["LATSTA"].values
#    altsta = sta["ZSTA"].values

    # Assign a unique number to each single separate station (index)
#    enum = enumerate(sta["NAMESTA"].values)
#    dictionary = dict((i, j) for j, i in enum)
#    # Define model (iasp91)
#    model = TauPyModel(model="iasp91")

#    stream = obspy.Stream()
#    all_rfs = glob.glob(path2rfs + '*.SAC')
#    print("|-----------------------------------------------|")
#    print("| Reading receiver functions...                 |")
#    for i, rf in enumerate(all_rfs):
#        trace_ = obspy.read(rf)
#        trace = trace_[0]
#        print(f"| Reading trace {i} of {len(all_rfs)}")
#        # Define station's index number
#        station_index = dictionary[trace.stats.station]
#        trace.station_index = station_index
        # Station's name
#        trace.kstnm = trace.stats.sac.kstnm
        # Delta value
#        trace.delta = trace.stats.sac.delta
#        trace.lon0 = lonsta[station_index]
#        trace.lat0 = latsta[station_index]
        # Station's latitude
#        trace.stla = trace.stats.sac.stla
        # Station's longitude
#        trace.stlo = trace.stats.sac.stlo
        # Station's altitude in km
#        trace.alt = sta["ZSTA"].values[station_index]
        # Back azimuth of station to eq
#        trace.baz = trace.stats.sac.baz
#        trace.stats.baz = trace.stats.sac.baz
        # Local back-azimuth (ori_prof is zero here)
#        trace.lbaz = trace.baz
#        trace.gcarc = trace.stats.sac.dist
#        trace.depth = trace.stats.sac.evdp

        # Should not have distances smaller than 30 degrees...
#        if trace.gcarc < 30:
            # Seconds per degrees
#            trace.prai = model.get_travel_times(source_depth_in_km=trace.depth,
#                                                distance_in_degree=trace.gcarc,
#                                                phase_list=["Pn"], )[0].ray_param_sec_degree
        # This is the distance range we use..
#        elif trace.gcarc >= 30 and trace.gcarc <= 95.1:
#            trace.prai = model.get_travel_times(source_depth_in_km=trace.depth,
#                                                distance_in_degree=trace.gcarc,
#                                                phase_list=["P"], )[0].ray_param_sec_degree
        # Should not have distances greater than 95 degrees...
#        elif trace.gcarc > 95.1:
#            try:
#                trace.prai = model.get_travel_times(source_depth_in_km=trace.depth,
#                                                    distance_in_degree=trace.gcarc,
#                                                    phase_list=["PKIKP"], )[0].ray_param_sec_degree
#            except Exception as e:
#                print(e)
#                trace.prai = -10

#        trace.time = range(len(trace.data)) * trace.delta - trace.stats.sac.a
#        trace.filename = rf
#        trace.rms = np.sqrt(np.mean(np.square(trace.data)))

#        stream.append(trace)
#    print("|-----------------------------------------------|")

#    return stream


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
        # iasp91 values
        # p_velocities.append(8.05)
        # s_velocities.append(4.45)
        # epcrust lower crust
        p_velocities.append(vp_lower[i])
        s_velocities.append(vs_lower[i])
        # Eighth point at the mantle...
        z_7 = 120
        point7 = [_, lat[i], z_7]
        points.append(point7)
        # iasp91 values
        # p_velocities.append(8.1)
        # s_velocities.append(4.5)
        # epcrust lower crust
        p_velocities.append(vp_lower[i] + 0.05)
        s_velocities.append(vs_lower[i] + 0.05)

    points = np.array(points)
    values_p = np.array(p_velocities)
    values_s = np.array(s_velocities)
    # rescale here is important for making the steps sharp (look at the following link:
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.LinearNDInterpolator.html
    liner_interpolation_of_velocities_p = LinearNDInterpolator(points, values_p, rescale=True)
    liner_interpolation_of_velocities_s = LinearNDInterpolator(points, values_s, rescale=True)

    return liner_interpolation_of_velocities_p, liner_interpolation_of_velocities_s


def get_zmodel_m60(min_lon=0, max_lon=32, min_lat=40, max_lat=55):
    """
    Retrieves P-wave, S-wave velocities and depths
    from ZMODEL_M60 velocity model (Zhu et al., 2015).

    Model Format: The model file ZMODEL_M60.dat is arranged as:
    Long(deg) / Lat(deg) / Depth(km) / Iso Vp(km/s) / Iso Vp perturbation(%) /
    Iso Vs(km/s) / Iso Vs perturbation(%) / radial anisotropy / Q value /

    Link: https://academic.oup.com/gji/article/201/1/18/724841#86405283

    :type : numpy.array
    :param : Numpy array of x values of the grid points.

    :returns: 3D interpolation of P-wave, S-wave velocities.
    """

    from scipy.interpolate import LinearNDInterpolator
    import os

    work_dir = os.getcwd()

    # Path to file

    path_zmodel_m60 = work_dir + '/data/ZMODEL_M60/'

    # Read x, y, z, etc .nps file of ZMODEL_M60 velocity model
    parameters = np.load(path_zmodel_m60 + 'z_model_m60.npz')
    parameters.items()
    longitudes = parameters["longitudes"].tolist()
    latitudes = parameters["latitudes"].tolist()
    depths = parameters["depths"].tolist()
    p_velocities = parameters["Vp"].tolist()
    s_velocities = parameters["Vs"].tolist()

    points_list = []
    vp_values = []
    vs_values = []
    for i, lon in enumerate(longitudes):
        if lon < max_lon and lon > min_lon and latitudes[i] > min_lat and latitudes[i] < max_lat:

            point = [lon, latitudes[i], depths[i]]
            if depths[i] == 10.0:
                # Extend velocity model on at least 5 km above sea level
                # IMPORTANT TO NOTE
                additional_point = [lon, latitudes[i], -5]
                points_list.append(additional_point)
                vp_values.append(p_velocities[i])
                vs_values.append(s_velocities[i])

            points_list.append(point)
            vp_values.append(p_velocities[i])
            vs_values.append(s_velocities[i])

    points = np.array(points_list)
    values_p = np.array(vp_values)
    values_s = np.array(vs_values)
    # rescale here is important for making the steps sharp (look at the following link:
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.LinearNDInterpolator.html
    print("|-----------------------------------------------|")
    print("| Interpolating 3D ZMODEL_M60...                |")
    print("| This might take a while...                    |")
    liner_interpolation_of_velocities_p = LinearNDInterpolator(points, values_p, rescale=True)
    liner_interpolation_of_velocities_s = LinearNDInterpolator(points, values_s, rescale=True)
    print("| Interpolated 3D ZMODEL_M60...                 |")

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




def tracing_3D_sphr_parallel(stream, migration_param_dict, P_vel, S_vel):
    """
    Function to calculate the theoretical ray paths of the receiver functions in spherical coordinates
    in three dimensions. Latest version uses multiple cpus running in parallel.

    :type stream: obspy.core.stream.Stream
    :param stream: Stream of traces.
    :type migration_param_dict: dict
    :param migration_param_dict: Dictionary of grid points for the migration.
    :type
    :param

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
    # TODO: add extra option here
    # Define the velocity values on each point of the grid

    # # EPcrust
    # if velocity_model == 'EPcrust':
    #     P_vel, S_vel = get_epcrust()
    # elif velocity_model == 'zmodel_m60':
    #     P_vel, S_vel = get_zmodel_m60()
    # if velocity_model == 'iasp91':
    #     zmoho = 35
    #     z_ = np.arange(minz, zmax + inc, inc)
    #     VP, VS = get_iasp91(x, y, z_, zmoho)
    #     # Interpolate
    #     P_vel_3D_grid = RegularGridInterpolator((x, y, z_), VP)
    #     S_vel_3D_grid = RegularGridInterpolator((x, y, z_), VS)
    # if velocity_model != 'EPcrust' and velocity_model != 'iasp91' and velocity_model != 'zmodel_m60':
    #     raise IOError('Velocity model should either be EPcrust, iasp91 or zmodel_m60!')

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
                # # IASP91
                # if velocity_model == 'iasp91':
                #     VPinterp[iz] = P_vel_3D_grid(pts)
                #     # print(z[iz], VPinterp[iz])
                # # EPcrust
                # if velocity_model == 'EPcrust':
                #     VPinterp[iz] = P_vel(pts)[0]
                #     # print(z[iz], VPinterp[iz])
                # # zmodel_m60
                # if velocity_model == 'zmodel_m60':
                #     VPinterp[iz] = P_vel(pts)[0]
                #     # print(z[iz], VPinterp[iz])
                VPinterp[iz] = P_vel(pts)[0]


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
                #Tp[iz + 1] = Tp[iz] + (inc / np.cos(id_p)) / VPinterp[iz]  #### WRONG OLD ######
                rlinep=np.sin(np.radians(delta_p))*(r_earth-z[iz+1])/np.sin(id_p)
                Tp[iz + 1] = Tp[iz] + rlinep / VPinterp[iz]

                # Same as above for S wave
                # IASP91
                # if velocity_model == 'iasp91':
                #     VSinterp[iz] = S_vel_3D_grid(pts)
                    # print(z[iz], VSinterp[iz])

                # EPcrust
                # if velocity_model == 'EPcrust':
                #     VSinterp[iz] = S_vel(pts)[0]
                    # print(z[iz], VSinterp[iz])
                # zmodel_m60
                # if velocity_model == 'zmodel_m60':
                #     VSinterp[iz] = S_vel(pts)[0]
                    # print(z[iz], VPinterp[iz])
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
                #Ts[iz + 1] = Ts[iz] + (inc / np.cos(id_s)) / VSinterp[iz]  #### WRONG OLD ######
                rlines=np.sin(np.radians(delta_s))*(r_earth-z[iz+1])/np.sin(id_s)
                Ts[iz + 1] = Ts[iz] + rlines / VSinterp[iz]

            # ____________________end of 3D migration_______
            #D_r = np.sqrt(np.square(Xp - Xs) + np.square(Yp - Ys))  #### WRONG OLD ######
            #D_r = (np.arccos(np.cos(np.radians(Yp))*np.cos(np.radians(Ys))*np.cos(np.radians(Xp))-np.radians(Xs)+np.sin(np.radians(Yp))*np.sin(np.radians(Ys))))*(r_earth-z[iz])
            A= (np.cos(np.radians(Yp)))*(np.cos(np.radians(Ys)))*(np.cos(np.radians(Xp-Xs)))
            B= (np.sin(np.radians(Yp)))*(np.sin(np.radians(Ys)))
            D = (np.arccos(A+B))*((np.radians(360)*r_earth-z[iz])/np.radians(360))
            E = np.sqrt(np.square(Xp - Xp[0]) + np.square(Yp - Yp[0])) #####NOT USED######

            Td = D * p
            #Te = 2 * E * p

            tr.Z = z + (-1) * minz + tr.alt
            tr.Xp = Xp
            tr.Yp = Yp
            tr.Xs = Xs
            tr.Ys = Ys
            tr.Ts = Ts
            tr.Tp = Tp

            interp = interp1d(tr.time, tr.data, bounds_error=False, fill_value=np.nan)
            tps = -Tp + Ts + Td
            #tpps = Tp + Ts + Td - Te
            #tpss = 2 * Ts + 2 * Td - Te
            tr.amp_ps = interp(tps)
            #tr.amp_pps = interp(tpps)
            #tr.amp_pss = interp(tpss)

            # Theoretical traces
            #interp = interpolate.interp1d(tpps, tr.amp_ps, bounds_error=False, fill_value=np.nan)
            #tr.amp_pps_theo = interp(tps)
            #interp = interpolate.interp1d(tpss, tr.amp_ps, bounds_error=False, fill_value=np.nan)
            #tr.amp_pss_theo = interp(tps)

            tr.tps = tps
            #tr.tpps = tpps
            #tr.tpss = tpss

        else:
            print("prai: ", tr.prai)
            tr.Xp = -1
            tr.Yp = -1
            tr.Xs = -1
            tr.Ys = -1
            tr.Z = -1
            tr.amp_ps = -1
            #tr.amp_pps = -1
            #tr.amp_pss = -1
    print("| End of 3D ray tracing...                      |")
    print("|-----------------------------------------------|")

    return st

import numpy as np
import obspy
from obspy import Trace, Stream
from obspy.geodetics import gps2dist_azimuth
import matplotlib.pyplot as plt

def tracing_3D_sphr1(stream, migration_param_dict, velocity_model):
    """
    Function to calculate the theoretical ray paths of the receiver functions in spherical coordinates
    in three dimensions.

    :type stream: obspy.core.stream.Stream
    :param stream: Stream of traces.
    :type migration_param_dict: dict
    :param migration_param_dict: Dictionary of grid points for the migration.
    :type velocity_model: str
    :param velocity model: Velocity model to be used (options include 'iasp91','zmodel_m60' and 'EPcrust')

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

    # Velocity model
    x = np.arange(minx, maxx, pasx)
    y = np.arange(miny, maxy, pasy)
    z = np.arange(minz, zmax + (2*minz) + inc, inc)

    # Initialize an empty Stream to store the ray paths
    ray_path_stream = Stream()

    # Ray tracing
    for tr in stream:
        if tr.prai > -1:
            # Ray parameter
            p = tr.prai / 111.19

            # Initialize arrays to store ray path coordinates
            latitudes = []
            longitudes = []
            depths = []

            # Starting point coordinates
            start_lat = tr.stats.sac.evla
            start_lon = tr.stats.sac.evlo
            if start_lat < -90:
                start_lat = -90
            elif start_lat > 90:
                start_lat = 90

            latitudes.append(start_lat)
            longitudes.append(start_lon)
            depths.append(0)

            # Ray tracing loop
            for iz in range(len(z) - 1):
                # Calculate arrival coordinates
                arrival_lat, arrival_lon, _ = gps2dist_azimuth(
                    latitudes[-1], longitudes[-1], tr.stats.sac.evla, tr.stats.sac.evlo)

                # Ensure arrival latitude is within valid range
                if arrival_lat < -90:
                    arrival_lat = -90
                elif arrival_lat > 90:
                    arrival_lat = 90

                # Append arrival coordinates
                latitudes.append(arrival_lat)
                longitudes.append(arrival_lon)
                depths.append(z[iz])

            # Create Trace objects containing the ray path coordinates
            lat_trace = Trace(data=np.array(latitudes))
            lon_trace = Trace(data=np.array(longitudes))
            depth_trace = Trace(data=np.array(depths))

            # Set the stats for each Trace
            lat_trace.stats.station = tr.stats.station
            lat_trace.stats.channel = "latitude"
            lon_trace.stats.station = tr.stats.station
            lon_trace.stats.channel = "longitude"
            depth_trace.stats.station = tr.stats.station
            depth_trace.stats.channel = "depth"

            # Append Trace objects to Stream
            ray_path_stream.append(lat_trace)
            ray_path_stream.append(lon_trace)
            ray_path_stream.append(depth_trace)

    return ray_path_stream

def tracing_3D_sphr(stream, migration_param_dict, velocity_model):
    """
    Function to calculate the theoretical ray paths of the receiver functions in spherical coordinates
    in three dimensions. Original version uses 1 cpu for each calculation.

    :type stream: obspy.core.stream.Stream
    :param stream: Stream of traces.
    :type migration_param_dict: dict
    :param migration_param_dict: Dictionary of grid points for the migration.
    :type velocity_model: str
    :param velocity model: Velocity model to be used (options include 'iasp91','zmodel_m60' and 'EPcrust')

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
    z = np.arange(minz, zmax + (2*minz) + inc, inc) ###!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # TODO: add extra option here
    # Define the velocity values on each point of the grid

    # # EPcrust
    if velocity_model == 'EPcrust':
        P_vel, S_vel = get_epcrust()
    elif velocity_model == 'zmodel_m60':
        P_vel, S_vel = get_zmodel_m60()
    if velocity_model == 'iasp91':
        zmoho = 45
        z_ = np.arange(minz, zmax + inc, inc)
        VP, VS = get_iasp91(x, y, z_, zmoho)
        # Interpolate
        P_vel = RegularGridInterpolator((x, y, z_), VP)
        S_vel = RegularGridInterpolator((x, y, z_), VS)
    if velocity_model != 'EPcrust' and velocity_model != 'iasp91' and velocity_model != 'zmodel_m60':
        raise IOError('Velocity model should either be EPcrust, iasp91 or zmodel_m60!')

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
                # # IASP91
                if velocity_model == 'iasp91':
                    # print(pts)
                    VPinterp[iz] = P_vel(pts)
                # # EPcrust
                if velocity_model == 'EPcrust':
                    VPinterp[iz] = P_vel(pts)[0]
                #     # print(z[iz], VPinterp[iz])
                # # zmodel_m60
                if velocity_model == 'zmodel_m60':
                    VPinterp[iz] = P_vel(pts)[0]
                #     # print(z[iz], VPinterp[iz])
                VPinterp[iz] = P_vel(pts)[0]

                #print (VPinterp[iz])
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
                #Tp[iz + 1] = Tp[iz] + (inc / np.cos(id_p)) / VPinterp[iz]  #### WRONG OLD ######
                rlinep=np.sin(np.radians(delta_p))*(r_earth-z[iz+1])/np.sin(id_p)
                Tp[iz + 1] = Tp[iz] + rlinep / VPinterp[iz]
                # Same as above for S wave
                # IASP91
                if velocity_model == 'iasp91':
                    VSinterp[iz] = S_vel(pts)
                    # print(z[iz], VSinterp[iz])

                # EPcrust
                if velocity_model == 'EPcrust':
                    VSinterp[iz] = S_vel(pts)[0]
                    # print(z[iz], VSinterp[iz])
                # zmodel_m60
                if velocity_model == 'zmodel_m60':
                    VSinterp[iz] = S_vel(pts)[0]
                    # print(z[iz], VPinterp[iz])
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
                #Ts[iz + 1] = Ts[iz] + (inc / np.cos(id_s)) / VSinterp[iz]  #### WRONG OLD ######
                rlines=np.sin(np.radians(delta_s))*(r_earth-z[iz+1])/np.sin(id_s)
                Ts[iz + 1] = Ts[iz] + rlines / VSinterp[iz]
            # ____________________end of 3D migration_______
            #D_r = np.sqrt(np.square(Xp - Xs) + np.square(Yp - Ys))  #### WRONG OLD ######
            #D_r = (np.arccos(np.cos(np.radians(Yp))*np.cos(np.radians(Ys))*np.cos(np.radians(Xp))-np.radians(Xs)+np.sin(np.radians(Yp))*np.sin(np.radians(Ys))))*(r_earth-z[iz])
            A= (np.cos(np.radians(Yp)))*(np.cos(np.radians(Ys)))*(np.cos(np.radians(Xp-Xs)))
            B= (np.sin(np.radians(Yp)))*(np.sin(np.radians(Ys)))
            D = (np.arccos(A+B))*((np.radians(360)*r_earth-z[iz])/np.radians(360))
            #E = np.sqrt(np.square(Xp - Xp[0]) + np.square(Yp - Yp[0])) #####NOT USED######
            Td = D * p
            #Te = 2 * E * p

            tr.Z = z + (-1) * minz + tr.alt
            tr.Xp = Xp
            tr.Yp = Yp
            tr.Xs = Xs
            tr.Ys = Ys
            tr.Ts = Ts
            tr.Tp = Tp

            interp = interp1d(tr.time, tr.data, bounds_error=False, fill_value=np.nan)
            tps = -Tp + Ts + Td
            #tpps = Tp + Ts + Td - Te
            #tpss = 2 * Ts + 2 * Td - Te
            tr.amp_ps = interp(tps)
            #fig, ax = plt.subplots()
            #ax.plot (tr.amp_ps, color = 'red')
            #plt.show()
            #tr.amp_pps = interp(tpps)
            #tr.amp_pss = interp(tpss)

            # Theoretical traces
            #interp = interpolate.interp1d(tpps, tr.amp_ps, bounds_error=False, fill_value=np.nan)
            #tr.amp_pps_theo = interp(tps)
            #interp = interpolate.interp1d(tpss, tr.amp_ps, bounds_error=False, fill_value=np.nan)
            #tr.amp_pss_theo = interp(tps)

            tr.tps = tps
            #tr.tpps = tpps
            #tr.tpss = tpss

        else:
            print("prai: ", tr.prai)
            tr.Xp = -1
            tr.Yp = -1
            tr.Xs = -1
            tr.Ys = -1
            tr.Z = -1
            tr.amp_ps = -1
            #tr.amp_pps = -1
            #tr.amp_pss = -1
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
            try:
                #print(tr.Xs)
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
            except Exception as e:
                print(e)
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
    # l0 = 1
    l0 = 1./111.19
    # dl = 1000000
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

    nbm = 5
    b, a = 3, 1.5
    sigma = 1.0
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
    pal_col = '/home/kalmar/rfmpy/data/colormaps/vik.txt'
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
    n_extra_points = num_of_points  * 2 # number of these points
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
        # 1) add and divide by the number of stacks
        # G = np.divide(amps_matrix_temp, nG)
        # amps.append(G.tolist())

        # 2) Just add (stack) - GO WITH THIS -
        # Whether we stack or add and divide with the number of cells doesn't matter
        # as long as the swaths are the same for all cross-sections (jul 13 2022)
        amps.append(amps_matrix_temp.tolist())
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

    # added this to only keep stations within the swath
    # NOTE Only works for N to S cross sections...
    # TODO: fix this...
    print(temp_lon[0],temp_lon[-1] )
    for index, row in sta.iterrows():
        if row[2] < temp_lon[0] or row[2] > temp_lon[-1]:
            sta = sta.drop(index=index)

    # Grid preparation
    xx = np.arange(0, profile_len, profile_len / n_extra_points)
    zz = np.arange(minz, maxz + pasz, pasz)

    return G2, sta, xx, zz


def create_2d_profile_4_moho_picker(G3, migration_param_dict, profile_points, sta, swath=300, plot=True):
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

    # Is it a N-S or a E-W cross section?
    if lat0 == lat1:
        orientation = 'W-E'
    elif lon0 == lon1:
        orientation = 'S-N'
    else:
        assert False, (
                'OH NO! Cross-section is neither E-W or S-W!'
                'moho picker only supports East to West or South to North orientations for the time being...')

    profile_swath = swath
    # Profile azimuth
    geoid = pyproj.Geod(ellps='WGS84')
    profile_az, back_azimuth, profile_len_ = geoid.inv(lon0, lat0, lon1, lat1)
    # Profile length (km)
    profile_len_gc = profile_len_ / 1000
    if orientation == 'S-N':
        profile_len = degrees2kilometers(lat1 - lat0)
    elif orientation == 'W-E':
        profile_len = degrees2kilometers(lon1 - lon0)

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
        # 1) add and divide by the number of stacks
        # G = np.divide(amps_matrix_temp, nG)
        # amps.append(G.tolist())

        # 2) Just add (stack) - GO WITH THIS -
        # Whether we stack or add and divide with the number of cells doesn't matter
        # as long as the swaths are the same for all cross-sections (jul 13 2022)
        amps.append(amps_matrix_temp.tolist())
    G2 = np.array(amps)  # 2 dimensions
    print("Number of points perpendicular to the profile: ", n_extra_points_, " Swath: ", swath)

    sta_, dxSta, dySta = project_stations(sta=sta, ori_prof=profile_az,
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

    xx = np.arange(0, profile_len, profile_len / n_extra_points)
    xx_gc = np.arange(0, profile_len_gc, profile_len_gc / n_extra_points)
    zz = np.arange(minz, maxz + pasz, pasz)
    # added this to only keep stations within the swath
    if orientation == 'W-E':
        for index, row in sta_.iterrows():
            if row[1] < temp_lat[0] or row[1] > temp_lat[-1]:
               sta_ = sta_.drop(index=index)
        for index, row in sta_.iterrows():
            if row[5] < 0.0 or row[5] > xx[-1]:
                sta_ = sta_.drop(index=index)

    elif orientation == 'S-N':
        for index, row in sta_.iterrows():
            if row[2] < temp_lon[-1] or row[2] > temp_lon[0]:
                sta_ = sta_.drop(index=index)
        for index, row in sta_.iterrows():
            if row[5] < 0.0 or row[5] > xx[-1]:
                sta_ = sta_.drop(index=index)

    return G2, sta_, xx, zz


def plot_migration_profile(Gp, xx, zz, migration_param_dict, sta, work_directory, filename=False, plot_title=None):
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
    pal_col = '/home/kalmar/rfmpy/data/colormaps/vik.txt'
    pal_col = pd.read_csv(pal_col, header=None, index_col=False, sep="\s+", names=["R", "G", "B"])
    cm = LinearSegmentedColormap.from_list("blue2red", pal_col.values, len(pal_col))
    c = np.min([np.max(Gp), 0.1])
    c = 0.15
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
    if plot_title:
        ax.set_title(plot_title)

    plt.tight_layout()
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

# TODO: add docstring
def moho_picker(Gp, xx, zz, migration_param_dict, sta, work_directory, profile, profile_name, path4file):
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

    pal_col = '/home/kalmar/rfmpy/data/colormaps/vik.txt'
    pal_col = pd.read_csv(pal_col, header=None, index_col=False, sep="\s+", names=["R", "G", "B"])
    cm = LinearSegmentedColormap.from_list("blue2red", pal_col.values, len(pal_col))
    c = np.min([np.max(Gp), 0.1])
    c = 0.15
    CL = 2

    plt.close('all')
    # PLOT

    f = plt.figure(1, figsize=[10, 8])
    gs0 = gridspec.GridSpec(nrows=1, ncols=1, figure=f,
                            hspace=0.08, right=0.91, left=0.09, bottom=0.08, top=1, )
    ax = f.add_subplot(gs0[0])  # Ray tracing
    # ax.scatter(XX, ZZ, c=Gp.T, cmap=cm, s=50, vmin=-c / CL, vmax=c / CL, alpha=.5,
    #            zorder=1, picker=True, edgecolors=None, marker='8')
    m = ax.pcolormesh(XX, ZZ, Gp.T, cmap=cm, vmin=-c / CL, vmax=c / CL,
                      zorder=1, shading="auto")
    add_colorbar(ax, m)
    # ax.colorbar(m)
    ax.scatter(sta["XSTA"].values, sta["ZSTA"].values,
               markersize, facecolors="red", edgecolors="k",
               marker="v", lw=0.95, zorder=3, clip_on=False,
               label="Station", )

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
    # ax.set_ylim([80, 0])
    ax.set_xlim([xx[0] - 10 , xx[-1] + 40])
    ax.set_xticks(np.arange(xx[0] , xx[-1], 50))


    ax.tick_params(axis="both", which="major", labelsize=fontsize)
    ax.tick_params(axis="both", which="minor", labelsize=fontsize)

    print("|-----------------------------------------------|\n"
          "|              Moho picker manual               |\n"
          "|-----------------------------------------------|\n"
          "| Make your picks using mouse left button and:  |\n"
          "| - the button m for a certain Moho depth,      |\n"
          "| - the button u for an uncertain pick.         |\n"
          "|-----------------------------------------------|")
    # Is it a N-S or a E-W cross section?
    if profile[0][1] == profile[1][1]:
        orientation = 'W-E'
    elif profile[0][0] == profile[1][0]:
        orientation = 'S-N'
    else:
        assert False, (
                'OH NO! Cross-section is neither E-W or S-W!'
                'moho picker only supports East to West or South to North orientations for the time being...')

    def onkey(event):
        # print(event.key)
        if event.key == 'm':
            if event.xdata is not None and event.ydata is not None:
                print('Dist:', event.xdata, 'Moho:', event.ydata)
                if orientation == 'S-N':
                    lat = profile[0][1] + kilometers2degrees(event.xdata)
                    print('Lon: ', profile[0][0], 'Lat: ', lat, 'Moho:', event.ydata)
                    lon = profile[0][0]
                else:
                    lon = profile[0][0] + kilometers2degrees(event.xdata)
                    print('Lon: ', lon, 'Lat: ', profile[0][1], 'Moho:', event.ydata)
                    lat = profile[0][1]
                # write moho depths
                # this line is for the cross-sections great circle distance
                gc_dist = kilometers2degrees(great_circle(profile[0][0],profile[0][1],lon,lat))
                with open(path4file + '/moho_depths_' + profile_name + '.txt', 'a') as of:
                    of.write('{}, {}, {}, {}\n'.
                             format(lon, lat, event.ydata, gc_dist))
                ax.plot(event.xdata, event.ydata, label='Moho depth',
                        color='black', marker='D',markerfacecolor='white',linestyle='',
                        markersize=7, linewidth=2, alpha=1)
                f.canvas.draw()
        elif event.key == 'u':
            if event.xdata is not None and event.ydata is not None:
                print('Dist:', event.xdata, 'Uncertain Moho:', event.ydata)
                if orientation == 'S-N':
                    lat = profile[0][1] + kilometers2degrees(event.xdata)
                    print('Lon: ', profile[0][0], 'Lat: ', lat, 'Uncertain Moho:', event.ydata)
                    lon = profile[0][0]
                else:
                    lon = profile[0][0] + kilometers2degrees(event.xdata)
                    print('Lon: ', lon, 'Lat: ', profile[0][1], 'Uncertain Moho:', event.ydata)
                    lat = profile[0][1]
                # Write moho depths
                # this line is for the cross-sections great circle distance
                gc_dist = kilometers2degrees(great_circle(profile[0][0],profile[0][1],lon,lat))
                with open(path4file + '/unc_moho_depths_' + profile_name + '.txt', 'a') as of:
                    of.write('{}, {}, {}, {}\n'.
                             format(lon, lat, event.ydata, gc_dist))
                ax.plot(event.xdata, event.ydata, markeredgecolor='black', marker='D',
                        markerfacecolor='gray', linestyle='',markersize=7,
                        linewidth=2, alpha=1, label='Moho')

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
    ax.plot(-20, 10, label='Moho', color='black', marker='D',
            markerfacecolor='white', linestyle='',
            markersize=7, linewidth=2, alpha=1)
    ax.plot(-20, 10, markeredgecolor='black', marker='D',
            markerfacecolor='gray', linestyle='', markersize=7,
            linewidth=2, alpha=1, label='Unc.')
    ax.legend(loc="lower right")


    # f.canvas.mpl_connect('pick_event', onkey)
    cid2 = f.canvas.mpl_connect('key_press_event', onkey)

    plt.show()

    return


def write_files_4_piercing_points_and_raypaths(st, sta, piercing_depth=35, plot=True):
    """..."""

    piercing_lon = []
    piercing_lat = []
    for i, tr in enumerate(st):
        # this try statement gets rid of RFs with no ray paths (ones with tr.Z = -1; tr.Xp = -1 etc)
        try:
            depth = tr.Z[0]
        except:
            continue
        for j, z in enumerate(tr.Z):
            if z > piercing_depth and z < piercing_depth + 1:
                piercing_lon.append(tr.Xp[j])
                piercing_lat.append(tr.Yp[j])
                with open('piercing_points.txt', 'a') as of:
                    of.write('{}, {}\n'.format(tr.Xp[j], tr.Yp[j]))

    # wav_p_lon = []
    # wav_p_lat = []
    # wav_p_dep = []
    # for i, tr in enumerate(st):
    #     try:
    #         depth = tr.Z[0]
    #     except:
    #         continue
    #     for j, z in enumerate(tr.Z):
    #             wav_p_lon.append(tr.Xp[j])
    #             wav_p_lat.append(tr.Yp[j])
    #             wav_p_dep.append(z)
    #             with open('ray_path.txt', 'a') as of:
    #                 of.write('{}, {}, {}\n'.
    #                         format(tr.Xp[j], tr.Yp[j], z))

    if plot:
        # # Plot raypaths
        # plt.scatter(wav_p_lon, wav_p_lat, alpha=0.5,
        #             c=wav_p_dep, marker='.', edgecolor=None, s=1)
        # plt.scatter(sta["LONSTA"], sta["LATSTA"],
        #             c='r', marker='v', edgecolor='k', s=100)
        # plt.show()

        # Plot piercing points
        plt.scatter(piercing_lon, piercing_lat, alpha=.3,
                    c='gray', marker='x', edgecolor='gray', s=50)
        plt.scatter(sta["LONSTA"], sta["LATSTA"],
                    c='r', marker='v', edgecolor='k', s=100)
        plt.show()

    return




