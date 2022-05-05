"""
Functions for calculating 3D migration of RFs in cartesian coordinates.

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


# TODO: finish documentation...
def project(station_lats, station_lons, point_lat, point_lon, angle):
    """
    Projects stations coordinates to a given point (lon, lat) in respect to an angle to the north.

    NOTE: Takes station coordinates and projects them with respect to the center of the profile and the angle
          of the profile with respect to the North direction.
          Output is in [km] for x,y coordinates with respect to lono and lato

    :type station_lats:
    :param station_lats:
    :type station_lons:
    :param station_lons:
    :type point_lat:
    :param point_lat:
    :type point_lon:
    :param point_lon:
    :type angle:
    :param angle:

    :returns:
    """

    ylat = (station_lats - point_lat) * 111.19
    xlon = (station_lons - point_lon) * 111.19 * np.cos(np.radians(station_lats))

    M = np.array([[np.cos(np.radians(angle)), np.sin(np.radians(angle))],
                  [-np.sin(np.radians(angle)), np.cos(np.radians(angle))],])
    R = np.dot(np.column_stack((xlon, ylat)), M)

    distx = R[:, 1]
    disty = R[:, 0]

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
        print(f"| Trace {i} of {len(all_rfs)}")
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

        # Making sure depths are in km (but in a quite weird way...)
        # if trace.stats.sac.evdp > 800:
        #     trace.depth = trace.stats.sac.evdp / 1000
        # else:
        #     trace.depth = trace.stats.sac.evdp
        # depths are in km
        trace.depth = trace.stats.sac.evdp

        # Should not have distances smaller than 30 degrees...
        if trace.gcarc < 30:
            # Seconds per degrees
            trace.prai = model.get_travel_times(source_depth_in_km=trace.depth,
                                                distance_in_degree=trace.gcarc,
                                                phase_list=["Pn"], )[0].ray_param_sec_degree
        # This is the distance range we use!
        elif trace.gcarc >= 30 and trace.gcarc <= 95.1:
            trace.prai = model.get_travel_times(source_depth_in_km=trace.depth,
                                                distance_in_degree=trace.gcarc,
                                                phase_list=["P"], )[0].ray_param_sec_degree
        # Should not have distances greater than 95 degrees...
        elif trace.gcarc > 95.1:
            trace.prai = model.get_travel_times(source_depth_in_km=trace.depth,
                                                distance_in_degree=trace.gcarc,
                                                phase_list=["PKIKP"], )[0].ray_param_sec_degree

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

    :type z
    :param z
    :type step: float
    :param step: Incremental step to increase depth values.
    :type zmoho: int
    :param zmoho: Moho depth in km.

    :rtype: numpy.ndarrays
    :returns: Array of P-wave, S-wave velocities and their depths.
    """

    RE = 6371  # Earth's radius
    x = (RE - z) / RE
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


def tracing_3D_sphr(stream, migration_param_dict, zMoho):
    # TODO: organise this function...

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
    # --------------#
    if maxz < zmax:
        print("Problem: maxz < zmax !!")
        quit()
    if pasz < inc:
        print("Problem: pasz < inc !!")
        quit()

    # --------------#
    # Velocity model
    x = np.arange(minx, maxx, pasx)
    y = np.arange(miny, maxy, pasy)
    z = np.arange(-5, zmax + inc, inc)
    # Define the velocity values on each point of the grid
    # TODO: Create another function that read the epCrust models!!!
    VP, VS = get_iasp91(x, y, z, zMoho)
    # Interpolate
    P_vel_3D_grid = RegularGridInterpolator((x, y, z), VP)
    S_vel_3D_grid = RegularGridInterpolator((x, y, z), VS)

    # Ray tracing
    st = stream.copy()
    st_len = len(st)
    print("|-----------------------------------------------|")
    print("| 3D ray tracing...                             |")
    for i, tr in enumerate(st):
        if tr.prai > -1:
            print('| Trace ' + str(i + 1) + ' of ' + str(st_len))
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
            _, _, baz_p[0] = gps2dist(tr.stats.sac.evla, tr.stats.sac.evlo, tr.lon0, tr.lat0)
            baz_s = np.zeros(len(z))
            _, _, baz_s[0] = gps2dist(tr.stats.sac.evla, tr.stats.sac.evlo, tr.lon0, tr.lat0)
            # Time it takes for waves to travel through each layer
            Tp = np.zeros(len(z))
            Ts = np.zeros(len(z))
            # -------------------------------
            # Migrate with 3-D velocity model
            # -------------------------------
            for iz in range(len(z) - 1):
                # Loop through the z layers moving downwards (Use departing level’s velocity interpolated from 3D model)
                pts = np.array([Xp[iz], Yp[iz], z[iz]])
                VPinterp[iz] = P_vel_3D_grid(pts)
                r_earth = 6371
                # Calculate departing incidence angle of the ray (p = r_earth * sin(incidence_angle) / V)
                id_p = np.arcsin(p * VPinterp[iz])
                id_degrees_p = np.rad2deg(id_p)
                # Calculate great - circle distance travelled delta_i - 1 (delta)
                ia_i_p = np.arcsin((np.sin(id_p)) / (r_earth -  (z[iz+1]+ (-1) * tr.alt)) * (r_earth - (z[iz]+ (-1) * tr.alt)))
                # 180 - this
                ia_i_degrees_p = 180 - np.rad2deg(ia_i_p)
                # Angle ...
                delta_p = 180 - id_degrees_p - ia_i_degrees_p
                # Distance from A to B in km
                gc_dist_p = np.radians(delta_p) * (r_earth - (z[iz]+ (-1) * tr.alt))
                # Location of B
                lat_2, lon_2 = get_end_point(Yp[iz], Xp[iz], baz_p[iz], gc_dist_p )
                Yp[iz + 1] = lat_2
                Xp[iz + 1] = lon_2
                _, _, baz_p[iz + 1] = gps2dist(tr.stats.sac.evla, tr.stats.sac.evlo, Xp[iz + 1], Yp[iz + 1])
                Tp[iz + 1] = Tp[iz] + (inc / np.cos(id_p)) / VPinterp[iz]

                # Same as above for S wave
                VSinterp[iz] = S_vel_3D_grid(pts)
                # Calculate departing incidence angle of the ray (p = r_earth * sin(incidence_angle) / V)
                id_s = np.arcsin(p * VSinterp[iz])
                id_degrees_s = np.rad2deg(id_s)
                # Calculate great - circle distance travelled delta_i - 1 (delta)
                ia_i_s = np.arcsin(
                    (np.sin(id_s)) / (r_earth - (z[iz + 1] + (-1) * tr.alt)) * (r_earth - (z[iz] + (-1) * tr.alt)))
                # 180 - this
                ia_i_degrees_s = 180 - np.rad2deg(ia_i_s)
                # Angle ...
                delta_s = 180 - id_degrees_s - ia_i_degrees_s
                # Distance from A to B in km
                gc_dist_s = np.radians(delta_s) * (r_earth - (z[iz] + (-1) * tr.alt))

                lat_2, lon_2 = get_end_point(Ys[iz], Xs[iz], baz_s[iz], gc_dist_s)
                Ys[iz + 1] = lat_2
                Xs[iz + 1] = lon_2
                _, _, baz_s[iz + 1] = gps2dist(tr.stats.sac.evla, tr.stats.sac.evlo, Xs[iz + 1], Ys[iz + 1])
                Ts[iz + 1] = Ts[iz] + (inc / np.cos(id_s)) / VSinterp[iz]
            # ____________________end of 3D migration_______
            D = np.sqrt(np.square(Xp - Xs) + np.square(Yp - Ys))
            E = np.sqrt(np.square(Xp - Xp[0]) + np.square(Yp - Yp[0]))

            Td = D * p
            Te = 2 * E * p

            tr.Z = z + tr.alt
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



def ccpm_3d(st, migration_param_dict, output_file, phase="PS"):

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
    rms = np.zeros(len(st), dtype="float")
    for k in range(len(st)):
        rms[k] = st[k].rms
    rms = np.sort(rms)
    i_rms = int(np.floor(len(st) * 0.98) - 1)
    rms_max = rms[i_rms]

    # Amplitude matrix
    G = np.zeros((len(x), len(y), len(z)))
    # Number of amplitudes in each cell of the matrix
    nG = np.zeros((len(x), len(y), len(z))) + 1e-8
    # Longitudes and latitudes of the matrix

    for i, tr in enumerate(st):
        if tr.prai >-1 and tr.rms <= rms_max:
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
            # TODO: figure out what that i1z is???
            elif phase == "PPS":
                G[ix, iy, iz] = G[ix, iy, iz] + st[i].amp_pps[:i1z]
            elif phase == "PSS":
                G[ix, iy, iz] = G[ix, iy, iz] - st[i].amp_pss[:i1z]
            # Number of observations in each cell
            nG[ix, iy, iz] = nG[ix, iy, iz] + 1
        else:
            print(f'Removing trace, {tr}, because of high rms-value.')
    # G2 = np.squeeze((np.sum(G, axis=1)))
    # nG2 = np.squeeze((np.sum(nG, axis=1)))
    # G2 = G2 / nG2
    # Get the average number of the amplitudes
    G = np.divide(G, nG)

    # TODO: Find a place to store G as a numpy array!
    np.save(output_file, G, allow_pickle=False)

    return G







# TODO: move this in the plotting part after the G matrix becomes 2 dimensional...
def ccp_smooth(G2, migration_param_dict):

    # Parameters

    # DEPTH SMOOTHING

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

    # l0, dl parameters decide from which depth the smoothing becomes important and how much

    zbegin_lisse = -2
    l0 = 1
    dl = 100

    with np.errstate(divide="warn"):
        G3 = G2
        for iz in range(G2.shape[1]):
            if zz[iz] < zbegin_lisse:
                G3[:, iz] = G2[:, iz]
            else:
                sigmal = (zz[iz] / dl + l0) / pasx
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

    # CONVOLUTION WITH A GAUSSIAN BELL FOR LOCAL SMOOTH

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
