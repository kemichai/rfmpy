"""
Functions for calculating 3D migration of RFs in cartesian coordinates.

Note: Based on codes originally written by Matteo Scarponi.

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
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from scipy import signal


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


def project_stations(sta, ori_prof, point_lat, point_lon):
    """
    # TODO: finish documentation...

    :type sta: Pandas DataFrames.
    :param sta: Station details.
    :type ori_prof:
    :param ori_prof:
    :type point_lat:
    :param point_lat:
    :type point_lon:
    :param point_lon:

    :return:
    """

    xsta, ysta = project(sta["LATSTA"].values, sta["LONSTA"].values, point_lat, point_lon, ori_prof)

    dx, dy = 0, 0
    sta["XSTA"] = xsta + dx
    sta["YSTA"] = ysta + dy
    # Set elevation with negative numbers in km
    sta["ZSTA"] = (-1) * sta["ALTSTA"].values / 1000

    return sta, dx, dy

# TODO: finish documentation...
def read_stations_from_sac(path2rfs):
    """
    ...

    :type path2rfs: str
    :param path2rfs: Path to the stored RF SAC files.
    :type ori_prof:
    :param ori_prof:

    :return:
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

    print(sta)

    return sta


def read_traces(path2rfs, sta, ori_prof):
    """
    Read receiver functions.
    # TODO: finish documentation...

    """

    # Stations x and y values
    lonsta = sta["LONSTA"].values
    latsta = sta["LATSTA"].values
    xsta = sta["XSTA"].values
    ysta = sta["YSTA"].values
    altsta = sta["ZSTA"].values

    # Assign a unique number to each single separate station (index)
    enum = enumerate(sta["NAMESTA"].values)
    dictionary = dict((i, j) for j, i in enum)
    # Define model (iasp91)
    model = TauPyModel(model="iasp91")

    stream = obspy.Stream()
    all_rfs = glob.glob(path2rfs + '*.SAC')
    for rf in all_rfs:
        trace_ = obspy.read(rf)
        trace = trace_[0]
        # Define station's index number
        station_index = dictionary[trace.stats.station]
        trace.station_index = station_index
        # Station's name
        trace.kstnm = trace.stats.sac.kstnm
        # Delta value
        trace.delta = trace.stats.sac.delta
        # Projected x value
        trace.lon0 = lonsta[station_index]
        trace.x0 = xsta[station_index]
        # Projected y value
        trace.lat0 = latsta[station_index]
        trace.y0 = ysta[station_index]
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
        trace.lbaz = trace.baz - ori_prof
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
        elif trace.gcarc >= 30 and trace.gcarc <= 95:
            trace.prai = model.get_travel_times(source_depth_in_km=trace.depth,
                                                distance_in_degree=trace.gcarc,
                                                phase_list=["P"], )[0].ray_param_sec_degree
        # Should not have distances greater than 95 degrees...
        elif trace.gcarc > 95:
            trace.prai = model.get_travel_times(source_depth_in_km=trace.depth,
                                                distance_in_degree=trace.gcarc,
                                                phase_list=["PKIKP"], )[0].ray_param_sec_degree

        trace.time = range(len(trace.data)) * trace.delta - trace.stats.sac.a
        trace.filename = rf
        trace.rms = np.sqrt(np.mean(np.square(trace.data)))

        stream.append(trace)

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


def tracing_3D(stream, ori_prof, migration_param_dict, point_lon, point_lat, zMoho):
    # TODO: sort this function...

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
    # TODO: somewhere here I should convert the lon, lat of a given velocity model to x,y
    x = np.arange(minx, maxx, pasx)
    y = np.arange(miny, maxy, pasy)
    z = np.arange(0, zmax + inc, inc)

    # Project x, y the same way we project the stations (Cartesian coordinates)
    # x_proj, y_proj = project(y, x, lat_c, lon_c, ori_prof)
    # Generate a 3D grid
    xg, yg, zg = np.meshgrid(x, y, z)
    # Define the velocity values on each point of the grid
    zMoho = zMoho
    VP, VS = get_iasp91(x, y, z, zMoho)
    # TODO: read epCrust models
    # Interpolate the values
    # For example VP[8.8, 46.2, -2.55] won't work here...
    P_vel_3D_grid = RegularGridInterpolator((x, y, z), VP)
    S_vel_3D_grid = RegularGridInterpolator((x, y, z), VS)
    # Define any location within the grid (in x, y, z and obtain velocity)
    pts = np.array([-101, 46.2, 22.55])
    # P_velocity = P_vel_3D_grid(pts)

    #
    # Creating dataset
    # Z = np.concatenate(([0], z), axis=0)
    # print(Z.shape)
    # print(VS.shape)
    # print(VP.shape)

    # Ray tracing
    # NOTE: This ray-tracing is slightly approximate
    #       but it works for the sake of the time-to-depth migration.
    #       The idea is to start from the seismic station and propagate the ray backwards!
    #       The initial condition to propagate the ray backwards is given by the
    #       seismic ray-parameter and the back-azimuth
    st = stream.copy()
    st_len = len(st)
    print("3-D Ray Tracing")
    for i, tr in enumerate(st):
        if tr.prai > -1:
            print('>>> Calculating migration for trace: ' + str(i) + '/' + str(st_len))
            # ray parameter
            p = tr.prai / 111.19
            # Local back-azimuth Y-component
            coslbaz = np.cos(tr.lbaz * np.pi / 180.0)
            # Local back-azimuth X-component
            sinlbaz = np.sin(tr.lbaz * np.pi / 180.0)

            VPinterp = np.zeros(len(z))
            VSinterp = np.zeros(len(z))

            # S-ray parameter at surface longitude
            Xs = np.zeros(len(z))
            Xs[0] = tr.x0
            # S-ray parameter at surface latitude
            Ys = np.zeros(len(z))
            Ys[0] = tr.y0
            # P-ray parameter at surface longitude
            Xp = np.zeros(len(z))
            Xp[0] = tr.x0
            # P-ray parameter at surface latitude
            Yp = np.zeros(len(z))
            Yp[0] = tr.y0

            Tp = np.zeros(len(z))
            Ts = np.zeros(len(z))

            vpvs = np.ones(VPinterp.shape) * 1.73
            ivpvs = np.argmin(np.abs(z - 10))
            vpvs[ivpvs:] = 1.78
            # -------------------------------
            # Migrate with 3-D velocity model
            # -------------------------------
            # P and S incidence-angle matrix
            incidp = np.arcsin(p * VP)
            incids = np.arcsin(p * VS)

            # for each layer: compute next one???
            # Needs to be z here...
            for iz in range(len(z) - 1):

                # TODO: local baz has to be updated within loop

                ###############################################
                # MATTEO's code...
                pts = np.array([Xp[iz], Yp[iz], z[iz]])
                VPinterp[iz] = P_vel_3D_grid(pts)
                VSinterp[iz] = S_vel_3D_grid(pts)

                incidp = np.arcsin(p * VPinterp[iz])
                incids = np.arcsin(p * VSinterp[iz])
                Ss = np.tan(incids) * inc
                Sp = np.tan(incidp) * inc
                Xs[iz + 1] = Xs[iz] + coslbaz * Ss
                Ys[iz + 1] = Ys[iz] + sinlbaz * Ss
                Xp[iz + 1] = Xp[iz] + coslbaz * Sp
                Yp[iz + 1] = Yp[iz] + sinlbaz * Sp

                Tp[iz + 1] = Tp[iz] + (inc / np.cos(incidp)) / VPinterp[iz]
                Ts[iz + 1] = Ts[iz] + (inc / np.cos(incids)) / VSinterp[iz]
                ###############################################
                # # MATLAB code that crushes... at some point....
                # x and y need to be the same size...
                # # Find neighbouring indices for all directions
                # yok = np.argwhere((Yp[iz] < y))
                # yok2 = yok[0]
                # yok1 = yok2 - 1
                # zok = np.argwhere((z[iz] < z))
                # zok2 = zok[0]
                # zok1 = zok2 - 1
                # xok = np.argwhere((Xp[iz] < x))
                # xok2 = xok[0]
                # xok1 = xok2 - 1
                #
                # # Myy = y, ...
                # # take the distance-weighted mean of the four neighbours
                # d1 = (Yp[iz] - y[yok1]) # distances
                # d2 = (y[yok2] - Yp[iz])
                # d3 = (z[iz] - z[zok1])
                # d4 = (z[zok2] - z[iz])
                # d5 = (Xp[iz] - x[xok1])
                # d6 = (x[xok2] - Xp[iz])
                #
                # # Weights
                # w = np.zeros((2, 2, 2))
                # w[0, 0, 0] = d1**2 + d3**2 + d5**2
                # w[0, 1, 0] = d2**2 + d3**2 + d5**2
                # w[1, 0, 0] = d1**2 + d4**2 + d5**2
                # w[1, 1, 0] = d2**2 + d4**2 + d5**2
                # w[0, 0, 1] = d1**2 + d3**2 + d6**2
                # w[0, 1, 1] = d2**2 + d3**2 + d6**2
                # w[1, 0, 1] = d1**2 + d4**2 + d6**2
                # w[1, 1, 1] = d2**2 + d4**2 + d6**2
                #
                # # distances under the root
                # w = np.sqrt(w)
                # # weight by inverse distance
                # w = 1./w
                #
                # # if the point falls on the grid
                # if d3 == 0:
                #     w[1,:,:] = 0
                # elif d4 == 0:
                #     w[0,:,:] = 0
                # elif d1 == 0:
                #     w[:,1,:] = 0
                # elif d2 == 0:
                #     w[:,0,:] = 0
                # elif d5 == 0:
                #     w[:,:,1] = 0
                # elif d6 == 0:
                #     w[:,:,0] = 0
                # # normalization
                # # Why so many sums???
                # w=w/sum(sum(sum(w)))
                #
                # # Horizontal displacement at this step
                # Ss = (w[0, 0, 0] * np.tan(incids[xok1, yok1, zok1]) + w[0, 1, 0] * np.tan(incids[xok1, yok2, zok1]) +
                #       w[1, 0, 0] * np.tan(incids[xok2, yok1, zok1]) + w[1, 1, 0] * np.tan(incids[xok2, yok2, zok1]) +
                #       w[0, 0, 1] * np.tan(incids[xok1, yok1, zok2]) + w[0, 1, 1] * np.tan(incids[xok1, yok2, zok2]) +
                #       w[1, 0, 1] * np.tan(incids[xok2, yok1, zok2]) + w[1, 1, 1] * np.tan(incids[xok2, yok2, zok2])) * inc
                #
                # Sp = (w[0, 0, 0] * np.tan(incidp[xok1, yok1, zok1]) + w[0, 1, 0] * np.tan(incidp[xok1, yok2, zok1]) +
                #       w[1, 0, 0] * np.tan(incidp[xok2, yok1, zok1]) + w[1, 1, 0] * np.tan(incidp[xok2, yok2, zok1]) +
                #       w[0, 0, 1] * np.tan(incidp[xok1, yok1, zok2]) + w[0, 1, 1] * np.tan(incidp[xok1, yok2, zok2]) +
                #       w[1, 0, 0] * np.tan(incidp[xok2, yok1, zok2]) + w[1, 1, 1] * np.tan(incidp[xok2, yok2, zok2])) * inc
                #
                # # Position on the next layer
                # Xs[iz+1] = Xs[iz] + sinlbaz * Ss
                # Ys[iz+1] = Ys[iz] + coslbaz * Ss
                # Xp[iz+1] = Xp[iz] + sinlbaz * Sp
                # Yp[iz+1] = Yp[iz] + coslbaz * Sp
                #
                # # following the above scheme...:
                # cosincS =  w[0, 0, 0] * np.cos(incids[xok1,yok1,zok1]) + w[0, 1, 0] * np.cos(incids[xok1,yok2,zok1]) +\
                #            w[1, 0, 0] * np.cos(incids[xok2,yok1,zok1]) + w[1, 1, 0] * np.cos(incids[xok2,yok2,zok1]) +\
                #            w[0, 0, 1] * np.cos(incids[xok1,yok1,zok2]) + w[0, 1, 1] * np.cos(incids[xok1,yok2,zok2]) +\
                #            w[1, 0, 1] * np.cos(incids[xok2,yok1,zok2]) + w[1, 1, 1] * np.cos(incids[xok2,yok2,zok2])
                #
                # cosincP = w[0, 0, 0] * np.cos(incidp[xok1,yok1,zok1]) + w[0, 1, 0] * np.cos(incidp[xok1,yok2,zok1]) +\
                #           w[1, 0, 0] * np.cos(incidp[xok2,yok1,zok1]) + w[1, 1, 0] * np.cos(incidp[xok2,yok2,zok1]) +\
                #           w[0, 0, 1] * np.cos(incidp[xok1,yok1,zok2]) + w[0, 1, 1] * np.cos(incidp[xok1,yok2,zok2]) +\
                #           w[1, 0, 1] * np.cos(incidp[xok2,yok1,zok2]) + w[1, 1, 1] * np.cos(incidp[xok2,yok2,zok2])
                #
                # VSloc = w[0, 0, 0] * VS[xok1,yok1,zok1] + w[0, 1, 0]  * VS[xok1,yok2,zok1] +\
                #         w[1, 0, 0] * VS[xok2,yok1,zok1] + w[1, 1, 0]  * VS[xok2,yok2,zok1] +\
                #         w[0, 0, 1] * VS[xok1,yok1,zok2] + w[0, 1, 1]  * VS[xok1,yok2,zok2] +\
                #         w[1, 0, 1] * VS[xok2,yok1,zok2] + w[1, 1, 1]  * VS[xok2,yok2,zok2]
                # VPloc = w[0, 0, 0] * VP[xok1,yok1,zok1] + w[0, 1, 0]  * VP[xok1,yok2,zok1] +\
                #         w[1, 0, 0] * VP[xok2,yok1,zok1] + w[1, 1, 0]  * VP[xok2,yok2,zok1] +\
                #         w[0, 0, 1] * VP[xok1,yok1,zok2] + w[0, 1, 1]  * VP[xok1,yok2,zok2] +\
                #         w[1, 0, 1] * VP[xok2,yok1,zok2] + w[1, 1, 1]  * VP[xok2,yok2,zok2]
                # # traveltimes
                # Ts[iz + 1] = Ts[iz] + (inc / cosincS) / VSloc
                # Tp[iz + 1] = Tp[iz] + (inc / cosincP) / VPloc

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
            # ???
            tr.Ts = Ts
            tr.Tp = Tp

            from scipy.interpolate import interp1d

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

    return st


def ccpM(st, migration_param_dict, sta, phase="PS", stack=0, bazmean=180, dbaz=180):

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
    inc = migration_param_dict['inc']
    zmax = migration_param_dict['zmax']

    ##############
    # Parameters #
    ##############
    nbtr = len(st)
    represent = 0
    distmin = 30
    distmax = 95
    magnmin = -12345
    magnmax = 10
    depthmin = 0
    depthmax = 700

    # Back-azimuth stacking preparation
    if stack > 0:
        print("Stacking interval: ", str(stack))
    if stack < 0:
        print("* WRONG STACKING INTERVAL * \nStack set to 0")
        stack = 0
    if stack == 0:
        stack = dbaz * 2  # Only one interval is made

    # Traces selection based on back-azimuthal interval
    baz = np.zeros(nbtr, dtype="float")
    for i, tr in enumerate(st):
        baz[i] = tr.baz

    ibaz = []
    k = int(np.ceil(2 * dbaz / stack))
    for i in range(k):
        lbazm = bazmean - dbaz + (i + 0.5) * stack
        index1 = np.argwhere(np.abs(baz - lbazm) <= stack / 2)
        index2 = np.argwhere(np.abs(360 - baz + lbazm) <= stack / 2)
        index = np.squeeze(np.concatenate((index1, index2), axis=0))
        ibaz.append(index)

    # rms selection
    rms = np.zeros(len(st), dtype="float")
    for k in range(len(st)):
        rms[k] = st[k].rms
    rms = np.sort(rms)
    i_rms = int(np.floor(len(st) * 0.98) - 1)
    rms_max = rms[i_rms]

    # Grid preparation
    xx = np.arange(minx, maxx + pasx, pasx)
    yy = np.arange(miny, maxy + pasy, pasy)
    zz = np.arange(minz, maxz + pasz, pasz)
    XX, YY, ZZ = np.meshgrid(xx, yy, zz)

    # XX, YY, ZZ = xg, yg, zg

    # Looping on all the selected traces
    index = {y: x for x, y in enumerate(sta["NAMESTA"].values)}
    G2tmp = []
    nG2 = []
    ikeep = np.zeros(len(ibaz))
    for ni in range(len(ibaz)):
        # amplitudes
        G = np.zeros((len(xx), len(yy), len(zz)))
        # number of aamplitudes in each cell
        nG = np.zeros((len(xx), len(yy), len(zz))) + 1e-8
        ikeep[ni] = 0

        for i in ibaz[ni]:
            if (st[i].prai > -1 and st[i].rms <= rms_max and st[i].gcarc >= distmin
                and st[i].gcarc <= distmax and st[i].depth >= depthmin
                and st[i].depth <= depthmax):
                # Look for the correct grid voxel and stack amplitude value
                ikeep[ni] += 1
                # TODO: change to lon, lat, dep instead of xy
                ix = np.floor((st[i].Xs - minx) / pasx)
                iy = np.floor((st[i].Ys - miny) / pasy)
                iz = np.floor((st[i].Z - minz) / pasz)
                ix = np.array(ix, dtype="int")
                iy = np.array(iy, dtype="int")
                iz = np.array(iz, dtype="int")
                if phase == "PS":
                    G[ix, iy, iz] = G[ix, iy, iz] + st[i].amp_ps

                elif phase == "PPS":
                    G[ix, iy, iz] = G[ix, iy, iz] + st[i].amp_pps[:i1z]
                elif phase == "PSS":
                    G[ix, iy, iz] = G[ix, iy, iz] - st[i].amp_pss[:i1z]
                elif phase == "PSasPPS":
                    G[ix, iy, iz] = G[ix, iy, iz] + st[i].amp_pps_theo[:i1z]
                elif phase == "PSasPSS":
                    G[ix, iy, iz] = G[ix, iy, iz] - st[i].amp_pss_theo[:i1z]
                elif phase == "MU":
                    G[ix, iy, iz] = (G[ix, iy, iz] + st[i].amp_pps[:i1z] - st[i].amp_pss[:i1z])
                elif phase == "PSasMU":
                    G[ix, iy, iz] = (G[ix, iy, iz] + st[i].amp_pps_theo[:i1z] - st[i].amp_pss_theo[:i1z])
                nG[ix, iy, iz] = nG[ix, iy, iz] + 1

        # 2D transformation
        G2tmp.append(np.sum(G, axis=1))
        nG2.append(np.sum(nG, axis=1))
        G2tmp[ni] = G2tmp[ni] / nG2[ni]
        i1 = np.argwhere(nG2[ni] > 0)
        nG2[ni][i1] = 1

    # Stack baz-stacks
    G2 = np.zeros(G2tmp[0].shape)
    nG2all = np.zeros(nG2[0].shape)
    for ni in range(len(ibaz)):
        if ikeep[ni] != 0:
            G2 += G2tmp[ni]
            nG2all += nG2[ni]

    G2 = G2 / nG2all

    return G2


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


def tracing_2D(stream, ori_prof, path_velocity_model, migration_param_dict, lon_c, lat_c, dx=0, dy=0):

    # Performs ray-tracing necessary for Time-to-depth migration

    stream = stream.copy()

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
    # Main Program #
    # --------------#

    if maxz < zmax:
        print("Problem: maxz < zmax !!")
        quit()

    if pasz < inc:
        print("Problem: pasz < inc !!")
        quit()

    Z = np.arange(0, zmax + inc, inc)

    ###########################
    # Read 2-D Velocity Model #
    ###########################

    """ Here I read a 2D velocity model stored in my computer
        for the extension of the target profile along which I want to 
        perform the migration.
        In this case, I define an interpolator:
        f = interpolate.interp2d(...)
        so that I can extract the velocity values at the location I need during the 
        ray tracing """

    with open(path_velocity_model + "/LonProfile.txt", "r") as f:
        LonProfile = np.array(f.readline().split(), dtype="float")

    with open(path_velocity_model + "/zProfile.txt", "r") as f:
        zProfile = np.array(f.readline().split(), dtype="float")

    vProfile = pd.read_csv(path_velocity_model + "/vProfile.txt", header=None, index_col=None, sep="\s+")

    xProfile, yProfile = project(np.ones(LonProfile.shape) * lat_c, LonProfile, lat_c, lon_c, ori_prof)

    xProfile += dx
    yProfile += dy

    xGrid, zGrid = np.meshgrid(xProfile, zProfile)

    """ To define the interpolator I expect the x,z,v to have a 1D,1D,2D shape
        such as, for example:

        xProfile.shape == (466,)
        zProfile.shape == (201,)
        vProfile.shape == (201,466)

        Interpolator defined here: """

    f = interpolate.interp2d(xProfile, zProfile, vProfile, kind="linear")

    ###############
    # RAY TRACING #
    ###############

    """ This ray-tracing is slightly approximate 
        but it works for the sake of the time-to-depth migration.
        The idea is to start from the seismic station and propagate the ray backwards!
        The initial condition to propagate the ray backwards is given by the
        seismic ray-parameter and the back-azimuth """

    nbtr = len(stream)
    for i in range(nbtr):
        if stream[i].prai > -1:

            p = stream[i].prai / 111.19
            coslbaz = np.cos(stream[i].lbaz * np.pi / 180.0)
            sinlbaz = np.sin(stream[i].lbaz * np.pi / 180.0)

            VPinterp = np.zeros(len(Z))
            VSinterp = np.zeros(len(Z))

            Xs = np.zeros(len(Z))
            Ys = np.zeros(len(Z))
            Xp = np.zeros(len(Z))
            Yp = np.zeros(len(Z))

            Xs[0] = stream[i].x0
            Ys[0] = stream[i].y0
            Xp[0] = stream[i].x0
            Yp[0] = stream[i].y0

            Tp = np.zeros(len(Z))
            Ts = np.zeros(len(Z))

            vpvs = np.ones(VPinterp.shape) * 1.73
            ivpvs = np.argmin(np.abs(Z - 10))
            vpvs[ivpvs:] = 1.78

            # Start of the 2D migration for the given event-station
            """ N.B. The backword propagation is computed initially only for 
                standard P- and S-phases.
                The time associated with the propagation of converted-phases
                is computed just after in the subsequent steps """

            for j in range(len(Z) - 1):

                VPinterp[j] = f(Xp[j], Z[j])
                VSinterp[j] = VPinterp[j] / vpvs[j]

                incidp = np.arcsin(p * VPinterp[j])
                incids = np.arcsin(p * VSinterp[j])

                Ss = np.tan(incids) * inc
                Sp = np.tan(incidp) * inc

                Xs[j + 1] = Xs[j] + coslbaz * Ss
                Ys[j + 1] = Ys[j] + sinlbaz * Ss
                Xp[j + 1] = Xp[j] + coslbaz * Sp
                Yp[j + 1] = Yp[j] + sinlbaz * Sp

                Tp[j + 1] = Tp[j] + (inc / np.cos(incidp)) / VPinterp[j]
                Ts[j + 1] = Ts[j] + (inc / np.cos(incids)) / VSinterp[j]

            # End of 2D Migration
            """ Once that ray geometry is provided
                Propagation time is computed for all the converthed phases.
                This will be useful for the next stages of time to depth migration """

            D = np.sqrt(np.square(Xp - Xs) + np.square(Yp - Ys))
            E = np.sqrt(np.square(Xp - Xp[0]) + np.square(Yp - Yp[0]))

            Td = D * p
            Te = 2 * E * p

            stream[i].Z = Z + stream[i].alt
            stream[i].Xp = Xp
            stream[i].Yp = Yp
            stream[i].Xs = Xs
            stream[i].Ys = Ys
            stream[i].Ts = Ts
            stream[i].Tp = Tp

            interp = interpolate.interp1d(
                stream[i].time, stream[i].data, bounds_error=False, fill_value=np.nan
            )

            tps = -Tp + Ts + Td
            tpps = Tp + Ts + Td - Te
            tpss = 2 * Ts + 2 * Td - Te

            stream[i].amp_ps = interp(tps)
            stream[i].amp_pps = interp(tpps)
            stream[i].amp_pss = interp(tpss)

            # Theoretical traces

            interp = interpolate.interp1d(
                tpps, stream[i].amp_ps, bounds_error=False, fill_value=np.nan
            )

            stream[i].amp_pps_theo = interp(tps)

            interp = interpolate.interp1d(
                tpss, stream[i].amp_ps, bounds_error=False, fill_value=np.nan
            )

            stream[i].amp_pss_theo = interp(tps)

            stream[i].tps = tps
            stream[i].tpps = tpps
            stream[i].tpss = tpss

        else:
            print("prai: ", stream[i].prai)
            stream[i].Xp = -1
            stream[i].Yp = -1
            stream[i].Xs = -1
            stream[i].Ys = -1
            stream[i].Z = -1
            stream[i].amp_ps = -1
            stream[i].amp_pps = -1
            stream[i].amp_pss = -1

    return stream


def tracing_1D(tr, ori_prof, migration_param_dict, lon_c, lat_c, zMoho=50):
    """

    """

    # Performs ray-tracing necessary for Time-to-depth migration

    tr = tr.copy()

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
    # Main Program #
    # --------------#

    if maxz < zmax:
        print("Problem: maxz < zmax !!")
        quit()

    if pasz < inc:
        print("Problem: pasz < inc !!")
        quit()

    # 1D velocity model

    Z = np.arange(inc, zmax + inc, inc)

    VP, VS = get_iasp91(Z, zMoho)
    Z = np.concatenate(([0], Z), axis=0)

    print(Z.shape)
    print(VS.shape)
    print(VP.shape)

    ###############
    # RAY TRACING #
    ###############

    """ This ray-tracing is slightly approximate 
        but it works for the sake of the time-to-depth migration.
        The idea is to start from the seismic station and propagate the ray backwards!
        The initial condition to propagate the ray backwards is given by the
        seismic ray-parameter and the back-azimuth """

    print("1-D Ray Tracing")
    nbtr = len(tr)

    for i in range(nbtr):
        if tr[i].prai > -1:

            p = tr[i].prai / 111.19
            Xs = tr[i].x0
            Ys = tr[i].y0
            Xp = tr[i].x0
            Yp = tr[i].y0
            coslbaz = np.cos(tr[i].lbaz * np.pi / 180.0)
            sinlbaz = np.sin(tr[i].lbaz * np.pi / 180.0)

            # Migrate with 1-D velocity model
            """ N.B. The backword propagation is computed initially only for 
                standard P- and S-phases.
                The time associated with the propagation of converted-phases
                is computed just after in the subsequent steps """

            # P and S incidence-angle matrix
            incidp = np.arcsin(p * VP)
            incids = np.arcsin(p * VS)
            # horizontal displacement
            Ss = np.tan(incids) * inc
            Sp = np.tan(incidp) * inc

            # position on the next layer
            Xs = np.concatenate(([Xs], coslbaz * Ss), axis=0)
            Ys = np.concatenate(([Ys], sinlbaz * Ss), axis=0)
            Xp = np.concatenate(([Xp], coslbaz * Sp), axis=0)
            Yp = np.concatenate(([Yp], sinlbaz * Sp), axis=0)

            Xs = np.cumsum(Xs)
            Ys = np.cumsum(Ys)
            Xp = np.cumsum(Xp)
            Yp = np.cumsum(Yp)

            if not Xs.any():
                print("!!! All zero tracing")
            if not Z.any():
                print("Problem")
                quit()

            Tp = np.concatenate(([0], (inc / np.cos(incidp)) / VP))
            Ts = np.concatenate(([0], (inc / np.cos(incids)) / VS))
            Tp = np.cumsum(Tp)
            Ts = np.cumsum(Ts)

            # End of 1D Migration
            """ Once that ray geometry is provided
                Propagation time is computed for all the converthed phases.
                This will be useful for the next stages of time to depth migration """

            D = np.sqrt(np.square(Xp - Xs) + np.square(Yp - Ys))
            E = np.sqrt(np.square(Xp - Xp[0]) + np.square(Yp - Yp[0]))

            Td = D * p
            Te = 2 * E * p

            tr[i].Z = Z + tr[i].alt
            tr[i].Xp = Xp
            tr[i].Yp = Yp
            tr[i].Xs = Xs
            tr[i].Ys = Ys

            interp = interpolate.interp1d(tr[i].time, tr[i].data, bounds_error=False, fill_value=np.nan)

            tps = -Tp + Ts + Td
            tpps = Tp + Ts + Td - Te
            tpss = 2 * Ts + 2 * Td - Te

            tr[i].amp_ps = interp(tps)
            tr[i].amp_pps = interp(tpps)
            tr[i].amp_pss = interp(tpss)

            # Theoretical traces
            interp = interpolate.interp1d(tpps, tr[i].amp_ps)
            tr[i].amp_pps_theo = interp(tps)

            interp = interpolate.interp1d(tpss, tr[i].amp_ps)
            tr[i].amp_pss_theo = interp(tps)

            tr[i].tps = tps
            tr[i].tpps = tpps
            tr[i].tpss = tpss

        else:
            print("prai: ", tr[i].prai)
            tr[i].Xp = -1
            tr[i].Yp = -1
            tr[i].Xs = -1
            tr[i].Ys = -1
            tr[i].Z = -1
            tr[i].amp_ps = -1
            tr[i].amp_pps = -1
            tr[i].amp_pss = -1

    return tr









