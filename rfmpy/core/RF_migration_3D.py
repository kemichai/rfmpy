"""
Function for calculating 3D migration of RFs in spherical coordinates for the AlpArray

Note: Based on codes originally written by Matteo Scarponi.

Location: Chavannes-pres-renens, CH
Date: Mar 2022
Author: Konstantinos Michailos
"""

import obspy
import glob
import numpy as np
from rfmpy.utils import migration_utils_3D
import platform
import os
from scipy import interpolate
from scipy import signal  # from skimage import measure
from obspy.taup import TauPyModel
import obspy
import glob
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from scipy.interpolate import RegularGridInterpolator


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


def read_stations(path2rfs, ori_prof):
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

    """ Project their coordinates with respect
        to the centre of the target linear profile [km] """

    lon_c = sta["LONSTA"].mean()  # Center of the target linear profile
    lat_c = sta["LATSTA"].mean()  # Center of the target linear profile

    xsta, ysta = project(sta["LATSTA"].values, sta["LONSTA"].values, lat_c, lon_c, ori_prof)

    """ You can set dx, dy for a different coordinate origin
        than the center of the profile """

    dx, dy = 0, 0
    sta["XSTA"] = xsta + dx
    sta["YSTA"] = ysta + dy
    # Set elevation with negative numbers in km
    sta["ZSTA"] = (-1) * sta["ALTSTA"].values / 1000

    return sta, dx, dy

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

# Define paths
work_dir = os.getcwd()
path = work_dir + "/data/RF/"

ori_prof = 0
sta, dxSta, dySta = read_stations(path2rfs=path, ori_prof=ori_prof)
####################################################################################
####################################################################################
# Read RF files
# stream = migration_utils_3D.Read_Traces(path2rfs=path, sta=sta, ori_prof=ori_prof)
# Stations x and y values
lonsta = sta["LONSTA"].values
latsta = sta["LATSTA"].values
xsta = sta["XSTA"].values
ysta = sta["YSTA"].values
altsta = sta["ZSTA"].values

#
plt.scatter(lonsta, latsta, c='r')
plt.show()

# Assign a unique number to each single separate station (index)
enum = enumerate(sta["NAMESTA"].values)
dictionary = dict((i, j) for j, i in enum)
# Define model (iasp91)
model = TauPyModel(model="iasp91")

stream = obspy.Stream()
all_rfs = glob.glob(path + '*.SAC')
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
    # TODO: what does gcarc stand for???
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
                                            phase_list=["Pn"],)[0].ray_param_sec_degree
    # This is the distance range we use!
    elif trace.gcarc >= 30 and trace.gcarc <= 95:
        trace.prai = model.get_travel_times(source_depth_in_km=trace.depth,
                                            distance_in_degree=trace.gcarc,
                                            phase_list=["P"],)[0].ray_param_sec_degree
    # Should not have distances greater than 95 degrees...
    elif trace.gcarc > 95:
        trace.prai = model.get_travel_times(source_depth_in_km=trace.depth,
                                            distance_in_degree=trace.gcarc,
                                            phase_list=["PKIKP"],)[0].ray_param_sec_degree

    trace.time = range(len(trace.data)) * trace.delta - trace.stats.sac.a
    trace.filename = rf
    trace.rms = np.sqrt(np.mean(np.square(trace.data)))

    stream.append(trace)

####################################################################################
####################################################################################
# Define migration parameters
# Ray-tracing parameters
inc = 0.25
zmax = 100
# Determine study area (x -> perpendicular to the profile)
minx = -200 + dxSta
maxx = 200 + dxSta
pasx = 1
miny = -200 + dySta
maxy = 200 + dySta
pasy = 1
minz = -2
# maxz needs to be >= zmax
maxz = 100
pasz = 0.5
# Pass all the migration parameters in a dictionary to use them in functions
# m_params = {'minx': minx, 'maxx': maxx, 'pasx': pasx, 'pasy': maxy-miny, 'miny': miny, 'maxy': maxy,
#             'minz': minz, 'maxz': maxz, 'pasz': pasz, 'inc': inc, 'zmax': zmax}

####################################################################################
####################################################################################
# Ray tracing
# stream_ray_trace = migration_utils_3D.tracing_1D(tr=stream, ori_prof=ori_prof,
#                                               migration_param_dict=m_params,
#                                               lon_c=lon_c, lat_c=lat_c, zMoho=50,)
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
z = np.arange(0, zmax + inc, inc)

# Project x, y the same way we project the stations (Cartesian coordinates)
# x_proj, y_proj = project(y, x, lat_c, lon_c, ori_prof)
# Generate a 3D grid
xg, yg ,zg = np.meshgrid(x, y, z)
# Define the velocity values on each point of the grid
zMoho=50
VP, VS = get_iasp91(x, y, z, zMoho)
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
print("3-D Ray Tracing")
for i, tr in enumerate(st):
    if tr.prai > -1:
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
        for iz in range(len(z) -1):
            print(iz)

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
            # TODO: figure out why it crushes... Maybe x and y need to be the same size???
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


####################################################################################
####################################################################################
# Migration
# mObs = migration_utils_3D.ccpM(stream_ray_trace, m_params, sta, phase="PS",
#                                stack=0, dbaz=180, bazmean=180)

# Time to depth Migration

##############
# Parameters #
##############
stack=0
dbaz=180
bazmean=180
phase="PS"

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
# ibaz = ibaz_[0]

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

    G = np.zeros((len(xx), len(yy), len(zz)))
    nG = np.zeros((len(xx), len(yy), len(zz))) + 1e-8
    ikeep[ni] = 0

    for i in ibaz[ni]:
        if (st[i].prai > -1 and st[i].rms <= rms_max and st[i].gcarc >= distmin
                            and st[i].gcarc <= distmax and st[i].depth >= depthmin
                            and st[i].depth <= depthmax):
            # Look for the correct grid voxel and stack amplitude value
            ikeep[ni] += 1
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

####################################################################################
####################################################################################
# Pass all the migration parameters in a dictionary to use them in functions
m_params = {'minx': minx, 'maxx': maxx, 'pasx': pasx, 'pasy': maxy-miny, 'miny': miny, 'maxy': maxy,
            'minz': minz, 'maxz': maxz, 'pasz': pasz, 'inc': inc, 'zmax': zmax}
# Smoothing
mObs = migration_utils_3D.ccp_smooth(G2, m_params)
mObs[np.abs(mObs) < np.max(np.abs(mObs)) * 15 / 100] = 0
mObs = migration_utils_3D.ccpFilter(mObs)

####################################################################################
####################################################################################
# PLoting
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from collections import OrderedDict
import matplotlib.patches as patches


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
Gp = mObs

# YY, ZZ = np.meshgrid(yy, zz)
XX, ZZ = np.meshgrid(xx, zz)

# COLOR PALETTE AND COLORMAP
# cork, bam, broc, vik
pal_col = work_dir+ "/data/colormaps/vik.txt"
pal_col = pd.read_csv(pal_col, header=None, index_col=False, sep="\s+", names=["R", "G", "B"])
cm = LinearSegmentedColormap.from_list("blue2red", pal_col.values, len(pal_col))
c = np.min([np.max(Gp), 0.1])
c = 0.06
CL = 2

# PLOT
f = plt.figure(1, figsize=[10, 8])
gs0 = gridspec.GridSpec(nrows=1, ncols=1, figure=f,
                        hspace=0.08, right=0.91,left=0.09, bottom=0.08,top=0.96,)

ax = f.add_subplot(gs0[0])  # Ray tracing
# bwr, seismic, coolwarm, RdBu
m = ax.pcolormesh(XX, ZZ, Gp.T, cmap=cm, vmin=-c / CL, vmax=c / CL,
                  zorder=1, shading="auto")
add_colorbar(ax, m)
ax.scatter(sta["XSTA"].values, sta["ZSTA"].values,
           markersize,facecolors="grey", edgecolors="k",
           marker="v", lw=0.95, zorder=3, clip_on=False,
           label="Seismic stations",)

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
ax.set_yticks(np.arange(10, 140, 10))
ax.set_ylim([100, 0])

ax.tick_params(axis="both", which="major", labelsize=fontsize)
ax.tick_params(axis="both", which="minor", labelsize=fontsize)
plt.show()
plt.close()
# Plotting
# Migration(Gp=mObs, migration_param_dict=m_params, sta=sta,
#                           work_directory=work_dir,
#                           filename=False)
