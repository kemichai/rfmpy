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

# TODO: move the codes bellow in functions after figuring out how to perform the 3d migration

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

# TODO: remove the projecting stuff...
ori_prof = 0
prof_azimuth = 90
sta, dxSta, dySta = migration_utils_3D.read_stations(path2rfs=path, ori_prof=ori_prof)
lon_c = sta["LONSTA"].mean()  # Center of the profile
lat_c = sta["LATSTA"].mean()  # Center of the profile

####################################################################################
####################################################################################
# Read RF files
# stream = migration_utils_3D.Read_Traces(path2rfs=path, sta=sta, ori_prof=ori_prof)

# Stations x and y values
lonsta = sta["LONSTA"].values
latsta = sta["LATSTA"].values
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
    # Projected y value
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
inc = 1
zmax = 100
# Determine study area (x -> perpendicular to the profile)
minx = 7.0 + dxSta
maxx = 11.0 + dxSta
pasx = 1
miny = 45 + dySta
maxy = 49 + dySta
pasy = 1
minz = -2
# maxz needs to be >= zmax
maxz = 100
pasz = 1
# Pass all the migration parameters in a dictionary to use them in functions
# m_params = {'minx': minx, 'maxx': maxx, 'pasx': pasx, 'pasy': maxy-miny, 'miny': miny, 'maxy': maxy,
#             'minz': minz, 'maxz': maxz, 'pasz': pasz, 'inc': inc, 'zmax': zmax}

####################################################################################
####################################################################################
# Ray tracing
# stream_ray_trace = migration_utils_3D.tracing_1D(tr=stream, ori_prof=ori_prof,
#                                               migration_param_dict=m_params,
#                                               lon_c=lon_c, lat_c=lat_c, zMoho=50,)

st = stream.copy()

# --------------#
# Main Program #
# --------------#
# Sanity tests...
if maxz < zmax:
    print("Problem: maxz < zmax !!")
    quit()
if pasz < inc:
    print("Problem: pasz < inc !!")
    quit()


# 1D velocity model
x = np.arange(minx, maxx, pasx)
y = np.arange(miny, maxy, pasy)
z = np.arange(inc, zmax + inc, inc)
zMoho=50
VP, VS = get_iasp91(x, y, z, zMoho)
#
# # Check velocity model
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# min_val = VP.min()
# max_val = VP.max()
# n_x, n_y, n_z = VP.shape
# colormap = plt.cm.plasma
# cut = VP[0,:,:]
# Y, Z = np.mgrid[0:n_y, 0:n_z]
# X = np.zeros((n_y, n_z))
# ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=colormap((cut-min_val)/(max_val-min_val)), shade=False)
# cut = VP[:,-1,:]
# X, Z = np.mgrid[0:n_x, 0:n_z]
# Y = n_y - 1 + np.zeros((n_x, n_z))
# ax.plot_surface(X, Y, Z, rstride=10, cstride=10, facecolors=colormap((cut-min_val)/(max_val-min_val)), shade=False)
# cut = VP[:,:,-1]
# X, Y = np.mgrid[0:n_x, 0:n_y]
# Z = z[-1] + np.zeros((n_x, n_y))
# ax.plot_surface(X, Y, Z, rstride=10, cstride=10, facecolors=colormap((cut-min_val)/(max_val-min_val)), shade=False)
# ax.invert_zaxis()
# ax.set_title("Velocity model")
# ax.set_xlabel('X')
# ax.set_ylabel('Y')
# ax.set_zlabel('Z')
# fig.tight_layout()
# plt.show()

#
# Creating dataset
# Z = np.concatenate(([0], Z), axis=0)
# print(Z.shape)
print(VS.shape)
print(VP.shape)


# Ray tracing
# NOTE: This ray-tracing is slightly approximate
#       but it works for the sake of the time-to-depth migration.
#       The idea is to start from the seismic station and propagate the ray backwards!
#       The initial condition to propagate the ray backwards is given by the
#       seismic ray-parameter and the back-azimuth

print("1-D Ray Tracing")
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
        for iz in range(len(y) - 1):
            print(iz)

            # Find neighbouring indices for all directions
            yok = np.argwhere((Yp[iz] < y))
            yok2 = yok[0]
            yok1 = yok2 - 1
            zok = np.argwhere((z[iz] < z))
            zok2 = zok[0]
            zok1 = zok2 - 1
            xok = np.argwhere((Xp[iz] < x))
            xok2 = xok[0]
            xok1 = xok2 - 1

            # Myy = y, ...
            # take the distance-weighted mean of the four neighbours
            d1 = (Yp[iz] - y[yok1]) # distances
            d2 = (y[yok2] - Yp[iz])
            d3 = (z[iz] - z[zok1])
            d4 = (z[zok2] - z[iz])
            d5 = (Xp[iz] - x[xok1])
            d6 = (x[xok2] - Xp[iz])

            # Weights
            w = np.zeros((2, 2, 2))
            w[0, 0, 0] = d1**2 + d3**2 + d5**2
            w[0, 1, 0] = d2**2 + d3**2 + d5**2
            w[1, 0, 0] = d1**2 + d4**2 + d5**2
            w[1, 1, 0] = d2**2 + d4**2 + d5**2
            w[0, 0, 1] = d1**2 + d3**2 + d6**2
            w[0, 1, 1] = d2**2 + d3**2 + d6**2
            w[1, 0, 1] = d1**2 + d4**2 + d6**2
            w[1, 1, 1] = d2**2 + d4**2 + d6**2

            # distances under the root
            w = np.sqrt(w)
            # weight by inverse distance
            w = 1./w

            # if the point falls on the grid
            if d3 == 0:
                w[1,:,:] = 0
            elif d4 == 0:
                w[0,:,:] = 0
            elif d1 == 0:
                w[:,1,:] = 0
            elif d2 == 0:
                w[:,0,:] = 0
            elif d5 == 0:
                w[:,:,1] = 0
            elif d6 == 0:
                w[:,:,0] = 0
            # normalization
            # Why so many sums???
            w=w/sum(sum(sum(w)))

            # Horizontal displacement at this step
            Ss = (w[0, 0, 0] * np.tan(incids[xok1, yok1, zok1]) + w[0, 1, 0] * np.tan(incids[xok1, yok2, zok1]) +
                  w[1, 0, 0] * np.tan(incids[xok2, yok1, zok1]) + w[1, 1, 0] * np.tan(incids[xok2, yok2, zok1]) +
                  w[0, 0, 1] * np.tan(incids[xok1, yok1, zok2]) + w[0, 1, 1] * np.tan(incids[xok1, yok2, zok2]) +
                  w[1, 0, 1] * np.tan(incids[xok2, yok1, zok2]) + w[1, 1, 1] * np.tan(incids[xok2, yok2, zok2])) * inc

            Sp = (w[0, 0, 0] * np.tan(incidp[xok1, yok1, zok1]) + w[0, 1, 0] * np.tan(incidp[xok1, yok2, zok1]) +
                  w[1, 0, 0] * np.tan(incidp[xok2, yok1, zok1]) + w[1, 1, 0] * np.tan(incidp[xok2, yok2, zok1]) +
                  w[0, 0, 1] * np.tan(incidp[xok1, yok1, zok2]) + w[0, 1, 1] * np.tan(incidp[xok1, yok2, zok2]) +
            # TODO: check if w is correct here
                  w[1, 0, 0] * np.tan(incidp[xok2, yok1, zok2]) + w[1, 1, 1] * np.tan(incidp[xok2, yok2, zok2])) * inc

            # Position on the next layer
            Xs[iz+1] = Xs[iz] + sinlbaz * Ss;
            Ys[iz+1] = Ys[iz] + coslbaz * Ss;
            Xp[iz+1] = Xp[iz] + sinlbaz * Sp;
            Yp[iz+1] = Yp[iz] + coslbaz * Sp;






















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
                # TODO: somethings is up here...
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










# Smoothing
mObs = migration_utils_3D.ccp_smooth(mObs, m_params)
mObs[np.abs(mObs) < np.max(np.abs(mObs)) * 15 / 100] = 0
mObs = migration_utils_3D.ccpFilter(mObs)
# Plotting
migration_utils_3D.Migration(Gp=mObs, migration_param_dict=m_params, sta=sta,
                          work_directory=work_dir,
                          filename=False)



