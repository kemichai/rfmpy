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

# TODO: look at what this parameters stands for...
ori_prof = 0
prof_azimuth = 90

# TODO: remove the projecting stuff...
sta, dxSta, dySta = migration_utils_3D.read_stations(path2rfs=path, ori_prof=ori_prof)
lon_c = sta["LONSTA"].mean()  # Center of the profile
lat_c = sta["LATSTA"].mean()  # Center of the profile

####################################################################################
# Write the code first without functions and then move them to functions...
####################################################################################
# Read RF files
# stream = migration_utils_3D.Read_Traces(path2rfs=path, sta=sta, ori_prof=ori_prof)

# Stations x and y values
lonsta = sta["LONSTA"].values
latsta = sta["LATSTA"].values

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
    # TODO: see where this is used
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


# Define migration parameters
# Ray-tracing parameters
inc = 5
zmax = 100
# Determine study area (x -> perpendicular to the profile)
minx = 7.0 + dxSta
maxx = 10.0 + dxSta
pasx = 1.0
miny = 45 + dySta
maxy = 48 + dySta
pasy = 1.0
minz = -2
# maxz needs to be >= zmax
maxz = 100
pasz = 10
# Pass all the migration parameters in a dictionary to use them in functions
m_params = {'minx': minx, 'maxx': maxx, 'pasx': pasx, 'pasy': maxy-miny, 'miny': miny, 'maxy': maxy,
            'minz': minz, 'maxz': maxz, 'pasz': pasz, 'inc': inc, 'zmax': zmax}


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

if maxz < zmax:
    print("Problem: maxz < zmax !!")
    quit()

if pasz < inc:
    print("Problem: pasz < inc !!")
    quit()


# 1D velocity model
Z = np.arange(inc, zmax + inc, inc)
zMoho=50
VP, VS = migration_utils_3D.get_iasp91(Z, zMoho)
Z = np.concatenate(([0], Z), axis=0)

print(Z.shape)
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
        p = tr.prai / 111.19
        Xs = tr.lon0
        Ys = tr.lat0
        Xp = tr.lon0
        Yp = tr.lat0
        coslbaz = np.cos(tr.lbaz * np.pi / 180.0)
        sinlbaz = np.sin(tr.lbaz * np.pi / 180.0)

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

        tr.Z = Z + tr.alt
        tr.Xp = Xp
        tr.Yp = Yp
        tr.Xs = Xs
        tr.Ys = Ys

        interp = interpolate.interp1d(tr.time, tr.data, bounds_error=False, fill_value=np.nan)

        tps = -Tp + Ts + Td
        tpps = Tp + Ts + Td - Te
        tpss = 2 * Ts + 2 * Td - Te

        tr.amp_ps = interp(tps)
        tr.amp_pps = interp(tpps)
        tr.amp_pss = interp(tpss)

        # Theoretical traces
        interp = interpolate.interp1d(tpps, tr.amp_ps)
        tr.amp_pps_theo = interp(tps)

        interp = interpolate.interp1d(tpss, tr.amp_ps)
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



# Migration
# mObs = migration_utils_3D.ccpM(stream_ray_trace, m_params, sta, phase="PS", stack=0, dbaz=180, bazmean=180)

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














# Smoothing
mObs = migration_utils_3D.ccp_smooth(mObs, m_params)
mObs[np.abs(mObs) < np.max(np.abs(mObs)) * 15 / 100] = 0
mObs = migration_utils_3D.ccpFilter(mObs)
# Plotting
migration_utils_3D.Migration(Gp=mObs, migration_param_dict=m_params, sta=sta,
                          work_directory=work_dir,
                          filename=False)
