import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from collections.abc import Iterable
from scipy import interpolate
from scipy import signal  # from skimage import measure
from obspy.taup import TauPyModel
import obspy
from rfmpy.visuals import plotting as plt_rf


def project(lat, lon, lato, lono, alpha):

    """Takes station coordinates and projects them
    with respect to the  center of the profile and the angle
    of the profile with respect to the North direction.
    Output is in [km] for x,y coordinates with respect to lono and lato"""

    nbp = len(lat)
    ylat = (lat - lato) * 111.19
    xlon = (lon - lono) * 111.19 * np.cos(np.radians(lat))

    M = np.array(
        [
            [np.cos(np.radians(alpha)), np.sin(np.radians(alpha))],
            [-np.sin(np.radians(alpha)), np.cos(np.radians(alpha))],
        ]
    )

    R = np.dot(np.column_stack((xlon, ylat)), M)

    distx = R[:, 1]
    disty = R[:, 0]

    return distx, disty


def Read_Seismic_Stations(file, ori_prof):

    #################################
    # Read List of Seismic Stations #
    #################################

    sta = pd.read_csv(
        file,
        header=None,
        index_col=False,
        sep="\s+",
        names=["NAMESTA", "LATSTA", "LONSTA", "ALTSTA"],
    )

    """ Project their coordinates with respect
        to the centre of the target linear profile [km] """

    lon_c = sta["LONSTA"].mean()  # Center of the target linear profile
    lat_c = sta["LATSTA"].mean()  # Center of the target linear profile

    xsta, ysta = project(
        sta["LATSTA"].values, sta["LONSTA"].values, lat_c, lon_c, ori_prof
    )

    """ You can set dx, dy for a different coordinate origin
        than the center of the profile """

    dx, dy = 0, 0
    sta["XSTA"] = xsta + dx
    sta["YSTA"] = ysta + dy
    sta["ZSTA"] = -sta["ALTSTA"].values / 1000

    return sta, dx, dy


def Read_Traces(RF_list, sta, ori_prof):

    # Parameters

    is_cor_topo = 1  # 1 -> Includes station elevation information

    # Stations

    xsta = sta["XSTA"].values
    ysta = sta["YSTA"].values

    # Relative shift among stations

    enum = enumerate(sta["NAMESTA"].values)
    dictionary = dict((i, j) for j, i in enum)

    ################################
    # reading R receiver functions #
    ################################

    model = TauPyModel(model="iasp91")
    stream = obspy.Stream()

    with open(RF_list, "r") as file:
        for line in file:

            trace = obspy.read(line.rstrip())[0]

            ista = dictionary[trace.stats.station]
            trace.ista = ista

            trace.kstnm = trace.stats.sac.kstnm
            trace.delta = trace.stats.sac.delta
            trace.x0 = xsta[ista]
            trace.y0 = ysta[ista]
            trace.stla = trace.stats.sac.stla
            trace.stlo = trace.stats.sac.stlo

            if is_cor_topo == 1:
                # trace.z0=-sta['ALTSTA'].values[ista]/1000
                trace.z0 = sta["ZSTA"].values[ista]
            else:
                trace.z0 = 0

            trace.baz = trace.stats.sac.baz
            trace.stats.baz = trace.stats.sac.baz
            trace.lbaz = trace.baz - ori_prof
            trace.gcarc = trace.stats.sac.dist

            if trace.stats.sac.evdp > 800:
                trace.depth = trace.stats.sac.evdp / 1000
            else:
                trace.depth = trace.stats.sac.evdp

            if trace.gcarc < 15:
                trace.prai = model.get_travel_times(
                    source_depth_in_km=trace.depth,
                    distance_in_degree=trace.gcarc,
                    phase_list=["Pn"],
                )[
                    0
                ].ray_param_sec_degree  # SECONDS PER DEGREES

            elif trace.gcarc >= 15 and trace.gcarc <= 100:
                trace.prai = model.get_travel_times(
                    source_depth_in_km=trace.depth,
                    distance_in_degree=trace.gcarc,
                    phase_list=["P"],
                )[
                    0
                ].ray_param_sec_degree  # SECONDS PER DEGREES

            elif trace.gcarc > 100:
                trace.prai = model.get_travel_times(
                    source_depth_in_km=trace.depth,
                    distance_in_degree=trace.gcarc,
                    phase_list=["PKIKP"],
                )[
                    0
                ].ray_param_sec_degree  # SECONDS PER DEGREES

            trace.time = range(len(trace.data)) * trace.delta - trace.stats.sac.a
            trace.filename = line.rstrip()
            trace.rms = np.sqrt(np.mean(np.square(trace.data)))

            stream.append(trace)

    return stream


def tracing_2D(
    stream, ori_prof, path_velocity_model, parameters, lon_c, lat_c, dx=0, dy=0
):

    # Performs ray-tracing necessary for Time-to-depth migration

    stream = stream.copy()

    minx, maxx, pasx = parameters[:3]
    miny, maxy, pasy = parameters[3:6]
    minz, maxz, pasz = parameters[6:9]
    inc, zmax = parameters[9:11]

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

    with open(path_velocity_model + "LonProfile.txt", "r") as f:
        LonProfile = np.array(f.readline().split(), dtype="float")

    with open(path_velocity_model + "zProfile.txt", "r") as f:
        zProfile = np.array(f.readline().split(), dtype="float")

    vProfile = pd.read_csv(
        path_velocity_model + "vProfile.txt", header=None, index_col=None, sep="\s+"
    )

    xProfile, yProfile = project(
        np.ones(LonProfile.shape) * lat_c, LonProfile, lat_c, lon_c, ori_prof
    )

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

            stream[i].Z = Z + stream[i].z0
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


def tracing_1D(stream, ori_prof, parameters, lon_c, lat_c, zMoho=50):

    # Performs ray-tracing necessary for Time-to-depth migration

    stream = stream.copy()

    minx, maxx, pasx = parameters[:3]
    miny, maxy, pasy = parameters[3:6]
    minz, maxz, pasz = parameters[6:9]
    inc, zmax = parameters[9:11]

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

    Z, VP, VS = plt_rf.get_iasp91(zmax=200, step=0.25, zmoho=75)

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
    nbtr = len(stream)

    for i in range(nbtr):
        if stream[i].prai > -1:

            p = stream[i].prai / 111.19
            Xs = stream[i].x0
            Ys = stream[i].y0
            Xp = stream[i].x0
            Yp = stream[i].y0
            coslbaz = np.cos(stream[i].lbaz * np.pi / 180.0)
            sinlbaz = np.sin(stream[i].lbaz * np.pi / 180.0)

            # Migrate with 1-D velocity model
            """ N.B. The backword propagation is computed initially only for 
                standard P- and S-phases.
                The time associated with the propagation of converted-phases
                is computed just after in the subsequent steps """

            incidp = np.arcsin(p * VP)
            incids = np.arcsin(p * VS)

            Ss = np.tan(incids) * inc
            Sp = np.tan(incidp) * inc

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

            stream[i].Z = Z + stream[i].z0
            stream[i].Xp = Xp
            stream[i].Yp = Yp
            stream[i].Xs = Xs
            stream[i].Ys = Ys

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

            interp = interpolate.interp1d(tpps, stream[i].amp_ps)
            stream[i].amp_pps_theo = interp(tps)

            interp = interpolate.interp1d(tpss, stream[i].amp_ps)
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


def ccpM(stream, parameters, sta, phase="PS", stack=0, bazmean=180, dbaz=180):

    # Time to depth Migration

    ##############
    # Parameters #
    ##############

    nbtr = len(stream)

    represent = 0

    distmin = 30
    distmax = 95
    magnmin = -12345
    magnmax = 10
    depthmin = 0
    depthmax = 700

    minx, maxx, pasx = parameters[:3]
    miny, maxy, pasy = parameters[3:6]
    minz, maxz, pasz = parameters[6:9]

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
    for i in range(nbtr):
        baz[i] = stream[i].baz

    ibaz = []
    k = int(np.ceil(2 * dbaz / stack))
    for i in range(k):
        lbazm = bazmean - dbaz + (i + 0.5) * stack
        index1 = np.argwhere(np.abs(baz - lbazm) <= stack / 2)
        index2 = np.argwhere(np.abs(360 - baz + lbazm) <= stack / 2)
        index = np.squeeze(np.concatenate((index1, index2), axis=0))
        ibaz.append(index)

    # rms selection

    rms = np.zeros(len(stream), dtype="float")
    for k in range(len(stream)):
        rms[k] = stream[k].rms
    rms = np.sort(rms)
    i_rms = int(np.floor(len(stream) * 0.98) - 1)
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
            if (
                stream[i].prai > -1
                and stream[i].rms <= rms_max
                and stream[i].gcarc >= distmin
                and stream[i].gcarc <= distmax
                and stream[i].depth >= depthmin
                and stream[i].depth <= depthmax
            ):

                # Look for the correct grid voxel and stack amplitude value

                ikeep[ni] += 1
                ix = np.floor((stream[i].Xs - minx) / pasx)
                iy = np.floor((stream[i].Ys - miny) / pasy)
                iz = np.floor((stream[i].Z - minz) / pasz)
                ix = np.array(ix, dtype="int")
                iy = np.array(iy, dtype="int")
                iz = np.array(iz, dtype="int")

                if phase == "PS":
                    G[ix, iy, iz] = G[ix, iy, iz] + stream[i].amp_ps
                elif phase == "PPS":
                    G[ix, iy, iz] = G[ix, iy, iz] + stream[i].amp_pps[:i1z]
                elif phase == "PSS":
                    G[ix, iy, iz] = G[ix, iy, iz] - stream[i].amp_pss[:i1z]
                elif phase == "PSasPPS":
                    G[ix, iy, iz] = G[ix, iy, iz] + stream[i].amp_pps_theo[:i1z]
                elif phase == "PSasPSS":
                    G[ix, iy, iz] = G[ix, iy, iz] - stream[i].amp_pss_theo[:i1z]
                elif phase == "MU":
                    G[ix, iy, iz] = (
                        G[ix, iy, iz]
                        + stream[i].amp_pps[:i1z]
                        - stream[i].amp_pss[:i1z]
                    )
                elif phase == "PSasMU":
                    G[ix, iy, iz] = (
                        G[ix, iy, iz]
                        + stream[i].amp_pps_theo[:i1z]
                        - stream[i].amp_pss_theo[:i1z]
                    )

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


def ccp_smooth(G2, parameters):

    # Parameters

    # DEPTH SMOOTHING

    minx, maxx, pasx = parameters[:3]
    miny, maxy, pasy = parameters[3:6]
    minz, maxz, pasz = parameters[6:9]
    inc, zmax = parameters[9:11]

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



def rfmops(trace, pref, Z, Vp, Vs):

    #####################################
    # APPLIES NORMAL MOVEOUT CORRECTION #
    #####################################

    # This is useful correction when stacking RFs traces in time
    # It homogenizes the RF arrival times by squeezing and stretching the trace
    # To account for differences in earthquake distance and depth with respect to the same station
    # In fact, give a certain discontinuity below a seismic station
    # The travel time of the coverted P-to-s phase can changed based on the incident angle,
    # which in turn depends on the distance and depth of the earthquake events
    # This correction can help to better stack the traces constructively.

    # It does require a velocity model, which in this case is assumed to be 1D

    trace = trace.copy()

    tbefore = trace.stats.sac.a

    data = trace.data[int(np.floor(tbefore / trace.stats.sac.delta)) :]

    ii = np.arange(1, len(Z) + 1, 1)
    Vp = ii / np.cumsum(1 / Vp)
    Vs = ii / np.cumsum(1 / Vs)

    tps1 = Z * (
        np.sqrt(1 / (Vs ** 2) - trace.prai ** 2)
        - np.sqrt(1 / (Vp ** 2) - trace.prai ** 2)
    )
    Zi = np.interp(np.arange(len(data)) * trace.stats.sac.delta, tps1, Z)
    Vpi = np.interp(Zi, Z, Vp)
    Vsi = np.interp(Zi, Z, Vs)

    t2 = Zi * (
        np.sqrt(1 / (Vsi ** 2) - pref ** 2) - np.sqrt(1 / (Vpi ** 2) - pref ** 2)
    )
    tfinal = np.arange(0, t2[-1], trace.stats.sac.delta)
    data = np.interp(tfinal, t2, data)

    if tbefore > 0:
        data = np.concatenate(
            (trace.data[: int(np.floor(tbefore / trace.stats.sac.delta))], data), axis=0
        )

    if len(data) < len(trace.data):
        newdata = np.zeros(len(trace.data))
        newdata[: len(data)] = data[:]
    else:
        newdata = data[: len(trace.data)]

    return newdata
