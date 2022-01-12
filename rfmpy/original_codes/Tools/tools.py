import numpy as np


def pvelIASP(z, zmoho=35):

    # vP from IASPEI91 velocity model

    RE = 6371  # Earth's radius

    x = (RE - z) / RE

    VP = np.zeros(z.shape)

    for i in range(len(z)):
        if z[i] < 20:
            VP[i] = 5.8
        elif z[i] < zmoho:
            VP[i] = 6.5
        elif z[i] < 120.0:
            VP[i] = 8.78541 - 0.74953 * x[i]
        elif z[i] < 210.0:
            VP[i] = 25.41389 - 17.69722 * x[i]
        elif z[i] < 410.0:
            VP[i] = 30.78765 - 23.25415 * x[i]
        elif z[i] < 660.0:
            VP[i] = 29.38896 - 21.40656 * x[i]
        elif z[i] < 760.0:
            VP[i] = 25.96984 - 16.93412 * x[i]
        elif z[i] < 2740.0:
            VP[i] = (
                25.1486
                - 41.1538 * x[i]
                + 51.9932 * x[i] * x[i]
                - 26.6083 * x[i] * x[i] * x[i]
            )
        elif z[i] < 2889.0:
            VP[i] = 14.49470 - 1.47089 * x[i]
        elif z[i] < 5153.9:
            VP[i] = 10.03904 + 3.75665 * x[i] - 13.67046 * x[i] * x[i]
        elif z[i] < 6371.0:
            VP[i] = 11.24094 - 4.09689 * x[i] * x[i]
        else:
            VP[i] = -1

    return VP


def svelIASP(z, zmoho=35):

    # vS from IASPEI91 velocity model

    RE = 6371

    x = (RE - z) / RE

    VS = np.zeros(z.shape)

    for i in range(len(z)):

        if z[i] < 20.0:
            VS[i] = 3.36
        elif z[i] < zmoho:
            VS[i] = 3.75
        elif z[i] < 120.0:
            VS[i] = 6.706231 - 2.248585 * x[i]
        elif z[i] < 210.0:
            VS[i] = 5.75020 - 1.27420 * x[i]
        elif z[i] < 410.0:
            VS[i] = 15.24213 - 11.08552 * x[i]
        elif z[i] < 660.0:
            VS[i] = 17.70732 - 13.50652 * x[i]
        elif z[i] < 760.0:
            VS[i] = 20.76890 - 16.53147 * x[i]
        elif z[i] < 2740.0:
            VS[i] = (
                12.9303
                - 21.2590 * x[i]
                + 27.8988 * x[i] * x[i]
                - 14.1080 * x[i] * x[i] * x[i]
            )
        elif z[i] < 2889.0:
            VS[i] = 8.16616 - 1.58206 * x[i]
        elif z[i] < 5153.9:
            VS[i] = 1.0e-20
        elif z[i] < 6371.0:
            VS[i] = 3.56454 - 3.45241 * x[i] * x[i]
        else:
            VS[i] = -1.0

    return VS


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
