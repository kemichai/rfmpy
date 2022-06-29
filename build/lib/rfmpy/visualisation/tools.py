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