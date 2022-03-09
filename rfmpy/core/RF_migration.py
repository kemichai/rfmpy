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


def get_iasp91(zmax=200, step=0.25, zmoho=35):
    """
    Retrieves P-wave, S-wave velocities and depths
    from IASPEI91 global velocity model.

    :type zmax: float
    :param zmax: Maximum depth for obtaining velocity values.
    :type step: float
    :param step: Incremental step to increase depth values.
    :type zmoho: int
    :param zmoho: Moho depth in km.

    :rtype: numpy.ndarrays
    :returns: Array of P-wave, S-wave velocities and their depths.
    """
    z = np.arange(0, zmax + step, step)
    RE = 6371  # Earth's radius
    x = (RE - z) / RE
    VP = np.zeros(z.shape)
    VS = np.zeros(z.shape)
    for i in range(len(z)):
        if z[i] < 20:
            VP[i] = 5.8
            VS[i] = 3.36
        elif z[i] < zmoho:
            VP[i] = 6.5
            VS[i] = 3.75
        elif z[i] < 120.0:
            VP[i] = 8.78541 - 0.74953 * x[i]
            VS[i] = 6.706231 - 2.248585 * x[i]
        elif z[i] < 210.0:
            VS[i] = 5.75020 - 1.27420 * x[i]
            VP[i] = 25.41389 - 17.69722 * x[i]
        elif z[i] < 410.0:
            VS[i] = 15.24213 - 11.08552 * x[i]
            VP[i] = 30.78765 - 23.25415 * x[i]
        elif z[i] < 660.0:
            VP[i] = 29.38896 - 21.40656 * x[i]
            VS[i] = 17.70732 - 13.50652 * x[i]
        elif z[i] < 760.0:
            VP[i] = 25.96984 - 16.93412 * x[i]
            VS[i] = 20.76890 - 16.53147 * x[i]
        elif z[i] < 2740.0:
            VP[i] = (25.1486 - 41.1538 * x[i] + 51.9932 * x[i] * x[i] - 26.6083 * x[i] * x[i] * x[i])
            VS[i] = (12.9303 - 21.2590 * x[i] + 27.8988 * x[i] * x[i] - 14.1080 * x[i] * x[i] * x[i])
        elif z[i] < 2889.0:
            VP[i] = 14.49470 - 1.47089 * x[i]
            VS[i] = 8.16616 - 1.58206 * x[i]
        elif z[i] < 5153.9:
            VP[i] = 10.03904 + 3.75665 * x[i] - 13.67046 * x[i] * x[i]
            VS[i] = 1.0e-20
        elif z[i] < 6371.0:
            VP[i] = 11.24094 - 4.09689 * x[i] * x[i]
            VS[i] = 3.56454 - 3.45241 * x[i] * x[i]
        else:
            VP[i] = -1
            VS[i] = -1.0
    return z, VP, VS


# Define paths
path = "/home/kmichall/Desktop/RF_test/RF_Km/RF/"

# Read receiver functions
all_rfs = glob.glob(path + '*.SAC')


stream = obspy.Stream()
# for rf in all_rfs:
trace = obspy.read(rf)
