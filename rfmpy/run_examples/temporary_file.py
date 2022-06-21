"""
Create list of seismic sites used.
Plot seismic network.

Location: Chavannes-pres-renens, CH
Date: Jan 2022
Author: Konstantinos Michailos
"""
import rfmpy.utils.RF_Util as rf_util
import platform

# Set up parameters and paths
if platform.node().startswith('kmichailos-laptop'):
    data_root_dir = '/media/kmichailos/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/home/kmichailos/Desktop/codes/bitbucket'
    desktop_dir = '/home/kmichailos/Desktop'
    hard_drive_dir = '/media/kmichailos/SEISMIC_DATA/'
else:
    data_root_dir = '/media/kmichall/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/home/kmichall/Desktop/Codes/bitbucket'
    desktop_dir = '/home/kmichall/Desktop'
    hard_drive_dir = '/media/kmichall/SEISMIC_DATA/'


path_wavs = [
             # hard_drive_dir + 'RF_data/DATA_RFAA_part_1/SWISS/data/',
             hard_drive_dir + 'RF_data/DATA_RFAA_part_1/EASI/data/',]
             # hard_drive_dir + 'RF_data/DATA_RFAA_part_1/SLOVENIA/data/',
             # hard_drive_dir + 'RF_data/DATA_RFAA_part_2/OBS/data/',
             # hard_drive_dir + 'RF_data/DATA_RFAA_part_1/FRANCE/south_Fr_unsort/',
             # hard_drive_dir + 'RF_data/DATA_RFAA_part_1/FRANCE/data/',
             # hard_drive_dir + 'RF_data/DATA_RFAA_part_1/North_Italy/events_fri_ven/',
             # hard_drive_dir + 'RF_data/DATA_RFAA_part_2/Austria/data_AAA_corrected/',
             # hard_drive_dir + 'RF_data/DATA_RFAA_part_2/CIFAlps/data_YP2012/',
             # hard_drive_dir + 'RF_data/DATA_RFAA_part_2/data_DINAR/',
             # hard_drive_dir + 'RF_data/DATA_RFAA_part_2/HU_SK/data/',
             # hard_drive_dir + 'RF_data/DATA_RFAA_part_3/AARF/DATA_MOBST/data/',
             # hard_drive_dir + 'RF_data/DATA_RFAA_part_3/AARF/DATA_PERMST/data/',
             # hard_drive_dir + 'RF_data/DATA_RFAA_part_3/GERMANY/DE_AA_RF/DATA/data/',
             # hard_drive_dir + 'RF_data/CIFALPS/cifalps_unsort/',
             # hard_drive_dir + 'RF_data/INGV-Permanent-data/',
             # hard_drive_dir + 'RF_data/INGV-Temporary-data/data/',
# hard_drive_dir + 'RF_data/AAPA/data/',]



# List of unique seismic sites
sta1 = rf_util.get_station_info(path_wavs)

sta = sta1

# CALCUlATED RFS
path_wavs_list_part = [hard_drive_dir + 'RF_calculations/RF/']
sta = rf_util.get_station_info(path_wavs_list_part)
unique_all_sta = []
for s in sta:
    if s not in unique_all_sta:
        unique_all_sta.append(s)
# number of RFs on each station
for station in unique_all_sta:
    print(station, sta.count(station))
# using this for making the gmt plot

# For plotting see
# wiggle_bins functions in miscellaneous
path_wavs_list_part = ['/home/kmichall/Downloads/PA-test/data/']





"""
Download metadata information to doublecheck the horizontal components are oriented to
North and East. Will need to check this before I apply the rotations.

Location: Chavannes-pres-renens, CH
Date: Jan 2022
Author: Konstantinos Michailos
"""

from obspy.clients.fdsn import Client
import platform
from obspy import read_inventory, read_events, read

# Set up parameters and paths
if platform.node().startswith('kmichailos-laptop'):
    data_root_dir = '/media/kmichailos/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/home/kmichailos/Desktop/codes/bitbucket'
    desktop_dir = '/home/kmichailos/Desktop'
    hard_drive_dir = '/media/kmichailos/SEISMIC_DATA/'
else:
    data_root_dir = '/media/kmichall/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/home/kmichall/Desktop/Codes/bitbucket'
    desktop_dir = '/home/kmichall/Desktop'
    hard_drive_dir = '/media/kmichall/SEISMIC_DATA/'



azimuths = []
dips = []
with open('Z3.txt', 'r') as f:
    for line in f:
        if line.startswith('#'):
            print(line)
            continue
        else:
            ln = line.split('|')
            az = ln[8]
            dp = ln[9]
            dips.append(dp)
            azimuths.append(az)
            print(ln[1],ln[3], az, dp)


path_wavs = '/media/kmichall/SEISMIC_DATA/RF_data/DATA_RFAA_part_1/SWISS/data/'

from obspy.signal.rotate import rotate2zne

# for each station we need the azimuth and the dip...

stream = read(path_wavs + 'P_2015.002.08.21.55/' + '*ZUR*')
tr1 = stream[0]
tr2 = stream[1]
tr3 = stream[2]

tr_e = tr1.copy()
tr_n = tr2.copy()
tr_z = tr3.copy()
tr_e.data, tr_n.data, tr_z.data = rotate2zne(tr1.data, 30, -90, tr2.data, 90, 3, tr3.data, 92, 3, inverse=False)


#####################################
"""
Reorient part WIP

"""



import rfmpy.utils.RF_Util as rf_util
from rfmpy.utils import signal_processing
from rfmpy.utils.qc import *
from obspy import read_inventory, read_events, UTCDateTime as UTC
import itertools
from pathlib import Path
import glob
import obspy
import numpy as np
import rfmpy.core.RF_Main as RF
import platform
from obspy import read_inventory, read_events, UTCDateTime as UTC
import os

# Define working directory
work_dir = os.getcwd()
try:
    print('>>> Reading inventory...')
    inv = read_inventory(work_dir + '/data/metadata/*.xml')
    print('>>> Read inventory...')
except Exception as e:
    raise type(e)('>>> TYPE cd ... to move to the base directory of the repository!')



# path_wavs = '/media/kmichall/SEISMIC_DATA/RF_data/DATA_RFAA_part_1/SWISS/data/'
path_ev=path_wavs
all_event_dir = glob.glob(path_ev + '*')
event_dir = all_event_dir[0]

ds=30
c1=10
c2=10
c3=1
c4=1
max_frequency=1.0

import rfmpy.utils.RF_Util as rf_util
from rfmpy.utils import signal_processing
from rfmpy.utils import qc
from obspy import read_inventory, read_events, UTCDateTime as UTC
import itertools
from pathlib import Path
import glob
import obspy
import numpy as np
import os


path_ev = '/home/kmichall/Downloads/PA-test/data/'
all_event_dir = glob.glob(path_ev + '*')
for event_dir in all_event_dir:
    print('Calculating RF for event in: ', event_dir)
    wav_files = glob.glob(event_dir + '/*SAC')
    for wav in wav_files:
        st = obspy.read(wav)
        print(st[0].stats)
        st.plot()










# Test how we interpolate velocity in depth
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RegularGridInterpolator.html
import rfmpy.core.migration_sphr as rf_mig
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os

inc = 0.25
zmax = 100
minx = 0.0
maxx = 4.0
pasx = 0.5
miny = 0.0
maxy = 4.0
pasy = 0.5
minz = -2
maxz = 100
pasz = 2

x = np.arange(minx, maxx, pasx)
y = np.arange(miny, maxy, pasy)
z = np.arange(minz, zmax + inc, inc)

VP, VS = rf_mig.get_iasp91(x, y, z, 50)

P_vel_3D_grid = RegularGridInterpolator((x, y, z), VP, method='nearest')
P_vel_3D_grid_ = RegularGridInterpolator((x, y, z), VP, method='linear')

depths = np.linspace(0, 100, 250)
vel = []
for d in depths:
    pts = np.array([1, 1, d])
    vel.append(P_vel_3D_grid(pts))
    print(d, P_vel_3D_grid(pts)[0])

vel_ = []
for d in depths:
    pts = np.array([1, 1, d])
    vel_.append(P_vel_3D_grid_(pts))


# TODO: read EPcrust
# READ velocity model
with open('EPcrust_0_5.txt', 'r') as f:
    for line in f:
        if line.startswith('#'):
            print(line)
            continue
        else:
            ln = line.split()
            lon = float(ln[0])
            lat = float(ln[1])
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


#   LON       LAT        TOPO     THICK_SED   THICK_UPP   THICK_LOW   VP_SED   VP_UPP     VP_LOW    VS_SED    VS_UPP   VS_LOW    RHO_SED   RHO_UPP   RHO_LOW
#  -56.000    89.500    -4.121     2.173       4.750      4.795       2.580     5.000     6.860     1.050     3.011     3.933     2.116     2.535     2.928
x_ = np.arange(minx, maxx, pasx)
y = np.arange(miny, maxy, pasy)
z = np.arange(minz, zmax + inc, inc)

lower_bound = 120
R = 6371  # Earth's radius
x = (R - z) / R
VP_ = np.zeros((x_.size, y.size, z.size))
VS = np.zeros((x_.size, y.size, z.size))
# for j in range(x_.size):
#     print(j)
#     if x[j] >
for i in range(z.size):
    if z[i] < thick_sed:
        VP_[:, :, i] = vp_sed
    elif z[i] < thick_sed + thick_upp:
        VP_[:, :, i] = vp_upp
    elif z[i] < thick_sed + thick_upp + thick_low:
        VP_[:, :, i] = vp_low
    elif z[i] < lower_bound:
        VP_[:, :, i] = 8.78541 - 0.74953 * x[i]
    else:
        VP_[:, :, i] = -1

EPcrust_vel_3D_grid = RegularGridInterpolator((x_, y, z), VP_, method='nearest')

depths = np.linspace(0, 100, 250)
vel_epcrust = []
for d in depths:
    pts = np.array([1, 1, d])
    vel_epcrust.append(EPcrust_vel_3D_grid(pts))




ax1 = plt.subplot2grid((1, 2), (0, 0), colspan=2)
ax1.plot(vel, depths, zorder=2, color='k', linestyle='solid', label='iasp91')
ax1.plot(vel_epcrust, depths, zorder=2, color='k', linestyle='--', label='EpCrust')

# ax1.scatter(vel, depths, facecolor='white', alpha=1,
#             edgecolor='k', linewidth=1., zorder=3)

ax1.set_ylabel('Depth (km)', fontsize=18)
ax1.set_xlabel('Vp (km/s)', fontsize=18)
plt.legend(loc="lower left", markerscale=1., scatterpoints=1, fontsize=14)

plt.gca().invert_yaxis()
plt.show()




#
# longitudes = []
# latitudes = []
# thick_sediments = []
# thick_upper = []
# thick_lower = []
# vp_sediments = []
# vp_upper = []
# vp_lower = []
# vs_sediments = []
# vs_upper = []
# vs_lower = []
# # TODO: ALT read EPcrust and interpolate it
# # READ velocity model
# with open('EPcrust_0_5.txt', 'r') as f:
#     for line in f:
#         if line.startswith('#'):
#             print(line)
#             continue
#         else:
#             ln = line.split()
#             lon = float(ln[0])
#             lat = float(ln[1])
#             topo = float(ln[2])
#             thick_sed = float(ln[3])
#             thick_upp = float(ln[4])
#             thick_low = float(ln[5])
#             vp_sed = float(ln[6])
#             vp_upp = float(ln[7])
#             vp_low = float(ln[8])
#             vs_sed = float(ln[9])
#             vs_upp = float(ln[10])
#             vs_low = float(ln[11])
#             if lon < 5 and lon > 0 and lat > 40 and lat < 45:
#                 longitudes.append(lon)
#                 latitudes.append(lat)
#                 thick_sediments.append(thick_sed)
#                 thick_upper.append(thick_upp)
#                 thick_lower.append(thick_low)
#                 vp_sediments.append(vp_sed)
#                 vp_upper.append(vp_upp)
#                 vp_lower.append(vp_low)
#                 vs_sediments.append(vs_sed)
#                 vs_upper.append(vs_upp)
#                 vs_lower.append(vs_low)
#
# # longitudes, latitudes, vp_sediments, vp_upper, vp_lower, thick_sediments, thick_upper, thick_lower = [list(v) for v in zip(*sorted(zip(longitudes, latitudes, vp_sediments, vp_upper, vp_lower, thick_sediments, thick_upper, thick_lower)))]
#
#
# x_ = np.array(longitudes)
# y = np.array(latitudes)
#
# #
# # x_ = np.arange(min(x_), max(x_) + 0.5, 0.5)
# # y = np.arange(min(y), max(y) + 0.5, 0.5)
# z = np.arange(minz, zmax + inc, inc)
# vp_sediments = np.array(vp_sediments)
# vp_upper = np.array(vp_upper)
# vp_lower = np.array(vp_lower)
# thick_sediments = np.array(thick_sediments)
# thick_upper = np.array(thick_upper)
# thick_lower = np.array(thick_lower)
#
#
# lower_bound = 120
# R = 6371  # Earth's radius
# x = (R - z) / R
# VP_ = np.zeros((x_.size, y.size, z.size))
# for j in range(x_.size):
#     for k in range(y.size):
#         temp_x = x_[j]
#         temp_y = y[k]
#         for i in range(z.size):
#             if z[i] < thick_sediments[j]:
#                 VP_[j, k, i] = vp_sediments[j]
#             elif z[i] < thick_sediments[j] + thick_upper[j]:
#                 VP_[j, k, i] = vp_upper[j]
#             elif z[i] < thick_sediments[j] + thick_upper[j] + thick_lower[j]:
#                 VP_[j, k, i] = vp_lower[j]
#             elif z[i] < lower_bound:
#                 VP_[j, k, i] = 8.78541 - 0.74953 * x[i]
#             else:
#                 VP_[j, k, i] = -1
#
# from scipy.interpolate import LinearNDInterpolator
#
# P_vel = LinearNDInterpolator((x_, y, z), VP_)
#


# TODO: Read EPcrust and create x, y, z , Vp...


longitudes = []
latitudes = []
thick_sediments = []
thick_upper = []
thick_lower = []
vp_sediments = []
vp_upper = []
vp_lower = []
vs_sediments = []
vs_upper = []
vs_lower = []
# TODO: ALT read EPcrust and interpolate it
# READ velocity model
with open('EPcrust_0_5.txt', 'r') as f:
    for line in f:
        if line.startswith('#'):
            print(line)
            continue
        else:
            ln = line.split()
            lon = float(ln[0])
            lat = float(ln[1])
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
            if lon < 5 and lon > 0 and lat > 40 and lat < 45:
                longitudes.append(lon)
                latitudes.append(lat)
                thick_sediments.append(thick_sed)
                thick_upper.append(thick_upp)
                thick_lower.append(thick_low)
                vp_sediments.append(vp_sed)
                vp_upper.append(vp_upp)
                vp_lower.append(vp_low)
                vs_sediments.append(vs_sed)
                vs_upper.append(vs_upp)
                vs_lower.append(vs_low)

x_ = np.array(longitudes)
y = np.array(latitudes)
z = np.arange(0, zmax, inc)
vp_sediments = np.array(vp_sediments)
vp_upper = np.array(vp_upper)
vp_lower = np.array(vp_lower)
thick_sediments = np.array(thick_sediments)
thick_upper = np.array(thick_upper)
thick_lower = np.array(thick_lower)

#   LON       LAT        TOPO     THICK_SED   THICK_UPP   THICK_LOW   VP_SED   VP_UPP     VP_LOW    VS_SED    VS_UPP   VS_LOW    RHO_SED   RHO_UPP   RHO_LOW
#  -56.000    89.500    -4.121     2.173       4.750      4.795       2.580     5.000     6.860     1.050     3.011     3.933     2.116     2.535     2.928
points = []
values = []
xx = []
yy = []
zz = []
vp = []
for i, _ in enumerate(x_):
    #
    z_0 = 0.0
    point0 = [_, y[i], z_0]
    points.append(point0)
    values.append(vp_sediments[i])
    #
    z_1 = thick_sediments[i]
    point1 = [_, y[i], z_1]
    points.append(point1)
    values.append(vp_sediments[i])
    #
    z_2 = thick_sediments[i] + 0.1
    point2 = [_, y[i], z_2]
    points.append(point2)
    values.append(vp_upper[i])
    #
    z_3 = thick_sediments[i] + thick_upper[i]
    point3 = [_, y[i], z_3]
    points.append(point3)
    values.append(vp_upper[i])
    #
    z_4 = thick_sediments[i] + thick_upper[i] + 0.1
    point4 = [_, y[i], z_4]
    points.append(point4)
    values.append(vp_lower[i])
    #
    z_5 = thick_sediments[i] + thick_upper[i] + thick_lower[i]
    point5 = [_, y[i], z_5]
    points.append(point5)
    values.append(vp_lower[i])
    #
    z_6 = thick_sediments[i] + thick_upper[i] + thick_lower[i] + 0.1
    point6 = [_, y[i], z_6]
    points.append(point6)
    values.append(8.1)
    #
    z_7 = 120
    point7 = [_, y[i], z_7]
    points.append(point7)
    values.append(8.1)


points = np.array(points)
values = np.array(values)
# rescale here is important for making the steps sharp (look at the following link:
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.LinearNDInterpolator.html
# EPcrust link: http://eurorem.bo.ingv.it/EPcrust_solar/
linInter = LinearNDInterpolator(points, values, rescale=True)
# request = np.array([[3.5, 41, 10]])
# linInter(request)

depths = np.linspace(0, 100, 250)
vel_epcrust = []
for d in depths:
    pts = np.array([1.5, 42, d])
    vel_epcrust.append(linInter(pts)[0])

ax1 = plt.subplot2grid((1, 2), (0, 0), colspan=2)
ax1.plot(vel_epcrust, depths, zorder=2, color='k', linestyle='--', label='EpCrust')

# ax1.scatter(vel, depths, facecolor='white', alpha=1,
#             edgecolor='k', linewidth=1., zorder=3)

ax1.set_ylabel('Depth (km)', fontsize=18)
ax1.set_xlabel('Vp (km/s)', fontsize=18)
plt.legend(loc="lower left", markerscale=1., scatterpoints=1, fontsize=14)

plt.gca().invert_yaxis()
plt.show()





import rfmpy.core.migration_sphr as rf_mig
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os

def get_epcrust(min_lon=0, max_lon=15, min_lat=40, max_lat=55):
    """

    """
    from scipy.interpolate import LinearNDInterpolator
    # Read x, y, z, etc .txt file of EPcrust velocity model
    longitudes = []
    latitudes = []
    thick_sediments = []
    thick_upper = []
    thick_lower = []
    vp_sediments = []
    vp_upper = []
    vp_lower = []
    vs_sediments = []
    vs_upper = []
    vs_lower = []
    with open('/home/kmichailos/Desktop/codes/github/rfmpy/data/EPcrust/EPcrust_0_5.txt', 'r') as f:
        for line in f:
            if line.startswith('#'):
                print(line)
                continue
            else:
                ln = line.split()
                lon = float(ln[0])
                lat = float(ln[1])
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
                if lon < max_lon and lon > min_lon and lat > min_lat and lat < max_lat:
                    longitudes.append(lon)
                    latitudes.append(lat)
                    thick_sediments.append(thick_sed)
                    thick_upper.append(thick_upp)
                    thick_lower.append(thick_low)
                    vp_sediments.append(vp_sed)
                    vp_upper.append(vp_upp)
                    vp_lower.append(vp_low)
                    vs_sediments.append(vs_sed)
                    vs_upper.append(vs_upp)
                    vs_lower.append(vs_low)

    x_ = np.array(longitudes)
    y = np.array(latitudes)
    # Vp
    vp_sediments = np.array(vp_sediments)
    vp_upper = np.array(vp_upper)
    vp_lower = np.array(vp_lower)
    # Vs
    vs_sediments = np.array(vs_sediments)
    vs_upper = np.array(vs_upper)
    vs_lower = np.array(vs_lower)
    # Thickness
    thick_sediments = np.array(thick_sediments)
    thick_upper = np.array(thick_upper)
    thick_lower = np.array(thick_lower)

    # Define depth profiles for each EPcrust's grid points
    points = []
    p_velocities = []
    s_velocities = []
    for i, _ in enumerate(x_):

        # First point at Earth's surface. #TOdo: add topgrapphy thickness here and to the rest of the layers
        z_0 = 0.0
        point0 = [_, y[i], z_0]
        points.append(point0)
        p_velocities.append(vp_sediments[i])
        s_velocities.append(vs_sediments[i])
        # Second point at the lower limit of the sediments.
        z_1 = thick_sediments[i]
        point1 = [_, y[i], z_1]
        points.append(point1)
        p_velocities.append(vp_sediments[i])
        s_velocities.append(vs_sediments[i])
        # Third point at the lower limit of the sediments with the velocity below.
        z_2 = thick_sediments[i] + 0.1
        point2 = [_, y[i], z_2]
        points.append(point2)
        p_velocities.append(vp_upper[i])
        s_velocities.append(vs_upper[i])
        # Fourth point at the lower part of the upper crust.
        z_3 = thick_sediments[i] + thick_upper[i]
        point3 = [_, y[i], z_3]
        points.append(point3)
        p_velocities.append(vp_upper[i])
        s_velocities.append(vs_upper[i])
        # Fifth point at the lower part of the upper crust...
        z_4 = thick_sediments[i] + thick_upper[i] + 0.1
        point4 = [_, y[i], z_4]
        points.append(point4)
        p_velocities.append(vp_lower[i])
        s_velocities.append(vs_lower[i])
        # Sixth point at the bottom of the crust...
        z_5 = thick_sediments[i] + thick_upper[i] + thick_lower[i]
        point5 = [_, y[i], z_5]
        points.append(point5)
        p_velocities.append(vp_lower[i])
        s_velocities.append(vs_lower[i])
        # Seventh point at the bottom of the crust with mantle's velocity
        z_6 = thick_sediments[i] + thick_upper[i] + thick_lower[i] + 0.1
        point6 = [_, y[i], z_6]
        points.append(point6)
        p_velocities.append(8.1)
        s_velocities.append(6.7)
        # Eighth point at the mantle...
        z_7 = 120
        point7 = [_, y[i], z_7]
        points.append(point7)
        p_velocities.append(8.1)
        s_velocities.append(6.7)

    points = np.array(points)
    values_p = np.array(p_velocities)
    values_s = np.array(s_velocities)
    # rescale here is important for making the steps sharp (look at the following link:
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.LinearNDInterpolator.html
    # EPcrust link: http://eurorem.bo.ingv.it/EPcrust_solar/
    liner_interpolation_of_velocities_p = LinearNDInterpolator(points, values_p, rescale=True)
    liner_interpolation_of_velocities_s = LinearNDInterpolator(points, values_s, rescale=True)

    return liner_interpolation_of_velocities_p, liner_interpolation_of_velocities_s



p_int, s_int = get_epcrust()
# Todo: to make it at 0.1 precision... the vertical array of position for each 100 meters
depths = np.linspace(-4, 100, 200)
Vs = []
Vp = []
for d in depths:
    pts = np.array([11, 45.5, d])
    Vs.append(s_int(pts)[0])
    Vp.append(p_int(pts)[0])

ax1 = plt.subplot2grid((1, 2), (0, 0), colspan=2)
ax1.plot(Vs, depths, zorder=2, color='k', linestyle='--', label='Vs')
ax1.plot(Vp, depths, zorder=2, color='k', linestyle='-', label='Vp')
ax1.set_ylabel('Depth (km)', fontsize=18)
ax1.set_xlabel('Velocity (km/s)', fontsize=18)
plt.legend(loc="lower left", markerscale=1., scatterpoints=1, fontsize=14)
plt.gca().invert_yaxis()
plt.show()