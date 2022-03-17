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
inc = 0.250
zmax = 100
# Determine study area (x -> perpendicular to the profile)
minx = -120 + dxSta
maxx = 120 + dxSta
pasx = 0.50
miny = -100 + dySta
maxy = 100 + dySta
pasy = maxy - miny
minz = -2
# maxz needs to be >= zmax
maxz = 100
pasz = 0.50
# Pass all the migration parameters in a dictionary to use them in functions
m_params = {'minx': minx, 'maxx': maxx, 'pasx': pasx, 'pasy': maxy-miny, 'miny': miny, 'maxy': maxy,
            'minz': minz, 'maxz': maxz, 'pasz': pasz, 'inc': inc, 'zmax': zmax}


####################################################################################
####################################################################################
# Ray tracing
# stream_ray_trace = migration_utils_3D.tracing_1D(tr=stream, ori_prof=ori_prof,
#                                               migration_param_dict=m_params,
#                                               lon_c=lon_c, lat_c=lat_c, zMoho=50,)

tr = stream.copy()

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




# Migration
mObs = migration_utils_3D.ccpM(stream_ray_trace, m_params, sta, phase="PS", stack=0, dbaz=180, bazmean=180)
# Smoothing
mObs = migration_utils_3D.ccp_smooth(mObs, m_params)
mObs[np.abs(mObs) < np.max(np.abs(mObs)) * 15 / 100] = 0
mObs = migration_utils_3D.ccpFilter(mObs)
# Plotting
migration_utils_3D.Migration(Gp=mObs, migration_param_dict=m_params, sta=sta,
                          work_directory=work_dir,
                          filename=False)
