"""
Code for calculating 3D receiver function migrations...

=============================================
Requirements:
    * obspy
    * scipy
    * pandas
    * ...
=============================================

Note: Based on codes originally written by Matteo Scarponi.

Location: Chavannes-pres-renens, CH
Date: Mar 2022
Author: Konstantinos Michailos
"""


import rfmpy.core.migration as rf_mig
import rfmpy.utils.migration_plots as plot_migration
import numpy as np
import platform
import os


# Set up paths
if platform.node().startswith('kmichailos-laptop'):
    data_root_dir = '/media/kmichailos/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/home/kmichailos/Desktop/codes/github'
    desktop_dir = '/home/kmichailos/Desktop'
    hard_drive_dir = '/media/kmichailos/SEISMIC_DATA/'
else:
    data_root_dir = '/media/kmichall/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/github'
    desktop_dir = '/home/kmichall/Desktop'
    hard_drive_dir = '/media/kmichall/SEISMIC_DATA/'

# Define paths
work_dir = os.getcwd()
path = work_dir + "/data/RF/"

#################
# Read stations #
#################
# Read station coordinates from the rfs (sac files) in a pandas dataframe
sta = rf_mig.read_stations_from_sac(path2rfs=path)

# Define profile for migration
profile_az = 0
# Define point (lon, lat) to project the station's coordinates with respect
profile_lon = sta["LONSTA"].mean()
profile_lat = sta["LATSTA"].mean()

# Update the pandas dataframe with the projected values of the stations
sta, dxSta, dySta = rf_mig.project_stations(sta=sta, ori_prof=profile_az,
                                            point_lat=profile_lat, point_lon=profile_lon)
# TODO: plot stations and profile on a map...
# #
# plt.scatter(lonsta, latsta, c='r')
# plt.show()

################
# Read RFs     #
################
stream = rf_mig.read_traces(path2rfs=path, sta=sta, ori_prof=profile_az)

# Define MIGRATION parameters
# Ray-tracing parameters
inc = 0.25
zmax = 100
# Determine study area (x -> perpendicular to the profile)
minx = -200 + dxSta
maxx = 200 + dxSta
# TODO: step to update in lon lat
pasx = 1
miny = -200 + dySta
maxy = 200 + dySta
# TODO: step to update in lon lat
pasy = 1
minz = -2
# maxz needs to be >= zmax
maxz = 100
pasz = 0.5
# Pass all the migration parameters in a dictionary to use them in functions called below
m_params = {'minx': minx, 'maxx': maxx, 'pasx': pasx, 'pasy': pasy, 'miny': miny, 'maxy': maxy,
            'minz': minz, 'maxz': maxz, 'pasz': pasz, 'inc': inc, 'zmax': zmax}

################
# Ray tracing  #
################
stream_ray_trace = rf_mig.tracing_3D(stream=stream, ori_prof=profile_az,
                                     migration_param_dict=m_params,
                                     point_lon=profile_lon,
                                     point_lat=profile_lat, zMoho=50)
# stream_ray_trace = rf_mig.tracing_1D(tr=stream, ori_prof=ori_prof,
#                                               migration_param_dict=m_params,
#                                               lon_c=lon_c, lat_c=lat_c, zMoho=50,)
# stream_ray_trace = rf_mig.tracing_2D(stream=stream, ori_prof=ori_prof, path_velocity_model=work_dir,
#                           parameters=migration_params, lon_c=lon_c, lat_c=lat_c, dx=dxSta, dy=dySta,)
################
# Migration    #
################
mObs = rf_mig.ccpM(stream_ray_trace, m_params, sta, phase="PS",
                   stack=0, dbaz=180, bazmean=180)

################
# Smoothing    #
################
mObs = rf_mig.ccp_smooth(mObs, m_params)
mObs[np.abs(mObs) < np.max(np.abs(mObs)) * 15 / 100] = 0
mObs = rf_mig.ccpFilter(mObs)

################
# Plotting     #
################
plot_migration.plot_migration_profile(Gp=mObs, migration_param_dict=m_params, sta=sta,
                                      work_directory=work_dir, filename=False)


import numpy as np
# Transformation from cartesian (x, y, z) to spherical coordinates (r, theta, phi)
# a = [x, y, z]
# TODO: make sure to test if it works with negative x,y,z values...
a_cart = [13., 11., 22.]
x = a_cart[0]
y = a_cart[1]
z = a_cart[2]
print(x,y,z)

# cartesian2spherical
r = np.sqrt(x**2 + y**2 + z**2)
if x >= 0:
    phi = np.arctan(np.sqrt(x ** 2 + y ** 2)/z)
elif x < 0:
    phi = np.arctan(np.sqrt(x ** 2 + y ** 2)/z) + np.pi
theta = np.arctan(y / x)

# rad to deg --> * 180/np.pi
# deg to rad --> * np.pi/180

# spherical2cartesian
x_conv = r * np.sin(phi) * np.cos(theta)
y_conv = r * np.sin(phi) * np.sin(theta)
z_conv = r * np.cos(phi)
print(x_conv,y_conv,z_conv)

