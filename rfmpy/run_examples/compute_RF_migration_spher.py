"""
Code for calculating 3D receiver function migrations...
- IN CARTESIAN COORDINATES
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


import rfmpy.core.migration_sphr as rf_mig
import rfmpy.utils.migration_plots as plot_migration_sphr
import numpy as np
import platform
import os
import matplotlib.pyplot as plt


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

plt.scatter(sta["LONSTA"], sta["LATSTA"], c='r')
plt.show()



################
# Read RFs     #
################
stream = rf_mig.read_traces_sphr(path2rfs=path, sta=sta)

# Define MIGRATION parameters
# min lat=46.0
# max lat=48.0
# min lon= 7.0
# max lon = 10.0
# Ray-tracing parameters
inc = 2.5
zmax = 100
# Determine study area (x -> perpendicular to the profile)
minx = 5.0
maxx = 10.0
pasx = 0.5

miny = 45.0
maxy = 50.0
pasy = 0.5

minz = -2
# maxz needs to be >= zmax
maxz = 100
pasz = 5.0
# Pass all the migration parameters in a dictionary to use them in functions called below
m_params = {'minx': minx, 'maxx': maxx, 'pasx': pasx, 'pasy': pasy, 'miny': miny, 'maxy': maxy,
            'minz': minz, 'maxz': maxz, 'pasz': pasz, 'inc': inc, 'zmax': zmax}

################
# Ray tracing  #
################
stream_ray_trace = rf_mig.tracing_3D_sphr(stream=stream, migration_param_dict=m_params, zMoho=50)
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






# Define profile for migration
profile_az = 0
# Define point (lon, lat) to project the station's coordinates with respect
profile_lon = sta["LONSTA"].mean()
profile_lat = sta["LATSTA"].mean()






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

