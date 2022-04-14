"""
Code for calculating 3D receiver function migrations...
-
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
import rfmpy.utils.migration_plots_spher as plot_migration_sphr
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

plt.scatter(sta["LONSTA"], sta["LATSTA"], c='r', marker='v')
plt.show()



################
# Read RFs     #
################
stream = rf_mig.read_traces_sphr(path2rfs=path, sta=sta)

# Define MIGRATION parameters
# min lat=45.0
# max lat=50.0
# min lon= 5.0
# max lon = 10.0
# Ray-tracing parameters
inc = 0.25
zmax = 100
# Determine study area (x -> perpendicular to the profile)
minx = 5.0
maxx = 15.0
pasx = 0.5

miny = 45.0
maxy = 55.0
pasy = 0.5

minz = -2
# maxz needs to be >= zmax
maxz = 100
pasz = 2
# Pass all the migration parameters in a dictionary to use them in functions called below
m_params = {'minx': minx, 'maxx': maxx,
            'pasx': pasx, 'pasy': pasy, 'miny': miny, 'maxy': maxy,
            'minz': minz, 'maxz': maxz, 'pasz': pasz, 'inc': inc, 'zmax': zmax}

################
# Ray tracing  #
################
stream_ray_trace = rf_mig.tracing_3D_sphr(stream=stream, migration_param_dict=m_params,
                                          zMoho=50)
# Plot ray tracing...
plot_migration_sphr.plot_ray_tracing(stream_ray_trace)
################
# Migration    #
################
mObs = rf_mig.ccpm_3d(stream_ray_trace, m_params, phase="PS")

################
# Smoothing    #
################
mObs = rf_mig.ccp_smooth(mObs, m_params)
mObs[np.abs(mObs) < np.max(np.abs(mObs)) * 15 / 100] = 0
mObs = rf_mig.ccpFilter(mObs)

################
# Plotting     #
################
plot_migration_sphr.plot_migration_profile(Gp=mObs, migration_param_dict=m_params, sta=sta,
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



import math

def distance_on_unit_sphere(lat1, long1, lat2, long2):

    # Convert latitude and longitude to
    # spherical coordinates in radians.
    degrees_to_radians = math.pi/180.0

    # phi = 90 - latitude
    phi1 = (90.0 - lat1)*degrees_to_radians
    phi2 = (90.0 - lat2)*degrees_to_radians

    # theta = longitude
    theta1 = long1*degrees_to_radians
    theta2 = long2*degrees_to_radians

    # Compute spherical distance from spherical coordinates.

    # For two locations in spherical coordinates
    # (1, theta, phi) and (1, theta', phi')
    # cosine( arc length ) =
    # sin phi sin phi' cos(theta-theta') + cos phi cos phi'
    # distance = rho * arc length

    cos = (np.sin(phi1)*np.sin(phi2)*np.cos(theta1 - theta2) + np.cos(phi1)*np.cos(phi2))
    arc = np.arccos( cos )

    # Remember to multiply arc by the radius of the earth
    # in your favorite set of units to get length.
    return arc

# Misc

                # TODO: figure out how to use the baz here to find the exact location!!!
                # AS gc_dist IS NOT THE POINT...

                # # Local back-azimuth Y-component
                # coslbaz_s = np.cos(baz_s * np.pi / 180.0)
                # coslbaz_p = np.cos(baz_p * np.pi / 180.0)
                # # Local back-azimuth X-component
                # sinlbaz_s = np.sin(baz_s * np.pi / 180.0)
                # sinlbaz_p = np.sin(baz_p * np.pi / 180.0)
                # import math
                #
                # R = 6378.1 #Radius of the Earth
                # brng = 1.57 #Bearing is 90 degrees converted to radians.
                # d = 15 #Distance in km
                #
                # #lat2  52.20444 - the lat result I'm hoping for
                # #lon2  0.36056 - the long result I'm hoping for.
                #
                # lat1 = np.radians(Yp[iz]) #Current lat point converted to radians
                # lon1 = np.radians(Xp[iz]) #Current long point converted to radians
                #
                # lat2 = math.asin( math.sin(lat1)*math.cos(d/R) + math.cos(lat1)*math.sin(d/R)*math.cos(brng))
                #
                # lon2 = lon1 + math.atan2(math.sin(brng)*math.sin(d/R)*math.cos(lat1),
                #              math.cos(d/R)-math.sin(lat1)*math.sin(lat2))

                # lat2 = math.degrees(lat2)
                # lon2 = math.degrees(lon2)

def getEndpoint(lat1,lon1,bearing,d):
    R = 6371                     #Radius of the Earth
    brng = np.radians(bearing) #convert degrees to radians
    lat1 = np.radians(lat1)    #Current lat point converted to radians
    lon1 = np.radians(lon1)    #Current long point converted to radians
    lat2 = np.arcsin(np.sin(lat1)*np.cos(d/R) + np.cos(lat1)*np.sin(d/R)*np.cos(brng))
    lon2 = lon1 + np.arctan2(np.sin(brng)*np.sin(d/R)*np.cos(lat1),np.cos(d/R)-np.sin(lat1)*np.sin(lat2))
    lat2 = np.degrees(lat2)
    lon2 = np.degrees(lon2)
    return lat2,lon2
                # lat_2, lon_2 = getEndpoint(Yp[iz],Xp[iz],baz_p[iz]-180,gc_dist)
                # print(Yp[iz], Xp[iz])
                # print(lat_2, lon_2)
                #Todo: print distance... calc dist

