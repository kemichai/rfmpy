"""
Location: Chavannes-pres-renens, CH
Date: April 2022
Author: Konstantinos Michailos
"""

import rfmpy.core.migration_sphr as rf_mig
import rfmpy.utils.migration_plots_spher as plot_migration_sphr
import numpy as np
import platform
import os
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
from pyproj import Geod
import pyproj


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



def get_end_point(lat1, lon1, baz, d):
    """
    Calculates the end point in lon, lat given we know:
     1) the initial point,
     2) the distance and
     3) the back azimuth
    """
    # Radius of the Earth
    R = 6371
    # Convert degrees to radians
    baz_ = np.radians(baz)
    # Current lat point converted to radians
    lat1 = np.radians(lat1)
    # Current long point converted to radians
    lon1 = np.radians(lon1)

    lat2 = np.arcsin(np.sin(lat1) * np.cos(d / R) + np.cos(lat1) * np.sin(d / R) * np.cos(baz_))
    lon2 = lon1 + np.arctan2(np.sin(baz_) * np.sin(d / R) * np.cos(lat1),
                             np.cos(d / R) - np.sin(lat1) * np.sin(lat2))
    # Convert to degrees
    lat2 = np.degrees(lat2)
    lon2 = np.degrees(lon2)
    return lat2, lon2

def az2baz(angle):
    if angle >= 180.:
        baz = angle - 180.
    else:
        baz = angle + 180.
    return baz

def get_perpendicular_azimuth(az):
    if az > 270:
        perp_az1 = (az + 90) - 360
        perp_az2 = az - 90
    elif az < 90:
        perp_az1 = az + 90
        perp_az2 = az - 90 + 360
    else:
        perp_az1 = az + 90
        perp_az2 = az - 90
    return perp_az1, perp_az2



# Define paths
work_dir = os.getcwd()
path = work_dir + "/data/RF/"

# Define MIGRATION parameters
# Ray-tracing parameters
inc = 0.25
zmax = 100
# Determine study area (x -> perpendicular to the profile)
minx = 5.0
maxx = 12.0
pasx = 0.5
miny = 45.0
maxy = 52.0
pasy = 0.5
minz = -2
# maxz needs to be >= zmax
maxz = 100
pasz = 2
# Pass all the migration parameters in a dictionary to use them in functions called below
m_params = {'minx': minx, 'maxx': maxx,
            'pasx': pasx, 'pasy': pasy, 'miny': miny, 'maxy': maxy,
            'minz': minz, 'maxz': maxz, 'pasz': pasz, 'inc': inc, 'zmax': zmax}

#############################################
#### TODO: CREATE 2D grid profile from the 3d grid???
# TODO: and apply correction to the depths for the Earth's sphericity (SEE NOTES)

# Read the 3D numpy array of the RF amplitudes
with open('obs_amplitudes_matrix.npy', 'rb') as f:
    mObs_ = np.load(f)

# Define grid
grid_3d_a = mObs_
grid_3d_x = np.arange(minx, maxx + pasx, pasx)
grid_3d_y = np.arange(miny, maxy + pasy, pasy)
grid_3d_z = np.arange(minz, maxz + pasz, pasz)

# Interpolate in 3D
G_interpolated = RegularGridInterpolator((grid_3d_x, grid_3d_y, grid_3d_z), grid_3d_a)
pts = np.array([7.4, 46, 6])
VPinterp = G_interpolated(pts)

# Profile start and end
ref_pnt = np.array([[grid_3d_x[0], grid_3d_y[1]],
                    [grid_3d_x[-1], grid_3d_y[-5]]])
# lon0, lat0 = ref_pnt[0][0], ref_pnt[0][1]
lon0, lat0 = 8, 45.5
lon1, lat1 = 8, 48
# lon1, lat1 = ref_pnt[1][0], ref_pnt[1][1]

# Profile swath (km)?????
profile_swath = 100
# Profile azimuth
geoid = pyproj.Geod(ellps='WGS84')
profile_az, back_azimuth, profile_len_ = geoid.inv(lon0, lat0, lon1, lat1)
# Profile length (km)
profile_len = profile_len_/1000
# Two perpendicular azimuths
az1, az2 = get_perpendicular_azimuth(profile_az)


# Coordinates of the points along the profile knowing start and end of profile
n_extra_points = 100 # number of these points
geoid = Geod(ellps="WGS84")
extra_points = np.array(geoid.npts(lon0, lat0, lon1, lat1, n_extra_points))

# Create new lists of lon, lat, dep and amps (interpolated)
lon_points_along_prof = extra_points[:,0]
lat_points_along_prof = extra_points[:,1]

# For each point of the profile find the two end point perpendicular
# to them in a distance equal to the swath of the profile
amps = []
for i, lon in enumerate(lon_points_along_prof):
    # Two points perpendicular to the azimuth of the profile at each point of the profile
    lat_1, lon_1 = get_end_point(lat_points_along_prof[i], lon_points_along_prof[i], az1, profile_swath)
    lat_2, lon_2 = get_end_point(lat_points_along_prof[i], lon_points_along_prof[i], az2, profile_swath)
    n_extra_points_ = 5 # number of these points
    points_perpendicular_2_prof = np.array(geoid.npts(lon_1, lat_1, lon_2, lat_2, n_extra_points_))

    temp_lon = points_perpendicular_2_prof[:,0]
    temp_lat = points_perpendicular_2_prof[:,1]
    # interpolate for each of the 5 points perpendicular to the profile
    amps_matrix_temp = np.zeros((len(grid_3d_z)))
    nG = np.zeros((len(grid_3d_z)))

    for j, lon_ in enumerate(temp_lon):
        print(j)
        amps_temp = np.zeros((len(grid_3d_z)))
        for k, z in enumerate(grid_3d_z):
            print(z)
            point = np.array([lon_, temp_lat[j], z])
            VPinterp = G_interpolated(point)
            amps_temp[k] = VPinterp[0]
            print(VPinterp[0])
        amps_matrix_temp = amps_matrix_temp + amps_temp
        nG = nG + 1

    # add and divide by the number of stacks
    G = np.divide(amps_matrix_temp, nG)
    amps.append(G.tolist())
G2 = np.array(amps) # 2 dimensions

mObs = rf_mig.ccp_smooth(G2, m_params)
mObs[np.abs(mObs) < np.max(np.abs(mObs)) * 15 / 100] = 0
mObs = rf_mig.ccpFilter(mObs)
# plot_migration_sphr.plot_migration_profile(Gp=mObs, migration_param_dict=m_params, sta=sta,
#                                            work_directory=work_dir, filename=False)
# TOdO: continue from here and plot the 2D array


Gp = mObs
migration_param_dict = m_params
sta = rf_mig.read_stations_from_sac(path2rfs=path)
import rfmpy.core.migration as rf_mig_

sta, dxSta, dySta = rf_mig_.project_stations(sta=sta, ori_prof=profile_az,
                                            point_lat=lat0, point_lon=lon0)
work_directory=work_dir
filename=False

# Plot stations and profile
lons = [lon0,lon1]
lats = [lat0,lat1]
plt.plot(lons, lats, c='dodgerblue')
plt.scatter(lon0, lat0, c='dodgerblue', marker='s', edgecolor='k', s=50)
plt.scatter(lon1, lat1, c='dodgerblue', marker='o', edgecolor='k', s=50)
plt.scatter(sta["LONSTA"], sta["LATSTA"],
            c='r', marker='v', edgecolor='k', s=100)
plt.ylim(miny, maxy)
plt.xlim(minx, maxx)
plt.show()


xx = np.arange(0, profile_len, profile_len/n_extra_points)
zz = np.arange(minz, maxz + pasz, pasz)

XX, ZZ = np.meshgrid(xx, zz)

pal_col = work_directory + "/data/colormaps/vik.txt"
pal_col = pd.read_csv(pal_col, header=None, index_col=False, sep="\s+", names=["R", "G", "B"])
cm = LinearSegmentedColormap.from_list("blue2red", pal_col.values, len(pal_col))
c = np.min([np.max(Gp), 0.1])
c = 0.06
CL = 2

# PLOT
f = plt.figure(1, figsize=[10, 8])
gs0 = gridspec.GridSpec(nrows=1, ncols=1, figure=f,
                        hspace=0.08, right=0.91, left=0.09, bottom=0.08, top=0.96, )

ax = f.add_subplot(gs0[0])  # Ray tracing
# bwr, seismic, coolwarm, RdBu
m = ax.pcolormesh(XX, ZZ, Gp.T, cmap=cm, vmin=-c / CL, vmax=c / CL,
                  zorder=1, shading="auto")
add_colorbar(ax, m)
ax.scatter(sta["XSTA"].values, sta["ZSTA"].values,
           markersize, facecolors="grey", edgecolors="k",
           marker="v", lw=0.95, zorder=3, clip_on=False,
           label="Seismic stations", )

ax.set_aspect("equal")
ax.set_xlabel("x [km]", fontsize=fontsize)
ax.set_ylabel("z [km]", fontsize=fontsize)

majorLocator = MultipleLocator(5)
minorLocator = MultipleLocator(2.5)
ax.xaxis.set_major_locator(majorLocator)
ax.xaxis.set_minor_locator(minorLocator)
ax.set_xticks(np.arange(0, 300, 50))

majorLocator = MultipleLocator(5)
minorLocator = MultipleLocator(2.5)
ax.yaxis.set_major_locator(majorLocator)
# ax.yaxis.set_minor_locator(minorLocator)
ax.set_yticks(np.arange(0, 100, 20))
ax.set_ylim([100, 0])

ax.tick_params(axis="both", which="major", labelsize=fontsize)
ax.tick_params(axis="both", which="minor", labelsize=fontsize)
if filename:
    f.savefig(filename, dpi=200)

plt.show()
plt.close()

