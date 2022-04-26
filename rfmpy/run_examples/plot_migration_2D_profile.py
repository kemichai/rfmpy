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

# Define MIGRATION parameters
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



def dist_calc(loc1, loc2):
    """
    Function to calculate the distance in km between two points.

    Uses the flat Earth approximation. Better things are available for this,
    like `gdal <http://www.gdal.org/>`_.

    :type loc1: tuple
    :param loc1: Tuple of lat, lon, depth (in decimal degrees and km)
    :type loc2: tuple
    :param loc2: Tuple of lat, lon, depth (in decimal degrees and km)

    :returns: Distance between points in km.
    :rtype: float
    :author: Calum Chamberlain
    """
    R = 6371.009  # Radius of the Earth in km
    dlat = np.radians(abs(loc1[0] - loc2[0]))
    dlong = np.radians(abs(loc1[1] - loc2[1]))
    ddepth = abs(loc1[2] - loc2[2])
    mean_lat = np.radians((loc1[0] + loc2[0]) / 2)
    dist = R * np.sqrt(dlat ** 2 + (np.cos(mean_lat) * dlong) ** 2)
    dist = np.sqrt(dist ** 2 + ddepth ** 2)
    return dist

def get_azimuth(x1, y1, x2, y2):
  angle = 0.0;
  dx = x2 - x1
  dy = y2 - y1
  if x2 == x1:
    angle = np.pi / 2.0
    if y2 == y1 :
      angle = 0.0
    elif y2 < y1 :
      angle = 3.0 * np.pi / 2.0
  elif x2 > x1 and y2 > y1:
    angle = np.arctan(dx / dy)
  elif x2 > x1 and y2 < y1 :
    angle = np.pi / 2 + np.arctan(-dy / dx)
  elif x2 < x1 and y2 < y1 :
    angle = np.pi + np.arctan(dx / dy)
  elif x2 < x1 and y2 > y1 :
    angle = 3.0 * np.pi / 2.0 + np.arctan(dy / -dx)
  return (angle * 180 / np.pi)

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
    return perp_az1, perp_az2

#############################################
#### TODO: CREATE 2D grid profile from the 3d grid???
# TODO: and apply correction to the depths for the Earth's sphericity (SEE NOTES)

# Read the 3D numpy array of the RF amplitudes
with open('obs_amplitudes_matrix.npy', 'rb') as f:
    mObs = np.load(f)

# Define grid
grid_3d_a = mObs
grid_3d_x = np.arange(minx, maxx + pasx, pasx)
grid_3d_y = np.arange(miny, maxy + pasy, pasy)
grid_3d_z = np.arange(minz, maxz + pasz, pasz)

# Interpolate in 3D
G_interpolated = RegularGridInterpolator((grid_3d_x, grid_3d_y, grid_3d_z), grid_3d_a)
pts = np.array([7.4, 46, 6])
VPinterp = G_interpolated(pts)

# Profile start and end
ref_pnt = np.array([[grid_3d_x[0], grid_3d_y[0]],
                    [grid_3d_x[-1], grid_3d_y[-5]]])
lon0, lat0 = ref_pnt[0][0], ref_pnt[0][1]
lon1, lat1 = ref_pnt[1][0], ref_pnt[1][1]

# Profile swath (km)
profile_swath = 100
# Profile length (km)
profile_len = dist_calc([lon0, lat0, 0], [lon1, lat1, 0])
# Profile azimuth
profile_az = get_azimuth(lon0, lat0, lon1, lat1)
# Two perpendicular azimuths
az1, az2 = get_perpendicular_azimuth(profile_az)



# Coordinates of the points along the profile knowing start and end of profile
# TODO: double check why it crushes with more than 20...
n_extra_points = 10 # number of these points
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
    n_extra_points = 5 # number of these points
    points_perpendicular_2_prof = np.array(geoid.npts(lon_1, lat_1, lon_2, lat_2, n_extra_points))

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



# lon = extra_points[:,0]
# lat = extra_points[:,1]
# amps = []
# for i, ln in enumerate(lon):
#     amps_temp = []
#     for z in grid_3d_z:
#         print(z)
#         point = np.array([ln, lat[i], z])
#         VPinterp = G_interpolated(point)
#         amps_temp.append(VPinterp[0])
#     amps.append(amps_temp)
# amps = np.array(amps) # this is the G in two dimensions...



# this is the G in two dimensions...
#
# prof_dist, prof_dep = [], []
# cos_lat = np.cos(ref_pnt[0][1] * np.pi / 180)
# vec_ab = ref_pnt[1] - ref_pnt[0]
# vec_ab[0] *= cos_lat
# abs_ab = np.linalg.norm(vec_ab)
# for i in range(len(lon)):
#     loc_c = np.array([lon[i], lat[i]])
#     vec_ac = loc_c - ref_pnt[0]
#     vec_ac[0] *= cos_lat
#     abs_ac = np.linalg.norm(vec_ac)
#     cos = vec_ac.dot(vec_ab) / abs_ab / abs_ac
#     if abs_ac * (1 - cos ** 2) ** 0.5 > profile_width / 111.: continue
#     if cos < 0 or abs_ac * cos > abs_ab: continue
#     prof_dist.append(abs_ac * cos * 111)
#     prof_dep.append(grid_3d_z.tolist())
#
# prof_dep = np.array(prof_dep)


# TODO: need to find a way to stack multiple profiles
# TODO: the corrected depths need to be defined for each distance...


zz_corr = []
R = 6371.009
for i, dist in enumerate(prof_dist):
    correction = np.sqrt(dist**2 + R)
    print(correction)
    zz_corr = zz + correction

    XX, ZZ = np.meshgrid(prof_dist, zz_corr)
    # COLOR PALETTE AND COLORMAP
    # cork, bam, broc, vik
    pal_col = work_directory + "/data/colormaps/vik.txt"
    pal_col = pd.read_csv(pal_col, header=None, index_col=False, sep="\s+", names=["R", "G", "B"])
    cm = LinearSegmentedColormap.from_list("blue2red", pal_col.values, len(pal_col))
    c = np.min([np.max(Gp), 0.1])
    c = 0.06
    CL = 2

    # PLOT
    f = plt.figure(1, figsize=[10, 8])
    gs0 = gridspec.GridSpec(nrows=1, ncols=1, figure=f,
                            hspace=0.08, right=0.91, left=0.09, bottom=0.08, top=0.96,)

    ax = f.add_subplot(gs0[0])  # Ray tracing
    # bwr, seismic, coolwarm, RdBu
    m = ax.pcolormesh(XX, ZZ, Gp.T, cmap=cm, vmin=-c / CL, vmax=c / CL,
                      zorder=1, shading="auto")
    add_colorbar(ax, m)
    ax.scatter(sta["LONSTA"].values, sta["ZSTA"].values,
               markersize, facecolors="grey", edgecolors="k",
               marker="v", lw=0.95, zorder=3, clip_on=False,
               label="Seismic stations", )

    ax.set_aspect("equal")
    ax.set_xlabel("x [km]", fontsize=fontsize)
    ax.set_ylabel("z [km]", fontsize=fontsize)

    majorLocator = MultipleLocator(10)
    minorLocator = MultipleLocator(2.5)
    ax.xaxis.set_major_locator(majorLocator)
    ax.xaxis.set_minor_locator(minorLocator)
    # ax.set_xticks(np.arange(0, 140, 10))

    majorLocator = MultipleLocator(10)
    minorLocator = MultipleLocator(2.5)
    ax.yaxis.set_major_locator(majorLocator)
    ax.yaxis.set_minor_locator(minorLocator)
    ax.set_yticks(np.arange(10, 540, 100))
    ax.set_ylim([500, 0])

    ax.tick_params(axis="both", which="major", labelsize=fontsize)
    ax.tick_params(axis="both", which="minor", labelsize=fontsize)
    if filename:
        f.savefig(filename, dpi=200)

    plt.show()
    plt.close()






# Define profile for migration
profile_az = 0
# Define point (lon, lat) to project the station's coordinates with respect
profile_lon = sta["LONSTA"].mean()
profile_lat = sta["LATSTA"].mean()