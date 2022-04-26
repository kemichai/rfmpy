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

plt.scatter(sta["LONSTA"], sta["LATSTA"],
            c='r', marker='v', edgecolor='k', s=100)
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

with open('obs_amplitudes_matrix.npy', 'rb') as f:
    G = np.load(f)

################
# Smoothing    #
################
# mObs = rf_mig.ccp_smooth(mObs, m_params)
# mObs[np.abs(mObs) < np.max(np.abs(mObs)) * 15 / 100] = 0
# mObs = rf_mig.ccpFilter(mObs)
#
# ################
# # Plotting     #
# ################
plot_migration_sphr.plot_migration_profile(Gp=mObs, migration_param_dict=m_params, sta=sta,
                                      work_directory=work_dir, filename=False)



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




#############################################
#### TODO: CREATE 2D grid profile from the 3d grid???
# TODO: and apply correction to the depths for the Earth's sphericity (SEE NOTES)
# Amplitudes
grid_3d_a = mObs
grid_3d_x = np.arange(minx, maxx + pasx, pasx)
grid_3d_y = np.arange(miny, maxy + pasy, pasy)
grid_3d_z = np.arange(minz, maxz + pasz, pasz)

G_interpolated = RegularGridInterpolator((grid_3d_x, grid_3d_y, grid_3d_z), grid_3d_a)
pts = np.array([7.4, 46, 6])
VPinterp = G_interpolated(pts)

# Profile start and end
ref_pnt = np.array([[grid_3d_x[0], grid_3d_y[0]],
                    [grid_3d_x[-1], grid_3d_y[-5]]])
profile_width = 100
profile_len = dist_calc([ref_pnt[0][0],ref_pnt[0][1], 0],
                        [ref_pnt[1][0],ref_pnt[1][1], 0])
# Coordinates of profile knowing start and end of profile
from pyproj import Geod
lon0, lat0 = ref_pnt[0][0], ref_pnt[0][1]
lon1, lat1 = ref_pnt[1][0], ref_pnt[1][1]
n_extra_points = 100
geoid = Geod(ellps="WGS84")
extra_points = np.array(geoid.npts(lon0, lat0, lon1, lat1, n_extra_points))

# Create new lists of lon, lat, dep and amps (interpolated)

lon = extra_points[:,0]
lat = extra_points[:,1]
amps = []
for i, ln in enumerate(lon):
    amps_temp = []
    for z in grid_3d_z:
        print(z)
        point = np.array([ln, lat[i], z])
        VPinterp = G_interpolated(point)
        amps_temp.append(VPinterp[0])
    amps.append(amps_temp)
amps = np.array(amps) # this is the G in two dimensions...

prof_dist, prof_dep = [], []
cos_lat = np.cos(ref_pnt[0][1] * np.pi / 180)
vec_ab = ref_pnt[1] - ref_pnt[0]
vec_ab[0] *= cos_lat
abs_ab = np.linalg.norm(vec_ab)
for i in range(len(lon)):
    loc_c = np.array([lon[i], lat[i]])
    vec_ac = loc_c - ref_pnt[0]
    vec_ac[0] *= cos_lat
    abs_ac = np.linalg.norm(vec_ac)
    cos = vec_ac.dot(vec_ab) / abs_ab / abs_ac
    if abs_ac * (1 - cos ** 2) ** 0.5 > profile_width / 111.: continue
    if cos < 0 or abs_ac * cos > abs_ab: continue
    prof_dist.append(abs_ac * cos * 111)
    prof_dep.append(grid_3d_z.tolist())

prof_dep = np.array(prof_dep)


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