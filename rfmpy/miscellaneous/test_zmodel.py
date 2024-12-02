"""
Doublecheck codes read correctly the 3D velocity model of Zhu et al., 2015). + IASP91

Issue reported from DK and GH.

Potential issues:
- Are the Vp and Vs values in the model changed during the scan?
  We can check this by plotting the vp and vs values along the same line and then we can see how it goes.
- Some depth value might be missing in the model and this "squeezes" the migration,
  because we almost see that the depth scale is compressed.
- The shape of the ray tracing itself is wrong for something that is probably related to velocity.
  This can result in the discontinuity red spots not being connected, the discontinuity being thick.

"""
from scipy.interpolate import LinearNDInterpolator
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.interpolate import RegularGridInterpolator


font = {'family': 'normal',
        'weight': 'normal',
        'size': 18}
matplotlib.rc('font', **font)
# Set figure width to 12 and height to 9
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 10
fig_size[1] = 8
plt.rcParams["figure.figsize"] = fig_size


min_lon=0
max_lon=32
min_lat=40
max_lat=55


work_dir = os.getcwd()
path_ = work_dir.split('/')
path = '/'.join(path_)

# Path to file

path_zmodel_m60 = path + '/data/ZMODEL_M60/'

# Read x, y, z, etc .nps file of ZMODEL_M60 velocity model
parameters = np.load(path_zmodel_m60 + 'z_model_m60.npz')
parameters.items()
longitudes = parameters["longitudes"].tolist()
latitudes = parameters["latitudes"].tolist()
depths = parameters["depths"].tolist()
p_velocities = parameters["Vp"].tolist()
s_velocities = parameters["Vs"].tolist()

points_list = []
vp_values = []
vs_values = []
for i, lon in enumerate(longitudes):
    if lon < max_lon and lon > min_lon and latitudes[i] > min_lat and latitudes[i] < max_lat:

        point = [lon, latitudes[i], depths[i]]
        if depths[i] == 10.0:
            # Extend velocity model on at least 5 km above sea level
            # IMPORTANT TO NOTE
            additional_point = [lon, latitudes[i], -5]
            points_list.append(additional_point)
            vp_values.append(p_velocities[i])
            vs_values.append(s_velocities[i])


        points_list.append(point)
        vp_values.append(p_velocities[i])
        vs_values.append(s_velocities[i])


points = np.array(points_list)
values_p = np.array(vp_values)
values_s = np.array(vs_values)
# rescale here is important for making the steps sharp (look at the following link:
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.LinearNDInterpolator.html
print("|-----------------------------------------------|")
print("| Interpolating 3D ZMODEL_M60...                |")
print("| This might take a while...                    |")
liner_interpolation_of_velocities_p = LinearNDInterpolator(points, values_p, rescale=True)
liner_interpolation_of_velocities_s = LinearNDInterpolator(points, values_s, rescale=True)
print("| Interpolated 3D ZMODEL_M60...                 |")

from scipy.interpolate import LinearNDInterpolator
import os

# Read x, y, z, etc .txt file of EPcrust velocity model
longitudes = []
latitudes = []
topos = []
thick_sediments = []
thick_upper = []
thick_lower = []
vp_sediments = []
vp_upper = []
vp_lower = []
vs_sediments = []
vs_upper = []
vs_lower = []

work_dir = os.getcwd()
# Path to EPcrust file
path_epcrust = work_dir + '/data/EPcrust/'
try:
    with open(path_epcrust + 'EPcrust_0_5.txt', 'r') as f:
        for line in f:
            if line.startswith('#'):
                print('|Reading EPcrust velocity model...              |')
                continue
            else:
                ln = line.split()
                lon_ = float(ln[0])
                lat_ = float(ln[1])
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
                if lon_ < max_lon and lon_ > min_lon and lat_ > min_lat and lat_ < max_lat:
                    longitudes.append(lon_)
                    latitudes.append(lat_)
                    topos.append(topo)
                    thick_sediments.append(thick_sed)
                    thick_upper.append(thick_upp)
                    thick_lower.append(thick_low)
                    vp_sediments.append(vp_sed)
                    vp_upper.append(vp_upp)
                    vp_lower.append(vp_low)
                    vs_sediments.append(vs_sed)
                    vs_upper.append(vs_upp)
                    vs_lower.append(vs_low)
except Exception as e:
    raise type(e)('>>> Move to the top directory of the repository!')

lon = np.array(longitudes)
lat = np.array(latitudes)
ele = np.array(topos)

# P wave velocities
vp_sediments = np.array(vp_sediments)
vp_upper = np.array(vp_upper)
vp_lower = np.array(vp_lower)
# S wave velocities
vs_sediments = np.array(vs_sediments)
vs_upper = np.array(vs_upper)
vs_lower = np.array(vs_lower)
# Thickness of three layers
thick_sediments = np.array(thick_sediments)
thick_upper = np.array(thick_upper)
thick_lower = np.array(thick_lower)

# Define depth profiles for each EPcrust's grid points
points = []
p_velocities = []
s_velocities = []
for i, _ in enumerate(lon):

    # Extend velocity model on at least 5 km above sea level
    z_0_ = -5.0
    point0_ = [_, lat[i], z_0_]
    points.append(point0_)
    # Checked that there is always vp and vp sedimentary.
    p_velocities.append(vp_sediments[i])
    s_velocities.append(vs_sediments[i])

    # First point at Earth's surface (including the topography)
    # If we do not include the topography the Moho depth value will be affected.
    elevation = ele[i]
    # If elevation is above sea surface
    if elevation > 0.0:
        # We define the depth profile starting from sea level (z = 0) with positive values as we move deeper.
        # And negative values for the topography above sea level.
        # So when elevation has a positive value we change it to negative.
        z_0 = (-1) * elevation
    elif elevation < 0.0:
        # If elevation is below sea level it will be negative. So similarly to above we multiply with -1
        # to change the value to positive and get the point where the velocity depth profile begins.
        z_0 = (-1) * elevation
    else:
        z_0 = 0.0
    point0 = [_, lat[i], z_0]
    points.append(point0)
    p_velocities.append(vp_sediments[i])
    s_velocities.append(vs_sediments[i])
    # Second point at the lower limit of the sediments.
    z_1 = thick_sediments[i] + z_0
    point1 = [_, lat[i], z_1]
    points.append(point1)
    p_velocities.append(vp_sediments[i])
    s_velocities.append(vs_sediments[i])
    # Third point at the lower limit of the sediments with the velocity below.
    z_2 = thick_sediments[i] + 0.01 + z_0
    point2 = [_, lat[i], z_2]
    points.append(point2)
    p_velocities.append(vp_upper[i])
    s_velocities.append(vs_upper[i])
    # Fourth point at the lower part of the upper crust.
    z_3 = thick_sediments[i] + thick_upper[i] + z_0
    point3 = [_, lat[i], z_3]
    points.append(point3)
    p_velocities.append(vp_upper[i])
    s_velocities.append(vs_upper[i])
    # Fifth point at the lower part of the upper crust...
    z_4 = thick_sediments[i] + thick_upper[i] + 0.01 + z_0
    point4 = [_, lat[i], z_4]
    points.append(point4)
    p_velocities.append(vp_lower[i])
    s_velocities.append(vs_lower[i])
    # Sixth point at the bottom of the crust...
    z_5 = thick_sediments[i] + thick_upper[i] + thick_lower[i] + z_0
    point5 = [_, lat[i], z_5]
    points.append(point5)
    p_velocities.append(vp_lower[i])
    s_velocities.append(vs_lower[i])
    # Seventh point at the bottom of the crust with mantle's velocity
    z_6 = thick_sediments[i] + thick_upper[i] + thick_lower[i] + 0.01 + z_0
    point6 = [_, lat[i], z_6]
    points.append(point6)
    # iasp91 values
    # p_velocities.append(8.05)
    # s_velocities.append(4.45)
    # epcrust lower crust
    p_velocities.append(vp_lower[i])
    s_velocities.append(vs_lower[i])
    # Eighth point at the mantle...
    z_7 = 120
    point7 = [_, lat[i], z_7]
    points.append(point7)
    # iasp91 values
    # p_velocities.append(8.1)
    # s_velocities.append(4.5)
    # epcrust lower crust
    p_velocities.append(vp_lower[i] + 0.05)
    s_velocities.append(vs_lower[i] + 0.05)

points = np.array(points)
values_p = np.array(p_velocities)
values_s = np.array(s_velocities)
# rescale here is important for making the steps sharp (look at the following link:
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.LinearNDInterpolator.html
liner_interpolation_of_velocities_p_ep = LinearNDInterpolator(points, values_p, rescale=True)
liner_interpolation_of_velocities_s_ep = LinearNDInterpolator(points, values_s, rescale=True)


def get_iasp91(x_, y, z, zmoho):
    """
    Retrieves P-wave, S-wave velocities and depths
    from IASPEI91 global velocity model.

    :type x_: numpy.array
    :param x_: Numpy array of x values of the grid points.
    :type y_: numpy.array
    :param y_: Numpy array of y values of the grid points.
    :type z: numpy.array
    :param z: Numpy array of z values of the grid points.
    :type zmoho: int
    :param zmoho: Moho depth in km.

    :rtype: numpy.ndarrays
    :returns: Array of P-wave, S-wave velocities and their depths.
    """

    R = 6371  # Earth's radius
    x = (R - z) / R
    VP = np.zeros((x_.size, y.size, z.size))
    VS = np.zeros((x_.size, y.size, z.size))
    for i in range(z.size):
        if z[i] < 20:
            VP[:, :, i] = 5.8
            VS[:, :, i] = 3.36
        elif z[i] < zmoho:
            VP[:, :, i] = 6.5
            VS[:, :, i] = 3.75
        elif z[i] < 120.0:
            VP[:, :, i] = 8.78541 - 0.74953 * x[i]
            VS[:, :, i] = 6.706231 - 2.248585 * x[i]
        elif z[i] < 210.0:
            VS[:, :, i] = 5.75020 - 1.27420 * x[i]
            VP[:, :, i] = 25.41389 - 17.69722 * x[i]
        elif z[i] < 410.0:
            VS[:, :, i] = 15.24213 - 11.08552 * x[i]
            VP[:, :, i] = 30.78765 - 23.25415 * x[i]
        elif z[i] < 660.0:
            VP[:, :, i] = 29.38896 - 21.40656 * x[i]
            VS[:, :, i] = 17.70732 - 13.50652 * x[i]
        elif z[i] < 760.0:
            VP[:, :, i] = 25.96984 - 16.93412 * x[i]
            VS[:, :, i] = 20.76890 - 16.53147 * x[i]
        elif z[i] < 2740.0:
            VP[:, :, i] = (25.1486 - 41.1538 * x[i] + 51.9932 * x[i] * x[i] - 26.6083 * x[i] * x[i] * x[i])
            VS[:, :, i] = (12.9303 - 21.2590 * x[i] + 27.8988 * x[i] * x[i] - 14.1080 * x[i] * x[i] * x[i])
        elif z[i] < 2889.0:
            VP[:, :, i] = 14.49470 - 1.47089 * x[i]
            VS[:, :, i] = 8.16616 - 1.58206 * x[i]
        elif z[i] < 5153.9:
            VP[:, :, i] = 10.03904 + 3.75665 * x[i] - 13.67046 * x[i] * x[i]
            VS[:, :, i] = 1.0e-20
        elif z[i] < 6371.0:
            VP[:, :, i] = 11.24094 - 4.09689 * x[i] * x[i]
            VS[:, :, i] = 3.56454 - 3.45241 * x[i] * x[i]
        else:
            VP[:, :, i] = -1
            VS[:, :, i] = -1.0
    return VP, VS


zmoho = 35
minx = -7.0 # degrees 2.optional:2 or -4
maxx = 40.0 # degrees 2.optional:30 or 38
pasx = 0.38 # degrees
miny = 38.0 # degrees 2.optional:41 or 38
maxy = 60.0 # degrees 2.optional:51 or 54
pasy = 0.27 # degrees

x = np.arange(minx, maxx, pasx)
y = np.arange(miny, maxy, pasy)
z_ = np.linspace(-5, 1000, 50)
VP, VS = get_iasp91(x, y, z_, zmoho)
# Interpolate
P_vel_3D_grid = RegularGridInterpolator((x, y, z_), VP)
S_vel_3D_grid = RegularGridInterpolator((x, y, z_), VS)



depths = np.linspace(-5, 1000, 50)
vel_zhu_s, vel_zhu_p = [], []
vel_ep_s, vel_ep_p = [], []
vel_iasp_s, vel_iasp_p = [], []
for d in depths:
    pts = np.array([5, 45, d])
    vel_zhu_s.append(liner_interpolation_of_velocities_s(pts))
    vel_zhu_p.append(liner_interpolation_of_velocities_p(pts))
    vel_ep_s.append(liner_interpolation_of_velocities_s_ep(pts))
    vel_ep_p.append(liner_interpolation_of_velocities_p_ep(pts))
    vel_iasp_p.append(P_vel_3D_grid(pts))
    vel_iasp_s.append(S_vel_3D_grid(pts))




ax1 = plt.subplot2grid((1, 2), (0, 0), colspan=2)
ax1.plot(vel_zhu_p, depths, zorder=2, color='k', linestyle='-', label='P - Zhu et al. 2015')
ax1.plot(vel_zhu_s, depths, zorder=2, color='r', linestyle='-', label='S - Zhu et al. 2015')

ax1.plot(vel_iasp_p, depths, zorder=2, color='k', linestyle='--', label='P - IASP91')
ax1.plot(vel_iasp_s, depths, zorder=2, color='r', linestyle='--', label='S - IASP91')

ax1.set_ylabel('Depth (km)', fontsize=18)
ax1.set_xlabel('Velocity (km/s)', fontsize=18)
ax1.set_ylim([-5, 1000])
plt.legend(loc="upper right", markerscale=1., scatterpoints=1, fontsize=14)
plt.gca().invert_yaxis()
plt.show()


