"""
Script to create linear interpolation of the moho depths.
And then get values with the surface.
Meta aftes tis times mporoun na dior8w8oun mia mia...

Location: Chavannes-pres-renens, CH
Date: Feb 2023
Author: Konstantinos Michailos
"""
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


lon = []
lat = []
dep = []
with open('moho_depths_all.dat', 'r') as f:
        for line in f:
            if line.startswith('#'):
                print(line)
                continue
            else:
                ln = line.split(',')
                lon.append(float(ln[0]))
                lat.append(float(ln[1]))
                dep.append(float(ln[2]))


def interp_surf(x_, y_, z_, plot=False, **kwargs):
    """
    Gives the 2-D interpolated surface of exhumation rates.

    Use scipy interpolate to create a 2-D surface of a grid of
    observations.

    :type x_: list
    :param x_: list of coordinates (longitude)
    :type y_: list
    :param y_: list of coordinates (latitutes)
    :type z_: list
    :param z_: list of values to be interpolated
    :param kwargs: Any other arguments accepted by scipy.interpolate

    :returns: interpolated surface
    :rtype: function

    .. rubric:: Example        fig = plt.figure(figsize=(10, 6))
        ax = axes3d.Axes3D(fig)
        ax.plot_wireframe(x, y, z)
        ax.plot_surface(x, y, z, cmap=cm.viridis, alpha=0.2)
        ax.scatter3D(x, y, z, c='r')
        ax.set_zlim(.0, 12.0)
        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
        plt.show()

    >>> f = interp_surf(x,y,z)
    >>> znew = f(xnew, ynew)
    """
    # Convert to arrays
    x = np.asarray(x_)
    y = np.asarray(y_)
    z = np.asarray(z_)
    f = interpolate.interp2d(x, y, z, kind='linear', bounds_error=False)
    if plot:
        fig = plt.figure(figsize = (12,10))
        ax = plt.axes(projection='3d')

        X, Y = np.meshgrid(x, y)
        # Z =

        surf = ax.plot_surface(X, Y, Z, cmap = plt.cm.cividis)

        fig.colorbar(surf, shrink=0.5, aspect=8)
        plt.show()

    return f


def create_profiles(profile_points):
    from pyproj import Geod
    import pyproj
    # Profile start and end
    lon0, lat0 = profile_points[0][0], profile_points[0][1]
    lon1, lat1 = profile_points[1][0], profile_points[1][1]

    # Profile azimuth
    geoid = pyproj.Geod(ellps='WGS84')
    profile_az, back_azimuth, profile_len_ = geoid.inv(lon0, lat0, lon1, lat1)
    # Profile length (km)
    profile_len = profile_len_ / 1000

    # Coordinates of the points along the profile knowing start and end of profile
    # TODO: when I define a finer grid I won't need the * here!!!!!!
    n_extra_points = 25  # number of these points
    print("Number of points along the profile: ", n_extra_points, " Length of profile: ", profile_len)

    geoid = Geod(ellps="WGS84")
    extra_points = np.array(geoid.npts(lon0, lat0, lon1, lat1, n_extra_points))
    # Create new lists of lon, lat, dep and amps (interpolated)
    lon_points_along_prof = extra_points[:, 0]
    lat_points_along_prof = extra_points[:, 1]

    return lon_points_along_prof, lat_points_along_prof

profile_A = np.array([[6, 49], [11.5, 44]])
# profile_A = np.array([[3, 44.1], [9, 44.8]])
# profile_A = np.array([[11, 45.5], [22, 50]])

lons, lats = create_profiles(profile_A)


f = interp_surf(lon, lat, dep)
for i, xnew in enumerate(lons):
    znew = f(xnew, lats[i])
    print(xnew, lats[i], znew[0])
