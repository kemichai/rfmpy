"""
Location: Chavannes-pres-renens, CH
Date: April 2022
Author: Konstantinos Michailos
"""

import rfmpy.core.migration_sphr as rf_mig
import rfmpy.utils.migration_plots_spher as plot_migration_sphr
import numpy as np
import os
import matplotlib.pyplot as plt


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
# Read the 3D numpy array of the RF amplitudes
with open('obs_amplitudes_matrix.npy', 'rb') as f:
    G3_ = np.load(f)

# read stations
sta = rf_mig.read_stations_from_sac(path2rfs=path)

profile_A = np.array([[8, 45.5], [8.6, 48]])
G2, sta, xx, zz = plot_migration_sphr.create_2d_profile(G3_, m_params, profile_A, sta, swath=30, plot=True)


mObs = rf_mig.ccp_smooth(G2, m_params)
# mObs[np.abs(mObs) < np.max(np.abs(mObs)) * 15 / 100] = 0
mObs = rf_mig.ccpFilter(mObs)
plot_migration_sphr.plot_migration_profile(Gp=mObs, xx=xx, zz=zz, migration_param_dict=m_params, sta=sta,
                                           work_directory=work_dir, filename=False)


for i, x in enumerate(xx):
    for j, z in enumerate(zz):
        print(x, z, mObs[i,j])



# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.patches import Wedge
#
# def perp( a ) :
#     ##from https://stackoverflow.com/a/3252222/2454357
#     b = np.empty_like(a)
#     b[0] = -a[1]
#     b[1] = a[0]
#     return b
#
#
# def seq_intersect(a1,a2, b1,b2) :
#     ##from https://stackoverflow.com/a/3252222/2454357
#     da = a2-a1
#     db = b2-b1
#     dp = a1-b1
#     dap = perp(da)
#     denom = np.dot( dap, db)
#     num = np.dot( dap, dp )
#     return (num / denom.astype(float))*db + b1
#
# def angle(a1, a2, b1, b2):
#     ##from https://stackoverflow.com/a/16544330/2454357
#     x1, y1 = a2-a1
#     x2, y2 = b2-b1
#     dot = x1*x2 + y1*y2      # dot product between [x1, y1] and [x2, y2]
#     det = x1*y2 - y1*x2      # determinant
#     return np.arctan2(det, dot)  # atan2(y, x) or atan2(sin, cos)
#
#
# def draw_wedge(
#     ax, r_min = 0.3, r_max = 0.5, t_min = np.pi/4, t_max = 3*np.pi/4
#     ):
#
#     ##some data
#     R = np.random.rand(100)*(r_max-r_min)+r_min
#     T = np.random.rand(100)*(t_max-t_min)+t_min
#     ax.scatter(T,R)
#
#     ##compute the corner points of the wedge:
#     axtmin = 0
#
#     rs = np.array([r_min,  r_max,  r_min, r_max, r_min, r_max])
#     ts = np.array([axtmin, axtmin, t_min, t_min, t_max, t_max])
#
#     ##display them in a scatter plot
#     ax.scatter(ts, rs, color='r', marker='x', lw=5)
#
#     ##from https://matplotlib.org/users/transforms_tutorial.html
#     trans = ax.transData + ax.transAxes.inverted()
#
#     ##convert to figure cordinates, for a starter
#     xax, yax = trans.transform([(t,r) for t,r in zip(ts, rs)]).T
#
#     for i,(x,y) in enumerate(zip(xax, yax)):
#         ax.annotate(
#             str(i), (x,y), xytext = (x+0.1, y), xycoords = 'axes fraction',
#             arrowprops = dict(
#                 width=2,
#
#             ),
#         )
#
#
#     ##compute the angles of the wedge:
#     tstart = np.rad2deg(angle(*np.array((xax[[0,1,2,3]],yax[[0,1,2,3]])).T))
#     tend = np.rad2deg(angle(*np.array((xax[[0,1,4,5]],yax[[0,1,4,5]])).T))
#
#     ##the center is where the two wedge sides cross (maybe outside the axes)
#     center=seq_intersect(*np.array((xax[[2,3,4,5]],yax[[2,3,4,5]])).T)
#
#     ##compute the inner and outer radii of the wedge:
#     rinner = np.sqrt((xax[1]-center[0])**2+(yax[1]-center[1])**2)
#     router = np.sqrt((xax[2]-center[0])**2+(yax[2]-center[1])**2)
#
#     wedge = Wedge(center,
#                   router, tstart, tend,
#                   width=router-rinner,
#                   #0.6,tstart,tend,0.3,
#                   transform=ax.transAxes, linestyle='--', lw=3,
#                   fill=False, color='red')
#     ax.add_artist(wedge)
#
#
#
#
# fig = plt.figure(figsize=(8,4))
#
# # ax1 = fig.add_subplot(121, projection='polar')
# ax2 = fig.add_subplot(122, projection='polar')
#
# ##reducing the displayed theta and r ranges in second axes:
# ax2.set_thetamin(0)
# ax2.set_thetamax(40)
#
# ## ax.set_rmax() does not work as one would expect -- use ax.set_ylim() instead
# ## from https://stackoverflow.com/a/9231553/2454357
# ax2.set_ylim([0.0,1.0])
# ax2.set_rorigin(-0.2)
#
# #from https://stackoverflow.com/a/41823326/2454357
# fig.canvas.draw()
#
# # draw_wedge(ax1)
# draw_wedge(ax2, t_min=np.deg2rad(15), t_max=np.deg2rad(35))
#
# plt.show()
#
#
#

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import FuncFormatter

minTheta = 0.1; maxTheta = 0.55

fig = plt.figure()
ax = fig.add_subplot(111, projection='polar')

#create four tick labels (in radians) dynamically based on the theta range
ticks = np.linspace(minTheta, maxTheta, 4)
ax.set_xticks(ticks)

#disable or enable the following line to change the tick display format*
rad2fmt = lambda x,pos : f"{np.rad2deg(x):.2f}°"
ax.xaxis.set_major_formatter(FuncFormatter(rad2fmt))

#Adjust the sector window: you must call these AFTER setting the ticks, since setting the ticks
#actually adjusts the theta range window. And it must be in degrees.
ax.set_thetamin(np.rad2deg(minTheta))
ax.set_thetamax(np.rad2deg(maxTheta))

plt.show()



"""
Demo of the floating axes.

This demo shows features of functions in floating_axes:
    * Using scatter function and bar function with changing the
      shape of the plot.
    * Using GridHelperCurveLinear to rotate the plot and set the
      boundary of the plot.
    * Using FloatingSubplot to create a subplot using the return
      value from GridHelperCurveLinear.
    * Making sector plot by adding more features to GridHelperCurveLinear.
"""
from matplotlib.transforms import Affine2D
import mpl_toolkits.axisartist.floating_axes as floating_axes
import numpy as np
import mpl_toolkits.axisartist.angle_helper as angle_helper
from matplotlib.projections import PolarAxes
from mpl_toolkits.axisartist.grid_finder import (FixedLocator, MaxNLocator,
                                                 DictFormatter)
import matplotlib.pyplot as plt


def setup_axes2(fig, rect):
    """
    With custom locator and formatter.
    Note that the extreme values are swapped.
    """
    tr = PolarAxes.PolarTransform()

    pi = np.pi
    angle_ticks = [(0, r"$0$"),
                   (.25*pi, r"$\frac{1}{4}\pi$"),
                   (.5*pi, r"$\frac{1}{2}\pi$")]
    grid_locator1 = FixedLocator([v for v, s in angle_ticks])
    tick_formatter1 = DictFormatter(dict(angle_ticks))

    grid_locator2 = MaxNLocator(2)

    grid_helper = floating_axes.GridHelperCurveLinear(
        tr, extremes=(.5*pi, 0, 2, 1),
        grid_locator1=grid_locator1,
        grid_locator2=grid_locator2,
        tick_formatter1=tick_formatter1,
        tick_formatter2=None)

    ax1 = floating_axes.FloatingSubplot(fig, rect, grid_helper=grid_helper)
    fig.add_subplot(ax1)

    # create a parasite axes whose transData in RA, cz
    aux_ax = ax1.get_aux_axes(tr)

    aux_ax.patch = ax1.patch  # for aux_ax to have a clip path as in ax
    ax1.patch.zorder = 0.9  # but this has a side effect that the patch is
    # drawn twice, and possibly over some other
    # artists. So, we decrease the zorder a bit to
    # prevent this.

    return ax1, aux_ax




##########################################################
fig = plt.figure(1, figsize=(8, 4))
fig.subplots_adjust(wspace=0.3, left=0.05, right=0.95)


ax2, aux_ax2 = setup_axes2(fig, 132)
# theta = np.random.rand(10)*.5*np.pi
# radius = np.random.rand(10) + 1.
theta = (8 + np.random.rand(10)*(14 - 8))*15.  # in degrees
radius = np.random.rand(10)*14000.
aux_ax2.scatter(theta, radius)




plt.show()

