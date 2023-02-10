"""
Calculate residuals between our moho estimates and from previous studies
such as Spada et al. or Grad and Tira

Location: Chavannes-pres-renens, CH
Date: Feb 2023
Author: Konstantinos Michailos
"""
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Qt5Agg")

def dist_calc(loc1, loc2):
    """..."""
    R = 6371.009  # Radius of the Earth in km
    dlat = np.radians(abs(loc1[0] - loc2[0]))
    dlong = np.radians(abs(loc1[1] - loc2[1]))
    ddepth = abs(loc1[2] - loc2[2])
    mean_lat = np.radians((loc1[0] + loc2[0]) / 2)
    dist = R * np.sqrt(dlat ** 2 + (np.cos(mean_lat) * dlong) ** 2)
    dist = np.sqrt(dist ** 2 + ddepth ** 2)
    return dist



lons = []
lats = []
deps = []
with open("moho_depths_all.dat", "r") as f:
    for line in f:
        ln = line.split(',')
        lon = float(ln[0])
        lat = float(ln[1])
        mh = float(ln[2])
        lons.append(lon)
        lats.append(lat)
        deps.append(mh)

# Spada et al
lons_Spada = []
lats_Spada = []
deps_Spada = []
with open("Spada_moho.dat", "r") as f:
    for line in f:
        ln = line.split(' ')
        lon = float(ln[0])
        lat = float(ln[1])
        mh = float(ln[2])
        lons_Spada.append(lon)
        lats_Spada.append(lat)
        deps_Spada.append(mh)

# Spada et al
lons_Grad = []
lats_Grad = []
deps_Grad = []
with open("grad.csv", "r") as f:
    for line in f:
        ln = line.split(',')
        lon = float(ln[0])
        lat = float(ln[1])
        mh = float(ln[2])
        if lon > 0 and lon < 25 and lat < 55 and lat > 40:
            lons_Grad.append(lon)
            lats_Grad.append(lat)
            deps_Grad.append(mh)

loc_array_sp = []
for i, j in enumerate(lons_Spada):
    loc_array_sp.append([lons_Spada[i], lats_Spada[i], deps_Spada[i]])
loc_array_grad = []
for i, j in enumerate(lons_Grad):
    loc_array_grad.append([lons_Grad[i], lats_Grad[i], deps_Grad[i]])


plt.scatter(lons_Spada, lats_Spada, c=deps_Spada, cmap='viridis',
            s=50, edgecolors="None",
            marker="v", label="Spada et al.", )
plt.scatter(lons, lats, c=deps, cmap='viridis',
            s=50, facecolors="grey", edgecolors="k",
            marker="o", label="Our study", )
plt.colorbar()
plt.legend()
plt.show()



plt.scatter(lons_Grad, lats_Grad, c=deps_Grad, cmap='viridis',
            s=50, edgecolors="None",
            marker="v", label="Grad and Tira", )
plt.scatter(lons, lats, c=deps, cmap='viridis',
            s=50, facecolors="grey", edgecolors="k",
            marker="o", label="Our study", )
plt.colorbar()
plt.legend()
plt.show()


lons_sp = []
lats_sp = []
dep_diff_spada = []
radius = 30
for i, longitude in enumerate(lons):
    # set depth to zero to find the distance from the surface
    location = [longitude, lats[i], 0]
    for sp_location in loc_array_sp:
        # set depth to zero to find the distance from the surface
        sp_location_ = sp_location[0:2] + [0]
        dist = dist_calc(location, sp_location)
        if dist < radius:
            print(f"Nearby observations:  {deps[i]}, {sp_location[-1]} !")
            depth_difference = deps[i] - sp_location[-1]
            lons_sp.append(sp_location_[0])
            lats_sp.append(sp_location_[1])
            dep_diff_spada.append(depth_difference)

lons_grad = []
lats_grad = []
dep_diff_grad = []
radius = 30
for i, longitude in enumerate(lons):
    # set depth to zero to find the distance from the surface
    location = [longitude, lats[i], 0]
    for grad_location in loc_array_grad:
        # set depth to zero to find the distance from the surface
        sp_location_ = grad_location[0:2] + [0]
        dist = dist_calc(location, grad_location)
        if dist < radius:
            print(f"Nearby observations:  {deps[i]}, {grad_location[-1]} !")
            depth_difference = deps[i] - grad_location[-1]
            lons_grad.append(grad_location[0])
            lats_grad.append(grad_location[1])
            dep_diff_grad.append(depth_difference)



x_obs = np.asarray(lons)
y_obs = np.asarray(lats)
x = np.asarray(lons_sp)
y = np.asarray(lats_sp)
v = np.asarray(dep_diff_spada)
# TODO: make this a function!!!
# def make_map(x, y, z):
# Plotting part of the function
font = {'family': 'normal',
        'weight': 'normal',
        'size': 14}
matplotlib.rc('font', **font)
# Set figure width to 12 and height to 9
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 12
fig_size[1] = 15
plt.rcParams["figure.figsize"] = fig_size

min_lat, max_lat, min_lon, max_lon = (40, 50, 0, 25)

fig = plt.figure()
ax1 = plt.subplot2grid((2,2), (0, 0), colspan=2, rowspan=2)
bmap = Basemap(llcrnrlon=min_lon, llcrnrlat=min_lat, urcrnrlon=max_lon,
               urcrnrlat=max_lat, resolution='l', projection='merc',
               lat_0=min_lat, lon_0=min_lon, ax=ax1, fix_aspect=False)
bmap.drawcoastlines()
bmap.drawmapboundary(fill_color='white')
bmap.fillcontinents(color='lightgray', lake_color='white')
bmap.drawparallels(np.arange(min_lat, max_lat, 2), labels=[1, 0, 0, 0],
                   linewidth=0.5, dashes=[1, 10])
bmap.drawmeridians(np.arange(min_lon, max_lon, 2), labels=[0, 0, 0, 1],
                   linewidth=0.5, dashes=[1, 10])
xi, yi = bmap(x, y)
x_obs_, yxobs_ = bmap(x_obs, y_obs)

conf = bmap.scatter(xi, yi, c=dep_diff_spada, s=150, cmap='RdBu', ax=ax1,
                    vmin=-20, vmax=20, zorder=101)
bmap.scatter(x_obs_, yxobs_, c='black', s=10, zorder=101)
cax = fig.add_axes([0.6, 0.2, 0.2, 0.02])
cb = fig.colorbar(conf, cax=cax, orientation="horizontal",
                  extend='both',
                  ticks=[-20, -10, 0, 10, 20])
cb.set_label("Moho_Michailos et al - Moho_Spada et al (km)", fontsize=16)
plt.show()


# Set figure width to 12 and height to 9
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 12
fig_size[1] = 9
plt.rcParams["figure.figsize"] = fig_size
# Plot histogram of the absolute differences in Moho depths
bins = np.arange(-25.5, 25.5, 1)
ax2 = plt.subplot2grid((1, 1), (0, 0), colspan=1, rowspan=1)
plt.axvspan(0, 30, facecolor='gray',alpha=0.1)
plt.axvline(x=0, c='black', linewidth=0.5)
ax2.hist(dep_diff_grad, bins, histtype='step', orientation='vertical',
         color='black',facecolor='black', alpha=0.8, linewidth=2.,
         linestyle='-', edgecolor='black',fill=False, label='vs Grad and Tira')
ax2.hist(dep_diff_spada, bins, histtype='step', orientation='vertical',
         color='dodgerblue',facecolor='black', alpha=0.8,
         linewidth=2., linestyle='--', edgecolor='dodgerblue',
         fill=False, label="vs Spada et al")
ax2.grid('True', linestyle="-", color='gray', linewidth=0.1, alpha=0.5)
ax2.set_ylabel('Number of observations', fontsize=20)
ax2.set_xlabel('Moho depth residuals (km)', fontsize=20)
ax2.tick_params(bottom=True, top=True, left=True, right=True)
ax2.legend(loc='best')
ax2.set_xlim([-25, 25])
plt.text(.74, 0.98, 'Moho estimates in this study are relatively deeper',
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax2.transAxes)
plt.tight_layout()
plt.show()
