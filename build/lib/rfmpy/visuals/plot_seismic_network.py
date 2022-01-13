from mpl_toolkits.basemap import Basemap
import matplotlib
matplotlib.use("Qt5Agg")
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


font = {'family': 'normal',
        'weight': 'normal',
        'size': 18}
matplotlib.rc('font', **font)
# Set figure width to 12 and height to 9
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 13
fig_size[1] = 8
plt.rcParams["figure.figsize"] = fig_size
subplot_rect = {'left': 0.08, 'right': 0.92, 'bottom': 0.08, 'top': 0.95, 'wspace': 0.1, 'hspace': 0.1}


MIN_LAT, MAX_LAT, MIN_LON, MAX_LON = (43, 50, 2, 15)
def map_stations():
    """

    """
fig = plt.figure()
# Map view
ax1 = plt.subplot2grid((3, 3), (0, 0), colspan=3, rowspan=3)
# fig.sca(ax1)
# Define Basemap
# Note: You can define the resolution of the map you just created. Higher
# resolutions take longer to create.
#    'c' - crude
#    'l' - low
#    'i' - intermediate
#    'h' - high
#    'f' - full
bmap = Basemap(llcrnrlon=MIN_LON, llcrnrlat=MIN_LAT, urcrnrlon=MAX_LON,
               urcrnrlat=MAX_LAT, resolution='i', projection='merc',
               lat_0=MIN_LAT, lon_0=MIN_LON, ax=ax1, fix_aspect=False)
# Draw some map elements on the map
bmap.drawcoastlines()
bmap.drawmapboundary(fill_color='white')
bmap.drawcountries()
bmap.fillcontinents(color='white',lake_color='white')
bmap.drawparallels(np.arange(MIN_LAT, MAX_LAT, 2), labels=[0, 0, 0, 0],
                   linewidth=0.5, dashes=[1, 10])
bmap.drawmeridians(np.arange(MIN_LON, MAX_LON, 2), labels=[0, 0, 0, 0],
                   linewidth=0.5, dashes=[1, 10])
xx, yy = bmap(lon, lat)
conf = bmap.scatter(xx, yy, edgecolor="k", alpha=1, s=dot_size, marker='o',
             label='Tibet 1D', c=dep, cmap='viridis', ax=ax1, zorder=2)


# Define and plot Mount Everest
MEx, MEy = bmap(86.9250, 27.9881)
bmap.scatter(MEx, MEy,
             color='k',
             marker='x',
             edgecolor='#ffffff',
             s=100,
             ax=ax1,
             zorder=101)
# Seismic network
HCx, HCy = bmap(HC_lon, HC_lat)
bmap.scatter(HCx, HCy,
             color='dodgerblue',
             marker='^',
             edgecolor='k',
             s=100,
             alpha=0.8,
             ax=ax1,
             zorder=101)
bmap.drawmapscale(
    83.5, 26.5, 83.5, 26.5,
    100,
    units='km', fontsize=10,
    yoffset=None,
    barstyle='simple', labelstyle='simple',
    fillcolor1='w', fillcolor2='#000000',
    fontcolor='#000000',
    zorder=101)
cax = fig.add_axes([0.6, 0.2, 0.2, 0.02])
cb = fig.colorbar(conf, cax=cax, orientation="horizontal",
                  extend='both',
                  ticks=[40, 60, 80, 100])
cb.set_label("Hypocentral depth (km)", fontsize=16)
plt.subplots_adjust(**subplot_rect)
plt.savefig('Map.png', bbox_inches="tight", format='png', dpi=300)
plt.show()