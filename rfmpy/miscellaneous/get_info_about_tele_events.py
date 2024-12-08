"""
Reads quakeml catalog to obtain info about the events.

Location: Chavannes-pres-renens, CH
Date: Jan 2022
Author: Konstantinos Michailos
"""
import numpy as np
import platform
from obspy import read_inventory, read_events, UTCDateTime as UTC

# Set up parameters and paths
if platform.node().startswith('kmichailos-laptop'):
    data_root_dir = '/media/kmichailos/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/home/kmichailos/Desktop/codes/bitbucket'
    desktop_dir = '/home/kmichailos/Desktop'
    hard_drive_dir = '/media/kmichailos/SEISMIC_DATA/'
else:
    data_root_dir = '/media/kmichall/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/home/kmichall/Desktop/Codes/bitbucket'
    desktop_dir = '/home/kmichall/Desktop'
    hard_drive_dir = '/media/kmichall/SEISMIC_DATA/'

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

# STore for gmt
def qml2gmt(cat, output_name, sort_by='depth'):
    """Read catalog output .dat file for GMT"""
    outdir = codes_root_dir + '/him_seismicity/maps/'
    if sort_by == 'depth':
        cat.events.sort(key=lambda e: e.origins[-1].depth)
        print('Catalog is sorted by depth!!! Change arg sort_by if that is not cool...')
    else:
        cat.events.sort(key=lambda e: e.origins[-1].time)

    for ev in cat:
        lat = ev.origins[-1].latitude
        lon = ev.origins[-1].longitude
        dep = ev.origins[-1].depth / 1000
        mag = ev.magnitudes[-1].mag
        nsta = ev.origins[-1].quality.used_station_count
        az_gap = ev.origins[-1].quality.azimuthal_gap
        with open(outdir + output_name + '.dat', 'a') as f:
            f.write('{} {} {} {} {} {} \n'.format(lon, lat, dep, mag, nsta, az_gap))
    return


cat = read_events(desktop_dir + '/Codes/github/rfmpy/rfmpy/data/catalog/*.xml')

net_loc = [47.0, 10.0, 0.0]
for ev in cat:
    ev_loc = [ev.origins[0].latitude, ev.origins[0].longitude, ev.origins[0].depth/1000]
    dist = dist_calc(net_loc, ev_loc)
    if ev.magnitudes[0].mag < 6.0 and ev.origins[0].depth/1000 > 100 and dist < 10000:
        print(str(dist) + ', ' + str(round(ev.origins[0].depth/1000,2)) + ', ' + str(round(ev.magnitudes[0].mag,2)) + ', ' + str(ev.comments[0].text))
