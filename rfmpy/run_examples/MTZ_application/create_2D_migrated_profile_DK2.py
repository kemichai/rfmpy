"""
Location: Chavannes-pres-renens, CH
Date: April 2022
Author: Konstantinos Michailos
"""
import rfmpy.utils.migration_mtz as mtz
import numpy as np
import os
from obspy.geodetics import degrees2kilometers, kilometers2degrees
import rfmpy.utils.migration_plots_spher as plot_migration_sphr

# Define paths
work_dir = os.getcwd()

path2grid = '/home/' + work_dir.split('/')[2] + '/Desktop/mtz_example/'

# Read the 3D grid (epcrust.npy) of stacked migrated RF amplitudes.
with open(path2grid + 'example.npy', 'rb') as f:
    mObs = np.load(f)

## Define MIGRATION parameters
# Ray-tracing parameters
inc = 2  # km
zmax = 800 # km
# Determine study area (x -> perpendicular to the profile)
minx = -13.0 # degrees 2.optional:2 or -4
maxx = 46.0 # degrees 2.optional:30 or 38
pasx = 0.26 # degrees oldest 0.38
miny = 30.0 # degrees 2.optional:41 or 38
maxy = 64.0 # degrees 2.optional:51 or 54
pasy = 0.18 # degrees oldest 0.27
minz = -5 # km
# maxz needs to be >= zmax
maxz = 800 # km
pasz = 2 # km
# Pass all the migration parameters in a dictionary to use them in functions called below
m_params = {'minx': minx, 'maxx': maxx,
            'pasx': pasx, 'pasy': pasy, 'miny': miny, 'maxy': maxy,
            'minz': minz, 'maxz': maxz, 'pasz': pasz, 'inc': inc, 'zmax': zmax}

# Define path to RFs
path = '/home/' + work_dir.split('/')[2] + '/Desktop/mtz_example/RF/'

# Read station coordinates from the rfs (sac files) in a pandas dataframe
sta = mtz.read_stations_from_sac(path2rfs=path)

profile_A = np.array([[8, 50.5], [30, 45.2]])
prof_name = 'Cross-section_A_and_A'

G2_, sta_, xx, zz = plot_migration_sphr.create_2d_profile(mObs, m_params, profile_A, sta, swath=300, plot=True)
mObs = mtz.ccp_smooth(G2_, m_params)
mObs = mtz.ccpFilter(mObs)

# File for creating cross-sections with GMT
for i, x in enumerate(xx):
    for j, z in enumerate(zz):
        print(kilometers2degrees(x), z, mObs[i,j])
        with open(path2grid + prof_name + '.txt', 'a') as of:
            of.write('{} {} {} \n'.
                     format(kilometers2degrees(x), z, mObs[i, j]))
