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
with open('all_G3.npy', 'rb') as f:
    G3_ = np.load(f)

# read stations
sta = rf_mig.read_stations_from_sac(path2rfs=path)

profile_A = np.array([[8, 45.5], [8.6, 48]])
G2, sta, xx, zz = create_2d_profile(G3_, m_params, profile_A, sta,
                                    swath=300, plot=True)


mObs = rf_mig.ccp_smooth(G2, m_params)
# mObs[np.abs(mObs) < np.max(np.abs(mObs)) * 15 / 100] = 0
mObs = rf_mig.ccpFilter(mObs)
plot_migration_sphr.plot_migration_profile(Gp=mObs, xx=xx, zz=zz, migration_param_dict=m_params, sta=sta,
                                           work_directory=work_dir, filename=False)


for i, z_slice in enumerate(mObs[0:]):
    for z in z_slice:
        print(z)

