"""
Code for calculating 3D receiver function migrations...
-
=============================================
Requirements:
    * obspy
    * scipy
    * pandas

=============================================

Location: Chavannes-pres-renens, CH
Date: Mar 2022
Author: Konstantinos Michailos
"""
import rfmpy.utils.migration_mtz as mtz
import os
import time
import multiprocessing as mp
from obspy import Stream


# Start a timer to keep a track how long the calculations take
t_beg = time.time()

# Define working directory
work_dir = os.getcwd()

# Define path to RFs
path = '/home/' + work_dir.split('/')[2] + '/Desktop/mtz_example/RF/'

# Read station coordinates from the rfs (sac files) in a pandas dataframe
sta = mtz.read_stations_from_sac(path2rfs=path)

# Read RFs
stream = mtz.read_traces_sphr(path2rfs=path, sta=sta)


# Number of cpus to use
n_cpu = mp.cpu_count() 
# Number of traces to process
num_traces = len(stream)
# Number of substreams (set one per cpu available here)
num_substreams = n_cpu
# Calculate the number of traces per substream
traces_per_substream = num_traces // (num_substreams)
# Create a list to store substreams
substreams = []
# Divide the stream into substreams
for i in range(num_substreams):
    start_index = i * traces_per_substream
    end_index = (i + 1) * traces_per_substream if i < num_substreams - 1 else num_traces
    substream = Stream(traces=stream[start_index:end_index])
    substreams.append(substream)

# Print the number of traces in each substream
for i, substream in enumerate(substreams):
    print(f"Substream {i + 1} includes {len(substream)} traces.")

# =================================================== #
# Define MIGRATION parameters
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

# Read velocity model
Vp, Vs = mtz.read_vel_model(m_params, 'zmodel_m60')

def wrapper(args):
    return mtz.tracing_3D_sphr_parallel(*args)

# Parallel processing
pool = mp.Pool(processes=n_cpu)
args_list = [(sub_stream, m_params, Vp, Vs) for sub_stream in substreams]
result_list = pool.map(wrapper, args_list)

# Add all traces (stored in a list) to a Stream
stream_ray_trace = Stream()
for trace in result_list:
    stream_ray_trace += trace

# Write piercing points in a file
mtz.write_files_4_piercing_points_and_raypaths(stream_ray_trace, sta, piercing_depth=535, plot=False)

# Migration
mObs = mtz.ccpm_3d(stream_ray_trace, m_params, output_file="/home/" + work_dir.split('/')[2] + "/Desktop/mtz_example/example", phase="PS")


total_time = time.time() - t_beg
print('Ray tracing took ' + str(round(total_time)/60) + ' minutes in total.')
