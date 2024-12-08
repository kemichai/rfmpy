Mantle transition zone example
==============

Introduction
~~~~~~~~~~~~
This tutorial demonstrates how to calculate time-to-depth migration for 100 randomly selected receiver functions,
enabling the mapping of the mantle transition zone (e.g., the 410 km, 520 km, and 660 km discontinuities).

.. note::
    For the MTZ example, the tutorial begins with pre-calculated receiver functions.
    If you need guidance on starting from raw waveform data, refer to the previous example.

Download example receiver function dataset
~~~~~~~~~~~~
First we have to download a subset of receiver functions from a ZENODO
repository in our local computer.

1. Create a directory to store the waveform data:

.. code-block:: bash

   $ mkdir ~/Desktop/mtz_example


2. Download the receiver functions from ZENODO in that directory along with two files (plot_cross_section.sh,vk.cpt) that we will need later:

.. code-block:: bash

   $ wget https://zenodo.org/records/14286724/files/100_RF_traces.tar.xz -P ~/Desktop/mtz_example/

.. parsed-literal::

    [2024-12-06 16:52:40]  https://zenodo.org/records/14286724/files/100_RF_traces.tar.xz
    Resolving zenodo.org (zenodo.org)... 188.185.43.25, 188.185.45.92, 188.185.48.194, ...
    Connecting to zenodo.org (zenodo.org)|188.185.43.25|:443... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 680364 (664K) [application/octet-stream]
    Saving to: ‘~/Desktop/mtz_example/100_RF_traces.tar.xz’
    100_RF_traces.tar.xz    100%[==============================>] 664.42K   446KB/s    in 1.5s
    [2024-12-06 16:52:43] (446 KB/s) - ‘~/Desktop/mtz_example/100_RF_traces.tar.xz’ saved [680364/680364]

.. code-block:: bash

   $ wget https://zenodo.org/records/14286724/files/plot_cross_section.sh -P ~/Desktop/mtz_example/
   $ wget https://zenodo.org/records/14286724/files/vik.cpt -P ~/Desktop/mtz_example/


3. Create a directory to store RFs:

.. code-block:: bash

    $ mkdir ~/Desktop/mtz_example/RF


3. Extract files from the tar file we just downloaded:

.. code-block:: bash

   $ tar -xf ~/Desktop/mtz_example/100_RF_traces.tar.xz --directory ~/Desktop/mtz_example/RF


Calculate time-to-depth migration
~~~~~~~~~~~~
To compute the time-to-depth migration for these RF traces, use the following code snippet.

.. note::
    This process can take up to 20 minutes to complete—plenty of time to enjoy a coffee or tea while you wait!


.. code-block:: python

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
    minx = -13.0 # degrees
    maxx = 46.0 # degrees
    pasx = 0.26 # degrees
    miny = 30.0 # degrees
    maxy = 64.0 # degrees
    pasy = 0.18 # degrees
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


.. parsed-literal::


    |-----------------------------------------------|
    | Reading receiver functions...                 |
    | Reading trace 0 of 100
    | Reading trace 1 of 100
    | Reading trace 2 of 100
    | Reading trace 3 of 100
    | Reading trace 4 of 100
    ...

    | 100 of 100
    | End of 3D ray tracing...                      |
    |-----------------------------------------------|

    |-----------------------------------------------|
    | Start of common conversion point stacking...  |
    | 1 of 100
    ...
    | 98 of 100
    | 99 of 100
    | 100 of 100
    | End of common conversion point stacking...    |
    |-----------------------------------------------|
    Ray tracing took 19.25 minutes in total.


This provides us with a 3D grid (example.npy) of stacked migrated RF amplitudes.


Plot migrated cross-sections
~~~~~~~~~~~~
We will use this 3D grid to plot the cross-section using GMT6.
Before we do this, we need to create the cross-section.

.. code-block:: python

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

    # Define MIGRATION parameters
    # Ray-tracing parameters
    inc = 2  # km
    zmax = 800 # km
    # Determine study area (x -> perpendicular to the profile)
    minx = -13.0 # degrees
    maxx = 46.0 # degrees
    pasx = 0.26 # degrees
    miny = 30.0 # degrees
    maxy = 64.0 # degrees
    pasy = 0.18 # degrees
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


Using the following commands we can create the cross-section using the GMT6 code we downloaded.

.. code-block:: bash

    $ cd ~/Desktop/mtz_example/
    $ conda deactivate
    $ conda activate gmt6
    $ bash plot_cross_section.sh


.. figure:: images/test.png
    :alt: Example of migrated receiver-function cross-section.

    Example of migrated receiver-function cross-section.

.. warning::

    The image generated here is based on a small sample of the dataset.
    This tutorial showcases the functionality of the codes without reproducing
    the full figures, which would require significant processing time.

