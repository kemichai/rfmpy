Moho discontinuity example
==============

Introduction
~~~~~~~~~~~~
Here you can find a tutorial for calculating receiver functions and time-to-depth
calculation for a given subset of seismic waveform data from the EASI seismic network. Note
that this is only a small sample of all the available data and its only purpose
is to show the functionality of the codes. The example dataset is available in the
following `link <https://zenodo.org/record/7695125#.YxtWIdJByut>`__.

Here we start from continuous data cut around the arrival times of selected teleseismic events
and apply a systematic processing routine (see Figure 1 for details on the steps we follow).

.. figure:: images/RF_Migration_workflow.png
    :alt: Processing steps for Receiver Function and time-to-depth migration calculations.

    Processing steps for Receiver Function and time-to-depth migration calculations within ``rfmpy``.




Download example dataset
~~~~~~~~~~~~
First we need to download the seismic waveform data from a ZENODO
repository in our local computer.

1. Create a directory to store the waveform data:

.. code-block:: bash

   $ mkdir ~/Desktop/data_sample


2. Download the data sample from ZENODO in that directory along with two files (plot_EASI.sh,vk.cpt) that we will need later:

.. code-block:: bash

   $ wget https://zenodo.org/record/7292588/files/seismic_data.tar.xz -P ~/Desktop/data_sample/

.. parsed-literal::

    [2022-09-27 15:56:54]  https://zenodo.org/record/7292588/files/seismic_data.tar.xz
    Resolving zenodo.org (zenodo.org)... 188.184.117.155
    Connecting to zenodo.org (zenodo.org)|188.184.117.155|:443... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 141181064 (135M) [application/octet-stream]
    Saving to: ‘~/Desktop/data_sample/seismic_data.tar.xz’
    seismic_data.tar.xz.2             100%[==========================================================>] 134.64M  8.43MB/s    in 13s
    [2022-09-27 15:57:08] (10.2 MB/s) - ‘~/Desktop/data_sample/seismic_data.tar.xz’ saved [141181064/141181064]


.. code-block:: bash

   $ wget https://zenodo.org/record/7292588/files/plot_EASI.sh -P ~/Desktop/data_sample/
   $ wget https://zenodo.org/record/7292588/files/vik.cpt -P ~/Desktop/data_sample/


.. parsed-literal::

    Resolving zenodo.org (zenodo.org)... 188.184.117.155
    Connecting to zenodo.org (zenodo.org)|188.184.117.155|:443... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 602 [application/octet-stream]
    Saving to: ‘~/Desktop/data_sample/plot_EASI.sh.1’
    Saving to: ‘~/Desktop/data_sample/vik.cpt.1’
    plot_EASI.sh.1      100%[===================>]     602  --.-KB/s    in 0s
    vik.cpt.1           100%[===================>]   9.86K  --.-KB/s    in 0s




3. Extract files from the tar file we just downloaded:

.. code-block:: bash

   $ tar -xf ~/Desktop/data_sample/seismic_data.tar.xz --directory ~/Desktop/data_sample

4. Create a directory to store RFs:

.. code-block:: bash

    $ mkdir ~/Desktop/data_sample/RF
    $ mkdir ~/Desktop/data_sample/TRF


Calculate receiver functions
~~~~~~~~~~~~

Run the following, code snippet from the repository's top folder to compute receiver functions.


.. code-block:: python

    import rfmpy.core.RF_Main as RF
    from obspy import read_inventory, read_events, UTCDateTime as UTC
    import os
    import time

    # Define working directory
    work_dir = os.getcwd()

    # Path in which waveforms are stored
    path_wavs = ['/home/' + work_dir.split('/')[2] + '/Desktop/data_sample/EASI/data/']

    # Define path to store RFs
    path_out_RF = '/home/' + work_dir.split('/')[2] + '/Desktop/data_sample/'

    # Start a timer to keep a track how long the calculations take
    t_beg = time.time()

    # Path for StationXML files
    path_meta = work_dir + '/data/metadata/'
    try:
        print('>>> Reading inventory...')
        inv = read_inventory(path_meta + '/*.xml')
        print('>>> Read inventory...')
    except Exception as e:
        raise type(e)('>>> Move to the top directory of the repository!')

    # =================================================== #
    # Define parameters for calculating receiver functions
    # Define sta/lta parameters
    sta_lta_qc_parameters = {'sta': 3, 'lta': 50, 'high_cut': 1.0, 'threshold': 2.5}

    # Define pre-processing parameters
    pre_processing_parameters = {'low_cut': 0.05, 'high_cut': 1.0, 'order': 2,
                                 't_before': 40, 't_after': 60}
    for path_wav in path_wavs:
        print(path_wav)
        RF.calculate_rf(path_ev=path_wav, path_out=path_out_RF, inventory=inv, iterations=200,
                        ds=30, c1=10, c2=10, sta_lta_qc=sta_lta_qc_parameters,
                        pre_processing=pre_processing_parameters, max_frequency=1, save=True, plot=False)
    # ==================================================== #
    t_end = time.time()
    total_time = t_end - t_beg
    print('It took ' + str(round(total_time)/60) + ' minutes in total.')


.. parsed-literal::

    [2022-09-27 15:58:01] >>> Reading inventory...
    >>> Read inventory...
    /home/*/Desktop/data_sample/EASI/data/
    Calculating RF for event in: /home/*/Desktop/data_sample/EASI/data/P_2014.363.09.29.37
    ...
    >>> Station: XT.AAE50 - Failed on QC 2.
    [2022-09-27 16:57:08] It took 20 minutes in total.


This created 273 RF files in SAC format...


Calculate time-to-depth migration
~~~~~~~~~~~~
Now to compute time-to-depth migration for these RF traces we use the following
code snippet.


.. code-block:: python

    import rfmpy.core.migration_sphr as rf_mig
    import rfmpy.utils.migration_plots_spher as plot_migration_sphr
    import os
    import time

    # Start a timer to keep a track how long the calculations take
    t_beg = time.time()

    # Define working directory
    work_dir = os.getcwd()

    # Define path to RFs
    path = '/home/' + work_dir.split('/')[2] + '/Desktop/data_sample/RF/'

    # Read station coordinates from the rfs (sac files) in a pandas dataframe
    sta = rf_mig.read_stations_from_sac(path2rfs=path)

    # Read RFs
    stream = rf_mig.read_traces_sphr(path2rfs=path, sta=sta)

    # =================================================== #
    # Define MIGRATION parameters
    # Ray-tracing parameters
    inc = 0.25
    zmax = 100
    # Determine study area (x -> perpendicular to the profile)
    minx = 0.0
    maxx = 30.0
    pasx = 0.05
    miny = 30.0
    maxy = 60.0
    pasy = 0.05
    minz = -5
    # maxz needs to be >= zmax
    maxz = 100
    pasz = 0.5
    # Pass all the migration parameters in a dictionary to use them in functions called below
    m_params = {'minx': minx, 'maxx': maxx,
                'pasx': pasx, 'pasy': pasy, 'miny': miny, 'maxy': maxy,
                'minz': minz, 'maxz': maxz, 'pasz': pasz, 'inc': inc, 'zmax': zmax}


    # Ray tracing
    # Pick one of the two velocity models
    # 'EPcrust' or 'iasp91'
    # We use EPcrust velocity model here...
    stream_ray_trace = rf_mig.tracing_3D_sphr(stream=stream, migration_param_dict=m_params,
                                              velocity_model='EPcrust')

    # Write piercing points in a file
    plot_migration_sphr.write_files_4_piercing_points_and_raypaths(stream_ray_trace, sta, piercing_depth=35, plot=True)
    # Migration
    mObs = rf_mig.ccpm_3d(stream_ray_trace, m_params, output_file="/home/" + work_dir.split('/')[2] + "/Desktop/data_sample/epcrust", phase="PS")
    total_time = time.time() - t_beg
    print('Time-to-depth migration took ' + str(round(total_time)/60) + ' minutes in total.')



.. parsed-literal::

    |-----------------------------------------------|
    | Reading receiver functions...                 |
    | Reading trace 0 of 273
    ...
    | 273 of 273
    | End of 3D ray tracing...                      |
    |-----------------------------------------------|


.. figure:: images/piercing_points.png
    :alt: Map showing the piercing points (gray crosses)
          at 35 km depth computed for each seismic station
          (inverted red triangles) using the EPcrust velocity model (Molinari and Morelli, 2011).

    Map showing the piercing points (gray crosses)
    at 35 km depth computed for each seismic station (inverted red triangles) using the EPcrust velocity model (Molinari and Morelli, 2011).

.. parsed-literal::

    |-----------------------------------------------|
    | Start of common conversion point stacking...  |
    | 1 of 273
    ...
    | 273 of 273
    | End of common conversion point stacking...    |
    |-----------------------------------------------|
    Time-to-depth migration took 0.7 minutes in total.

This provides us with a 3D grid (epcrust.npy) of stacked migrated RF amplitudes.

Plot migrated cross-sections
~~~~~~~~~~~~
We will use this 3D grid to plot the cross-section using GMT6.
Before we do this, we need to create the cross-section

.. code-block:: python

    import rfmpy.core.migration_sphr as rf_mig
    import rfmpy.utils.migration_plots_spher as plot_migration_sphr
    import numpy as np
    import os
    import matplotlib.pyplot as plt
    from obspy.geodetics import degrees2kilometers, kilometers2degrees

    path2grid = '/home/' + work_dir.split('/')[2] + '/Desktop/data_sample/'

    # Read the 3D grid (epcrust.npy) of stacked migrated RF amplitudes.
    with open(path2grid + 'epcrust.npy', 'rb') as f:
        mObs_ep = np.load(f)

    profile = np.array([[13.35, 50.6], [13.35, 45.6]])
    profile_name = 'EASI'
    G2_, sta, xx, zz = plot_migration_sphr.create_2d_profile(mObs_ep, m_params, profile, sta, swath=37.5, plot=True)
    mObs = rf_mig.ccp_smooth(G2_, m_params)
    mObs = rf_mig.ccpFilter(mObs)


    # File for GMT plot
    for i, x in enumerate(xx):
        for j, z in enumerate(zz):
            print(kilometers2degrees(x), z, mObs[i,j])
            with open(path2grid + profile_name + '.txt', 'a') as of:
                of.write('{} {} {} \n'.
                         format(kilometers2degrees(x), z, mObs[i, j]))


Using the following commands we can create the cross-section using the GMT6 code we downloaded.

.. code-block:: bash

    $ cd ~/Desktop/data_sample/
    $ conda deactivate
    $ conda activate gmt6
    $ bash plot_EASI.sh


.. figure:: images/easi.png
    :alt: Migrated receiver-function cross-sections along the EASI seismic network.

    Migrated receiver-function cross-sections along the EASI seismic network.

.. warning::

    The image created here is based only on a small sample of the dataset available.

Mantle transition zone example
==============

Introduction
~~~~~~~~~~~~
Here you can find a tutorial for calculating time-to-depth migration for a given subset of seismic waveform data
from the **YYY** seismic network in order to map the thickness and phase boundaries of the mantle transition zone
(e.g., 410 km, 520 km, and 660 km discontinuities).

.. note::
    For the MTZ example we start the tutorial from having calculated the receiver functions (for how to start from
    raw waveform data have a look at the previous example).

Download example receiver function dataset
~~~~~~~~~~~~
First we have to download a subset of receiver functions from a ZENODO
repository in our local computer.

1. Create a directory to store the waveform data:

.. code-block:: bash

   $ mkdir ~/Desktop/mtz_example


2. Download the receiver functions from ZENODO in that directory:

.. code-block:: bash

   $ wget https://zenodo.org/records/14286133/files/RF_data.tar -P ~/Desktop/data_sample/

.. parsed-literal::

    [2024-12-06 10:59:41]  https://zenodo.org/records/14286133/files/RF_data.tar
    Resolving zenodo.org (zenodo.org)... 188.185.45.92, 188.185.43.25, 188.185.48.194, ...
    Connecting to zenodo.org (zenodo.org)|188.185.45.92|:443... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 111656960 (106M) [application/octet-stream]
    Saving to: ‘~/Desktop/mtz_example/RF_data.tar’
    RF_data.tar     100%[=================================================================>] 106.48M  3.63MB/s    in 32s
    [2024-12-06 11:00:14] (3.38 MB/s) - ‘~/Desktop/mtz_example/RF_data.tar’ saved [111656960/111656960]


3. Create a directory to store RFs:

.. code-block:: bash

    $ mkdir ~/Desktop/mtz_example/RF


3. Extract files from the tar file we just downloaded:

.. code-block:: bash

   $ tar -xf ~/Desktop/mtz_example/RF_data.tar --directory ~/Desktop/mtz_example/RF


Calculate time-to-depth migration
~~~~~~~~~~~~
To compute time-to-depth migration for these RF traces we use the following
code snippet. Grab a coffee or tea while running this as it will take a few minutes to run.

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
    mObs = mtz.ccpm_3d(stream_ray_trace, m_params, output_file="/home/" + work_dir.split('/')[2] + "/Desktop/mtz_example", phase="PS")


    total_time = time.time() - t_beg
    print('Ray tracing took ' + str(round(total_time)/60) + ' minutes in total.')

