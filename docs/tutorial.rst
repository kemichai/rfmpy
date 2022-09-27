Tutorial
==============

Introduction
~~~~~~~~~~~~
Here you can find a tutorial for calculating receiver functions and time-to-depth
calculation for a given subset of seismic waveform data from the EASI seismic network. Note
that this is only a small sample of all the available data and its only purpose
is to show the functionality of the codes. The example dataset is available in the
following `link <https://zenodo.org/record/7065029#.YxtWIdJByut>`__.

Here we start from continuous data cut around the arrival times of selected teleseismic events
and apply a systematic processing routine (see Figure 1 for details on the steps we follow).

.. figure:: images/RF_Migration_workflow.png
    :alt: Processing steps for Receiver Function and time-to-depth migration calculations.

    Processing steps for Receiver Function and time-to-depth migration calculations within ``rfmpy``.




Download example dataset
~~~~~~~~~~~~
First we need to download the seismic waveform data from a ZENODO
repository in our local computer.

1. Create a directory to store the data:

.. code:: bash

   $ mkdir ~/Desktop/data_sample

2. Download the data sample from ZENODO in that directory:

.. code:: bash

   $ wget https://zenodo.org/record/7065029/files/seismic_data.tar.xz -P ~/Desktop/data_sample/

3. Extract files from the tar file we just downloaded:

.. code:: bash

   $ tar -xf ~/Desktop/data_sample/seismic_data.tar.xz --directory ~/Desktop/data_sample



Calculate receiver functions
~~~~~~~~~~~~

Run the following, code snippet to compute receiver functions.


.. code:: python3

   $ import rfmpy.core.RF_Main as RF
   $ import platform
   $ from obspy import read_inventory, read_events, UTCDateTime as UTC
   $ import os
   $ import time

   $ # Path in which waveforms are stored
   $ path_wavs = []

   $ # Define working directory
   $ work_dir = os.getcwd()

   $ # Define path to store RFs
   $ path_out_RF = work_dir + '/data/RF/'

   $ # Start a timer to keep a track how long the calculations take
   $ t_beg = time.time()

   $ # Path for StationXML files
   $ path_meta = work_dir + '/data/metadata/'
   $ try:
   $     print('>>> Reading inventory...')
   $     inv = read_inventory(path_meta + '/*.xml')
   $     print('>>> Read inventory...')
   $ except Exception as e:
   $     raise type(e)('>>> Move to the top directory of the repository!')

   $ # =================================================== #
   $ # Define parameters for calculating receiver functions
   $ # Define sta/lta parameters
   $ sta_lta_qc_parameters = {'sta': 3, 'lta': 50, 'high_cut': 1.0, 'threshold': 2.5}

   $ # Define pre-processing parameters
   $ pre_processing_parameters = {'low_cut': 0.05, 'high_cut': 1.0, 'order': 2,
   $                              't_before': 40, 't_after': 60}
   $ for path_wav in path_wavs:
   $     print(path_wav)
   $     RF.calculate_rf(path_ev=path_wav, path_out=path_out_RF, inventory=inv, iterations=200,
   $                     ds=30, c1=10, c2=10, sta_lta_qc=sta_lta_qc_parameters,
   $                     pre_processing=pre_processing_parameters, max_frequency=1, save=True, plot=False)
   $ # ==================================================== #
   $ t_end = time.time()
   $ total_time = t_end - t_beg
   $ print('It took ' + str(round(total_time)) + ' seconds in total.')


Calculate time-to-depth migration
~~~~~~~~~~~~
Run the following code snippet to compute time to depth migrations.


.. code:: python3

   $ import os


(WIP)
