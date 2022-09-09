## Tutorial
Here you can find a tutorial for calculating receiver functions and time-to-depth
calculation for a given subset of seismic waveform data (EASI seismic network). Note 
that this is only a small sample of all the available data and its only purpose 
is to show the functionality of the codes. The data are available in the 
following [link](https://zenodo.org/record/7065029#.YxtWIdJByut).

<!---
INK TO DOWNLOAD A DATASET TO USE...
LOOK here for more ideas:...https://github.com/insarlab/MintPy/tree/main/docs
-->



First we need to download the seismic waveform data from a ZENODO 
repository in our local computer. 

1. Create a directory to store the data:
    ```bash
    mkdir ~/Desktop/data_sample
    ```
2. Download the data sample from ZENODO in that directory:
    ```bash
    wget https://zenodo.org/record/7065029/files/seismic_data.tar.xz -P ~/Desktop/data_sample/
    ```
3. Extract files from the tar file we just downloaded:
    ```bash
    tar -xf ~/Desktop/data_sample/seismic_data.tar.xz --directory ~/Desktop/data_sample
    ```
   




[Link for EPcrust](http://eurorem.bo.ingv.it/EPcrust_solar/)

![My Image](images/RF_Migration_workflow.png)
_Figure 1: Processing steps for Receiver Function and time-to-depth migration calculations._


#### 2.1 Routine workflow `compute_RF.py` ####
Run the following, code snippet to compute receiver functions.

```python
import rfmpy.core.RF_Main as RF
import platform
from obspy import read_inventory, read_events, UTCDateTime as UTC
import os
import time

# Path in which waveforms are stored
path_wavs = []

# Define working directory
work_dir = os.getcwd()

# Define path to store RFs
path_out_RF = work_dir + '/data/RF/'

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

# =================================================================================================================== #
# Define parameters for calculating receiver functions
# Define sta/lta parameters
sta_lta_qc_parameters = {'sta': 3, 'lta': 50, 'high_cut': 1.0, 'threshold': 2.5}
# Define pre-processing parameters
pre_processing_parameters = {'low_cut': 0.05, 'high_cut': 1.0, 'order': 2, 't_before': 40, 't_after': 60}
for path_wav in path_wavs:
    print(path_wav)
    RF.calculate_rf(path_ev=path_wav, path_out=path_out_RF,
                inventory=inv, iterations=200, ds=30,
                c1=10, c2=10,
                sta_lta_qc=sta_lta_qc_parameters,
                pre_processing=pre_processing_parameters,
                max_frequency=1, save=True, plot=False)
# =================================================================================================================== #

t_end = time.time()
total_time = t_end - t_beg
print('It took ' + str(round(total_time)) + ' seconds in total.')

```


#### 2.1 Routine workflow `compute_RF_migration_spher.py` ####
Run the following code snippet to compute time to depth migrations.

```python3
import os

```

(WIP)
