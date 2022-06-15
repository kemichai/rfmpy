<!---
LINK TO DOWNLOAD A DATASET TO USE...
LOOK here for more ideas:...https://github.com/insarlab/MintPy/tree/main/docs
--->
1) upload all XT data on ZENODO

```bash
wget https://zenodo.org/record/3952953/files/FernandinaSenDT128.tar.xz
tar -xvJf FernandinaSenDT128.tar.xz
cd FernandinaSenDT128/mintpy
smallbaselineApp.py ${MINTPY_HOME}/mintpy/data/input_files/FernandinaSenDT128.txt

2) store data on you pc

3) copy paste the RF calculation codes here 

4) add plots and figures for rfs

5) continue with migration and etc...
```

[Link for EPcrust](http://eurorem.bo.ingv.it/EPcrust_solar/)


#### 2.1 Routine workflow `compute_RF.py` ####

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


```

![My Image](../plots/rf_steps.jpg)

