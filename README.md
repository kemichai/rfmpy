## rfmpy - Receiver_Function_Migration_Py

### Description 
Python codes for calculating receiver functions (RF) and 
performing time to depth migration, in a 3D spherical coordinate system, to the RFs. 
We use these codes for RF calculations with the AlpArray Seismic Network. 

Original codes provided by Matteo Scarponi.
Modified by KM.

#### Requirements:
The simplest way to run the codes is using anaconda and a virtual environment.
This way you won't damage your system python.  
If you do not already have an anaconda or miniconda install go to the
conda-install instructions and get yourself either of the two [conda-install](https://docs.conda.io/en/latest/miniconda.html).

Once you have installed conda, create a new environment with the following dependencies using:
```bash
conda config --add channels conda-forge
conda create -n rfmpy python=3.6 pip obspy=1.2.1 matplotlib numpy pandas basemap cartopy shapely fortran-compiler
conda activate rfmpy
conda install -c anaconda ipython=7.13
```

Install from source:
```bash
git clone https://github.com/kemichai/rfmpy.git
cd rfmpy
```
Once you clone the project open a terminal in the
top directory (the one containing setup.py) and type the 
following to install the functions and make em
 available everywhere on your machine (within your environment).
```bash
pip install .
```


### Version
0.0.1

### Note
Codes are currently at development and are subject to 
change at any time.
 
### To do: ###
* Write tests []
* Code review/editing []
* Check that modified codes give same results as original codes MS sent [DONE]
* Rotate to Z-E-N before we rotate to R-T-Z []
* QC 1 parameters?
* Plot eq details in the deconv summary plots...
* Currently throwing away channels named 1, 2 (Need to orient them to Z-E-N and rename them)

