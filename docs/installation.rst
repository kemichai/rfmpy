## Installation:
A way to run these Python codes is with conda.

Before we get started, you need to install Anaconda.
Anaconda is cross-platform package manager software for scientific data analysis.
You can download the installation file based on your operating system and install Anaconda or
miniconda using the following [link](https://docs.conda.io/en/latest/miniconda.html)

Once you have installed conda, open a terminal (Linux)
create a new environment with the following dependencies using:
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


In case, you would like to run the gmt codes you can create a separate conda environment using the
commands bellow:
```bash
conda create --name gmt6
conda activate gmt6
conda config --prepend channels conda-forge
conda install python=3.9 gmt
```
