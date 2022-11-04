===================================================
Receiver functions and time-to-depth migration calculations
===================================================
``rfmpy`` is a set of Python scripts for 1)
calculating receiver functions (RF) and 2) performing time to depth, in a 3D spherical coordinate system.

The codes are actively developed on `GitHub <https://github.com/kemichai/rfmpy/>`_ and were used for
producing the results for the manuscript draft entitled **"Moho depths beneath the European Alps:a homogeneously processed map and receiver functions database"**
to be submitted to *Earth System Science Data (ESSD)*.


--------------

Installation
~~~~~~~~~~~~

``rfmpy`` is currently at development and are subject to change
at any time. Also please note that, at least at this stage,
the codes are designed to reproduce our results.
For different applications the codes will need to be modified. The codes are only tested on Linux OS.

It is recommended to install ``rfmpy`` inside a Conda environment to
preserve your root environment. You can download Conda at the
following `link <https://docs.conda.io/en/latest/miniconda.html>`__.


.. code:: bash

   $ conda config --add channels conda-forge
   $ conda create -n rfmpy python=3.6 pip obspy=1.2.1 matplotlib numpy pandas basemap cartopy shapely fortran-compiler
   $ conda activate rfmpy
   $ conda install -c anaconda ipython=7.13
   $ git clone https://github.com/kemichai/rfmpy.git
   $ cd rfmpy
   $ pip install .



How to Cite
~~~~~~~~~~~~
If you use ``rfmpy``, consider citing the related publication:
`Michailos et al. (in prep.) <https://doi.org/10.5194/egusphere-egu22-8174>`__.

Michailos, K., Scarponi, M., Stipčević, J., Hetényi, G.,
Hannemann, K., Kalmár, D., Mroczek, S., Paul, A., Plomerová, P.,
Tilmann, F., Vergne, J., and AlpArray Working Group:
Moho depths beneath the European Alps from
receiver functions of the AlpArray Seismic Network,
EGU General Assembly 2022, Vienna, Austria, 23–27 May
2022, EGU22-8174, https://doi.org/10.5194/egusphere-egu22-8174, 2022.


.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: How to install

   installation

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Tutorial

   tutorial
