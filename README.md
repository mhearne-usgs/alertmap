Introduction
------------

alertmap is (currently) a script for generating a map of probable earthquake early warning alert times, 
roughly based off of "Probabilistic Warning Times for Earthquake Ground Shaking in the San Francisco Bay Area", Allen, R.M., Seismological Research Letters Vol. 77, Number 3 pp. 371-376.

This code uses a local ShakeMap installation and assumes that the ShakeMap root path is "/home/shake/ShakeMap".


Installation and Dependencies
-----------------------------

This package depends on:
 * numpy, the fundamental package for scientific computing with Python. <a href="http://www.numpy.org/">http://www.numpy.org/</a>  
 * matplotlib, a Python 2D plotting library which produces publication quality figures. <a href="<a href="http://matplotlib.org/index.html">http://matplotlib.org/index.html</a>
 * scipy, a Python library which provides many user-friendly and efficient numerical routines such as routines for numerical integration and optimization. <a href="<a href="http://www.scipy.org/scipylib/index.html">http://www.scipy.org/scipylib/index.html</a>
 * obspy, a Python library for dealing with seismology data.
 * neicio, a Python library for reading/writing various spatial data formats (including ShakeMap grid.xml). 
 * rasterio, a Python library for reading/writing various spatial data formats.

The best way to install numpy,matplotlib,and scipy is to use one of the Python distributions described here:

<a href="http://www.scipy.org/install.html">http://www.scipy.org/install.html</a>

Anaconda and Enthought distributions have been successfully tested with smtools.

To install obspy:

pip install obspy

To install neicio:

pip install git+git://github.com/usgs/neicio.git

Installing rasterio can be complicated, depending on your setup.  See this page:
https://github.com/mapbox/rasterio

or, if using Anaconda, you can install rasterio using the following command at the time of this writing:

conda install -c https://conda.binstar.org/rsignell rasterio


To install alertmap:

pip install git+git://github.com/mhearne-usgs/alertmap.git

Uninstalling and Updating
-------------------------

To uninstall:

pip uninstall alertmap

To update:

pip install -U git+git://github.com/mhearne-usgs/alertmap.git

Command line usage
------------------

<pre>
usage: alertmap.py [-h] event

This script does the following:
    1) Find ShakeMap event folder from input ID.
    2) Parse alert.conf file in event data folder (i.e. /home/shake/ShakeMap/data/eventID/alert.conf)
    [FAULTS]
    lats = 32.1 32.2 32.3
    lons = -118.1 -118.2 -118.3
    [MAP]
    mmithresh = 6.0 #MMI value below which alert time will not be saved
    gmpe = AkkarBommer07
    gmice = WGRW11
    ipe = AW07_CEUS
    xmin = -119.0
    xmax = -117.0
    ymin = 31.0
    ymax = 33.0
    dx = 0.01
    dy = 0.01
    3) Replace grind.conf file in event data folder with one created using information supplied above
    4) For each epicenter in (lats,lons):
      a) Write new event.xml file
      b) Run ShakeMap grind program
      c) Find 4 nearest stations to epicenter, calculate P arrival times for each, return the slowest.
      d) Calculate S arrival times for each cell in ShakeMap, subtract SlowP+8.5 from each.
      e) Save grid of S arrival times
    5) TODO SCIENCE
    

positional arguments:
  event       Select the event ID

optional arguments:
  -h, --help  show this help message and exit
</pre>
