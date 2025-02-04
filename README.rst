MobilityApp
===================================
Python-based tool to analyse conductance vs gate voltage data and extract mobility and density.

See below for background and instructions on use.

Questions/assistance: Damon Carrad, damonc@dtu.dk or Christian Petersen cenpe@dtu.dk


Installation:
===================================
There are three options for installing, from easiest/most limited to hardest/greatest flexibility:



1. GUI only (Windows only)
-----------------------------------
If you are on windows and only want to use the GUI, the easiest thing is to download 'mobilityapp.exe' from the file list above.



2. GUI plus modules in python (all systems*)
-----------------------------------
If you need to embed the functions for calculating mobility within your own python code (e.g. automated analysis of large numbers of data sets), and/or prefer running the GUI from within your own python installation, simply use:

    pip install mobilityapp

To run the GUI, simply run from the command line:

    mobilityapp

The 'examples' Jupyter notebook runs through how to use each module of the MobilityApp. If Jupyter isn't already installed:

    pip install jupyterlab



3. Editable installation (all systems*)
-----------------------------------
To install an editable version using python installed from Anaconda:

In git bash:

    cd C:/git

    git clone https://github.com/djcarrad/MobilityApp


In Anaconda Prompt:

    conda create -n mobility python

    activate mobility

    pip install -e C:/git/MobilityApp

    pip install jupyterlab

The GUI can be run from Anaconda Prompt with:

    mobilityapp

*Note: Only tested on Windows, but the code is pure python with a tkinter gui, so should work fine on Mac and 'nix.


Using the app:
===================================
The best tutorial for using the app/code is the paper itself, especially the instruction box.
The second best tutorial is to hover the mouse over the different elements of the app, and you
may receive helpful hints.

In short; for any FET-like device, conductance vs gate voltage data can be transformed into 
mobility vs density data, assuming that scattering at low density is dominated by random impurities.
The key is to find V_0, the gate voltage at which the density tends towards zero; this process is at
the core of the app/code. In general, we do it by taking the derivative of G(V_g), finding the peak
co-ordinates, and using this information to extrapolate to an approximate V_0. We then fit G(V_g) at
large V_g to find the series resistance in the circuit. Finally, you, the user need to input the 
geometric properties of the device, and the app will return mobility and density.

Loading data
-----------------------------------
Currently only datasets compatible with numpy.loadtxt are supported; but these are most datasets.
Basically, if you can load it into Excel, you can load it into this app. If you have a fancy database
or something like that, you will have to reshape it into something that can be loaded into Excel.

The app is flexible towards data being measured in siemens or units of the conductance quantum.
If you only have current data (recorded in ampere) this is also fine; you just need to provide the
source-drain voltage V_sd used.

Fitting the derivative
-----------------------------------
Any method involving fitting the peak of a numerical derivative is going to be sensitive to 
uncertainties/errors. Finding both the peak value of dG/dVg and the Vg at which this occurs is 
critical to the success of the method. To make the process less immune to noise or quantum fluctuations,
the app uses a specific form of an asymmetric Lorentzian that conforms well enough to the
derivative of most G(V_g) data we have tried. The important point is that only the peak
coordinates are relevant! Therefore, it doesn't matter too much if the fit is poor at other Vg; as 
long as the fitting procedure finds the co-ordinates accurately, it's all good. But please check it's 
done a good job!! If it hasn't done a good job, play with the parameters until it does. If you really
give up, but it's super obvious where the peak is (the human eye is a really good fitting tool), 
just write the values manually below the plot.

Exported data
------------------------------------
The data can be exported either as a json database (most useful if you want to re-import to python),
or csv (most useful for re-import into other plotting programs). Since some of the data is processed
over a restricted range of V_g, the dataset is irregular.


Using the python code
====================================
The main advantages of using the jupyter notebook code are batch processing large numbers of datasets,
and troubleshooting. The notebook also allows you access to the full lmfit fit reports for each of the
fits, if uncertainty reporting and processing is necessary. I hope the example notebook contains
enough comments to make usage clear enough, otherwise please reach out to me, damonc@dtu.dk.