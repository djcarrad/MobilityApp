# MobilityApp
Python-based tool to analyse conductance vs gate voltage data and extract mobility and density.

There are three options for installing, from easiest/most limited to hardest/greatest flexibility:

1. GUI only (Windows only)

If you are on windows and only want to use the GUI, the easiest thing is to download 'mobilityapp.exe' from the file list above.

2. GUI plus modules in python (all systems)

If need to embed the functions for calculating mobility within your own python code (e.g. automated analysis of large numbers of data sets), and/or prefer running the GUI from within your own python installation, simply use:

pip install mobilityapp

Calling mobilityapp from the command line (e.g. Anaconda prompt) runs the GUI

The 'examples' Jupyter notebook runs through how to use each module of the MobilityApp

3. Editable installation

To install an editable version using python installed from Anaconda:

In git bash:

cd C:/git

git clone https://github.com/djcarrad/MobilityApp

In Anaconda Prompt:

conda create -n mobility python

activate mobility

pip install -e C:/git/MobilityApp


The GUI can be run from Anaconda Prompt with simply:

mobilityapp