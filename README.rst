MobilityApp
===================================
Python-based tool to analyse conductance vs gate voltage data and extract mobility and density.

This application provides access to the methods described in this paper: EVENTUAL LINK

In short; for any FET-like device, conductance, G, vs gate voltage, Vg, data can be transformed into 
mobility vs density data, assuming that scattering at low density is dominated by random impurities.
The key is to find V0, the gate voltage at which the density tends towards zero; this process is at
the core of the app/code. In general, we do it by taking the derivative of G(Vg), finding the peak
co-ordinates, and using this information to extrapolate to an approximate V0. We then fit G(Vg) at
large Vg to find the series resistance in the circuit. Finally, you, the user need to input the 
geometric properties of the device, and the app will return mobility and density.

Questions/assistance: Damon Carrad, damonc@dtu.dk or Christian Petersen cenpe@dtu.dk


Installation:
===================================
There are three options for installing, from easiest/most limited to hardest/greatest flexibility:



1. GUI only (Windows only)
-----------------------------------
If you are on windows and only want to use the GUI, you can download 'mobilityapp.exe' from the file list above.



2. GUI plus modules in python (all systems*)
-----------------------------------
If you need to embed the functions for calculating mobility within your own python code 
(e.g. automated analysis of large numbers of data sets), and/or prefer running the GUI from 
within your own python installation, you can install from pip:

    pip install mobilityapp

To run the GUI, run from e.g. Anaconda Prompt:

    mobilityapp

To rescale/resize the GUI to better fit your screen, you can specify a multiplier, e.g.

    mobilityapp 2

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

Loading data
-----------------------------------
Currently only datasets compatible with numpy.loadtxt are supported; but these are most datasets.
Basically, if you can load it into Excel, you can probably load it into this app.

The app is flexible towards data being measured in siemens or units of the conductance quantum.
If you only have current data (recorded in ampere) this is also fine; you just need to provide the
source-drain voltage V_sd used.

Importantly! For the procedure to work accurately, the G(Vg) measurement should be performed to sufficently low 
voltages -- i.e. below the threshold voltage -- and also sufficiently high voltages -- such that an accurate
value of series resistance can be obtained.

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
just write the values manually below the plot. However, it's not clear what the uncertainties should 
be in this case, so if you manually declaring the peak co-ordinates, the uncertainties in mobility and density are not calculated.

Exported data
------------------------------------
The data can be exported either as a json database (most useful if you want to re-import to python),
or csv (most useful for re-import into other plotting programs). Since some of the data is processed
over a restricted range of V_g, the dataset is irregular (hence _not_ using numpy.savetxt)


Using the python code
====================================
The main advantages of using the jupyter notebook code are batch processing large numbers of datasets,
and troubleshooting. The notebook also allows you access to the full lmfit fit reports for each of the
fits, if further uncertainty reporting and processing is necessary. We hope the example notebook contains
enough comments to make usage clear enough, otherwise please reach out to us, damonc@dtu.dk or 
cenpe@dtu.dk

Definitions
=====================================
G: Conductance, in units of S, or the conductance quantum

Vg: Gate voltage

dG/dVg: Derivative of G with respect to Vg. In the app the derivative is performed numerically.

V0: The gate voltage for which the electron density tends to zero in the Drude model, assuming a linear
capacative coupling between the gate and FET channel. V0 = Vth - 2*(Vg_infl - Vth)

Vth: Threshold voltage. The gate voltage for which conductance tends to zero

Vg_infl: The gate voltage at the inflection point in G(Vg); equivalently, the maxima of dG/dVg

Rs: The series resistance. The total measured resistance minus the resistance of the device itself.
The series resistance includes contact resistance and any resistance introduced by the external circuitry.

Capacitance, C: The capcitance between the gate and the channel. Usually calculated or simulated.

Length/width: Dimensions of the device, underneath the gate (i.e. excluding any un-gated region)

Cap per area: Capacitance per area

density: calculated as Cperarea*(Vg-V0)/e where e is electron charge.

mu_eff: The effective mobility. Calculated as mu_eff=length^2/(C*(Vg-V0)*((1/G)-Rs)). If found correctly,
the effective mobility is identical to the Hall mobility in the single-carrier limit.

mu_FET: Field effect mobility found by fitting 1/G = Rs + length^2/(mu_FET*C*(Vg-Vth)). The field effect 
mobility may approximate the Hall/effective mobility for a small range of density, but is in general a 
less accurate and less quantitative measure of material quality.