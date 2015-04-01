#PANTHEON
_**P**ython **An**alysis **T**ools for **H**igh-energy **E**vent data Manipulati**on**_
#####J.M.Court, 2015

_The new home of the former_ XTE L-Extract _project!_

_**PANTHEON**_ is an open-source suite of tools to help analyse and visualise raw X-Ray satellite data.  Comes with tools to sort relevant data out of pre-downloaded archives, construct a variety of products and manipulate these.  It's easy to use, with most instructions given in the form of on-screen prompts, and the various tools are able to function mostly independently of each other.

##Contents

1. Specifications
2. Overview
 1. Mission-Independent Tools
 2. Instrument Specific Tools
3. PlotDemon Manual
4. SpecAngel Manual

_____

##1. Specifications

* _FITSGenie_, _SpecAngel_ and _PlotDemon_ must be placed in the same directory as their associated library _Pan\_Lib_.
* The auxilliary function library relevant to the mission of interest (i.e. _SzkPan\_Lib_ or _XTEPan\_Lib_) must be placed in the same directory as _FITSGenie_.
* The following modules for Python are required:
  * astropy
  * cPickle
  * math
  * numba
  * numpy
  * os
  * PyAstronomy (for data folding only)
  * pylab
  * scipy
  * sys
  * warnings
* You must have the 'fkeyprint' and 'make\_se' tools from [HEASARC](http://heasarc.gsfc.nasa.gov/ftools/)'s FTOOLS package in order to use _Xte GoodXenon-Extractor_.  If you do not intend to work with XTE _GoodXenon_ files, or want to select the files manually, no HEASARC products are required.

_____
##2. Overview

##i. Mission-Independent Tools

###FITSGenie

>_fitsgenie.py_

The all-in-one extractor for event files!  Currently implemented for _GoodXenon\_2s_ and _E\_125us\_64M\_0\_1s_ data types from PCA on RXTE and .evt data types from XIS on SUZAKU (_Beta_).

#####Functions:

* Extracts and analyses data from FITS event files.
* Allows user to select which energy channels they want to analyse.
* Splits the data into chunks ('Fourier Windows') and creates power spectra of each chunk such that periodicity can be detected.
* Fourier Windows can be made discrete or sliding.
* Generates a number of statistics describing the data being procesed.
* Saves lightcurve information in a .plotd file to be read by the _**PANTHEON**_ programme _PlotDemon_.
* Saves power spectrum information in a .speca file to be read by the _**PANTHEON**_ programme _SpecAngel_.

#####Arguments:

Arguments can be given with the function call, else they will the user will be prompted to enter them by the software.  They are as follows:

* Minimum Channel: the lowest channel (or energy, depending on the data format input) of photons that will be used for analysis.
* Maximum Channel: the highest channel or energy to be used.
* Photon count bin-size: the size of bins to create the full-resolution histogram of photon events over time.
* Time per Fourier Spectrum: the size of the chunks the data will be divided into for Fourier analysis of saved data with _SpecAngel_.
* Estimate of background: an estimate of the background count rate in cts/s.

_____

###SpecAngel

>_specangel.py_

The power spectral analysis software for visualising time-domain variability!  See in-programme help for detailed instructions on specific commands.  Type 'help' to bring up a menu.

#####Functions:
* Opens .speca files and allows visualisation of power density spectra of data.
* Can create the average power density spectrum of an observation that has been split into chunks.
* Can display the individual power density spectrum for any given chunk.
* Creates spectrograms to show how the power density spectrum varies in time.
* Allows data to be clipped, renormalised or rebinned.
* Allows manipulation of the scale in the spectrogram to more easily see peaks and troughs.

#####Arguments:

Arguments can be given with the function call, else they will the user will be prompted to enter them by the software.  They are as follows:
* Logarithmic binning Factor: the value of _x_; two adjacent logarithmic bins in time will have their left-hand edges separated by no less than a multiplicative factor of 10^_x_.

_____

###PlotDemon

>_plotdemon.py_

The lightcurve and colour analysis software for hardness and intensity plotting!  See in-programme help for detailed instructions on specific commands.  Type 'help' to bring up a menu.

#####Functions:
* Opens up to 3 .plotd files and allows creation of a large number of plots comparing them.
* Allows plotting and manipulation of lightcurves, colour-colour diagrams and hardness-intensity diagrams.
* Lets users plot lightcurves of multiple different energy bands alongside each other to see how they correlate.
* Creates animations of how lightcurves change as the binning is increased.
* Allows folding of data over a period of the user's selection.
* Exports data to ASCII text files for compatibility with, for example, GnuPlot.

#####Arguments:

Arguments can be given with the function call, else they will the user will be prompted to enter them by the software.  They are as follows:
* Binsize: the size of the binning, in time, for lightcurves and colour curves.

_____

###DataFairy

>_datafairy.py_

Creates fake data, the form of which can be changed by the user, which is readable by _PlotDemon_.

#####Functions:
* The form of the data output is determined by functions in the _DataFairy_ code which can be edited by the user.
* Output data is in the form of three _PlotDemon_ .plotd files.
* Coming soon: output .speca files.

_____

###Pan_Lib

>_pan\_lib.py_

Contains useful mission-independent functions used in the other scripts.  See documentation in Pan_Lib for details on individual functions.

_____

##ii. Instrument-Specific Tools

###XTE Event-Extractor

>_xteevex.sh_

When run in a directory, locates all event files saved in a .evt format and copies these into a new subdirectory named event0.

_____

###XTE GoodXenon-Extractor

>_xtegxex.sh_

When run in a directory, locates both parts of all possible Good Xenon files, and then produces these combined Good Xenon files.  Copies these resultant files into a new subdirectory named gx0.

_____

###XTEPan_Lib

>_xtepan\_lib.py_

Contains useful XTE-specific functions used in the other scripts.  See documentation in XTEPan_Lib for details on individual functions.

_____

###SzkPan_Lib

>_szkpan\_lib.py_

Contains useful Suzaku-specific functions used in the other scripts.  See documentation in XTEPan_Lib for details on individual functions.

_____
##3. PlotDemon Manual

_Coming soon!_

_____
##4. SpecAngel Manual

_Coming soon!_

_____

_Have fun using **PANTHEON**!_
