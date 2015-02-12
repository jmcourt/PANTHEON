#PANTHEON
#####Python ANalysis Tools for High-energy Event data manipulatiON
#####-J.M.Court, 2015

The new home of the former _XTE L-Extract_ project!

_**PANTHEON**_ is an open-source suite of tools to help analyse and visualise raw X-Ray satellite data.  Comes with tools to sort relevant data out of pre-downloaded archives, construct a variety of products and manipulate these.  It's easy to use, with most instructions given in the form of on-screen prompts, and the various tools are able to function mostly independently of each other.

##Specifications

* _SpecAngel_, _PlotDemon_ and _XTEPan\_Lib_ must be placed in the same directory as their associated library _Pan\_Lib_.
* _XTE-Get_ must be placed in the same directory as the associated XTE library _XTEPan\_Lib_.
* The following modules for Python are required:
  * astropy
  * cPickle
  * math
  * numpy
  * os
  * pylab
  * scipy
  * sys
  * warnings
* You must have the 'fkeyprint' and 'make\_se' tools from [HEASARC](http://heasarc.gsfc.nasa.gov/ftools/) in order to use _Xte Good-Xenon Extractor_ product.  If you do not intend to work with XTE GoodXenon files, or want to select the files manually, no HEASARC products are required.

##Mission-Independent Tools

###SpecAngel

_specangel.py_

See in-programme help for detailed instructions on specific commands.  Type 'help' to bring up a menu.

#####Functions:
* Opens .speca files and allows visualisation of power density spectra of data.
* Can create the average power density spectrum of an observation that has been split into chunks.
* Can display the individual power density spectrum for any given chunk.
* Creates spectrograms to show how the power density spectrum varies in time.
* Allows data to be clipped or rebinned.
* Allows manipulation of the scale in the spectrogram to more easily see peaks and troughs.

###PlotDemon

_plotdemon.py_

See in-programme help for detailed instructions on specific commands.  Type 'help' to bring up a menu.

#####Functions:
* Opens .plotd files and allows creation of a large number of plots.
* Allows plotting and manipulation of lightcurves, colour-colour diagrams and hardness-intensity diagrams.
* Lets users plot lightcurves of multiple different energy bands alongside each other to see how they correlate.

###Pan_Lib

_pan\_lib.py_

Contains useful mission-independant functions used in the other scripts.  See documentation in Pan_Lib for details on individual functions.

##RXTE PCA Tools

###XTE-Get

_xteget.py_

The all-in-one extractor for RXTE event files!  Currently implemented for _GoodXenon\_2s_ and _E\_125us\_64M\_0\_1s_ data types.

#####Functions:

* Extracts and analyses data from _RXTE_ event files.
* Allows user to select which energy channels they want to analyse.
* Splits the data into chunks ('Fourier Resolution') and creates power spectra of each chunk such that periodicity can be detected.
* Saves lightcurve information in a .plotd file to be read by the _**PANTHEON**_ programme _PlotDemon_.
* Saves power spectrum information in a .speca file to be read by the _**PANTHEON**_ programme _SpecAngel_.

###XTE Event-Extractor

_xteevex.sh_

When run in a directory, locates all event files saved in a .evt format and copies these into a new subdirectory named event0.

###XTE GoodXenon-Extractor

_xtegxex.sh_

When run in a directory, locates both parts of all possible Good Xenon files, and then produces these combined Good Xenon files.  Copies these resultant files into a new subdirectory named gx0.

###XTE-Pan_Lib

_xtepan\_lib.py_

Contains useful XTE-specific functions used in the other scripts.  See documentation in XTEPan_Lib for details on individual functions.

##_Coming soon_:

###FoldGenie

* A new section of PlotDemon that allows the folding of data.
* Takes user input for period, and seeks a better fit.
* Allows for plotting of light curves.

_Have fun using **PANTHEON**!_
