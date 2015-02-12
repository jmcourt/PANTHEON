#PANTHEON

The new home of the former _XTE L-Extract_ project! Extracts data from astronomical observation 'event files' and allows this data to be manipulated and visualised. Current functionality of _PANTHEON_ is as such:

##XTE-Get

*Extracts and analyses data from _RXTE_ event files.
*Allows user to select which energy channels they want to analyse.
*Splits the data into chunks ('Fourier Resolution') and creates power spectra of each chunk such that periodicity can be detected.
*Saves lightcurve information in a .plotd file to be read by the _PANTHEON_ programme _PlotDemon_.
*Saves power spectrum information in a .speca file to be read by the _PANTHEON_ programme _SpecAngel_.
*Currently implemented for **GoodXenon\_2s** and **E\_125us\_64M\_0\_1s** data types.

##SpecAngel

*Opens .speca file and allows visualisation of power spectrum data.
*Can create the average power spectrum of an observation that has been split into chunks.
*Can display the individual power spectrum for any given chunk.
*Creates spectrograms to show how the power spectrum varies in time.
*Allows data to be clipped or rebinned.
*Allows manipulation of the scale in the spectrogram to more easily see peaks and troughs.

##PlotDemon

*PlotDemon from the XTE L-Extract project is being rewritten to work with .plotd files created by XTE-Get.
*Allows plotting and manipulation of lightcurves, colour-colour diagrams and hardness-intensity diagrams.
*Lets users plot lightcurves of multiple different energy bands alongside each other to see how they correlate.
*Folds data to a period of the user's choice to find the waveform of periodic data.
*Allows the user to search for a better fit period

##XTE L-Lib

The library from the old XTE L-Extract project. Contains useful functions used in the other scripts.

##_Coming soon_:


