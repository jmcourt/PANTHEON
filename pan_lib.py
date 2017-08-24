#! /usr/bin/env python

# |----------------------------------------------------------------------|
# |-------------------------------PAN_LIB--------------------------------|
# |----------------------------------------------------------------------|

# A selection of useful functions which are placed here to reduce clutter in the other files of
# PANTHEON.
#
# Contents:
#
#  ARGCHECK  - compares the list of arguments against a value given as the minimum allowed number of
#              of arguments.  If the list of arguments is too short, throw a warning and kill the script.
#
#  BINIFY    - takes a x-series with its associated y-axis data and y-axis errors.  Rebins the data
#              into larger linear bins with a width of the user's choosing, and returns the tuple.
#              x,y,y_error.
#
#  BOOLVAL   - takes a list of Boolean values and, interpreting it as binary, returns its integer value.
#
#  CIRCFOLD  - 'circularly folds' data by converting [time, data] points into polar-coordinate [phase, data]
#              pairs and returning these in cartesian or polar coordinates.
#
#  EQRANGE   -
#
#  EVAL_BURST-
#
#  FILENAMECHECK  - checks to see whether a proposed input file has the correct file extension.
#
#  FOLDIFY   - takes a time series with its associated y-axis data and y-axis errors.  Folds this data
#              over a time period of the user's choosing, and returns them as the tuple x,y,y_error.
#
#  FOLD_BURSTS - uses GET_BURSTS to obtain burst locations then interpolates to populate phase information
#              for all other points
#
#  GET_BURSTS- takes an array of data, looks for bursts and returns an array of tuples containing
#              the start and end points of these bursts.
#
#  GET_DIP   - returns the index of the lowest point between two user-defined flags in a dataset.
#
#  GTIMASK   - returns a data mask when given a time series and a GTI object
#
#  LBINIFY   - takes a linearly binned x-series with associated y-axis data and y-axis errors and rebins
#              them into bins of a constant width in logx space.  In places where the logarithmic bins
#              would be finer than the linear bins, the linear bins are retained.
#
#  LEAHYN    - takes the raw power spectrum output from the scipy FFT algorithm and normalises it using
#              Leahy normalisation.
#
#  LH2RMS    - takes a Leahy-normalised power spectrum and converts it to an (RMS/Mean)^2-normalised
#              power spectrum.
#
#  LHCONST   - returns the normalisation of the white noise component in a Leahy-normalised power spectrum
#              with no features in the range 1.5kHz - 4kHz.
#
#  MXREBIN   - takes a 2-dimensional set of data and corresponding errors linearly binned on the x-axis and
#              rebins them by an integer binning factor of the user's choice.
#
#  NONES     - like np.zeros, but with None.
#
#  PDCOLEX   - extracts colours from a set of 2 or 3 lightcurves
#
#  PLOTDLD   - load and unpickle a .plotd file and extract its data.
#
#  PLOTDSV   - collect a selection of data products as a library, pickle it and save as a .plotd file.
#
#  RMS_N     - takes the raw power spectrum output from the scipy FFT algorithm and normalises it using
#              (RMS/Mean)^2 normalisation.
#
#  SAFE_DIV  - Divides two arrays by each other, replacing NaNs that would be caused by div 0 errors with
#              zeroes.
#
#  SIGNOFF   - prints an dividing line with some space.  That's all it does.
#
#  SINFROMCOS- calculates the sines of an array of values when also passed their cosines.  If both sines
#              and cosines of the array are required, this method is faster than calling both trig functions.
#              Also contains function COSFROMSIN.
#
#  SLPLOT    - plots an x-y line plot of two sets of data, and then below plots the same data on another
#              set of axes in log-log space.
#
#  SPECALD   - load and unpickle a .speca file and extract its data.
#
#  SPECASV   - collect a selection of data products as a library, pickle it and save as a .speca file.
#
#  SRINR     - calculates whether a value given by a user is within an existant evenly-spaced array and,
#              if it is, returns the index value of the closest match to this value within the array.
#              Intended for validating subranges specified by user.
#
#  TNORM     - takes a list of times, and subtracts the lowest value from each entry such that a new
#              list starting with 0 is produced.  Large number subtraction errors are avoided by checking
#              that every entry is an integer number of time-resolution steps from zero.
#
#  UNIQFNAME - checks if a proposed filename is currently in use and, if so, proposes an alternative
#              filename to prevent overwrite.
#
#  XTRFILLOC - takes a filepath and outputs the file name and its absolute(ish) location
#
#

#-----Importing Modules------------------------------------------------------------------------------------------------

import os
try:
   import cPickle as pickle
except:
   import pickle
import pylab as pl
import warnings
import scipy.optimize as optm
import scipy.interpolate as intp
import scipy.signal as sgnl
import numpy as np
try:
   import numba as nb
   gotnumba=True
except:
   print 'Warning: numba module not found!  May run slow.'
   gotnumba=False
from math import pi


#-----Setup Conditional Wrapper to use jit when available--------------------------------------------------------------

class mjit(object):
    def __call__(self, f):
        if not gotnumba:
            return f
        else:
            return nb.jit(f)


#-----ArgCheck---------------------------------------------------------------------------------------------------------

def argcheck(x,y):

   '''Argument Checker

   Description:

    Takes a list (of arguments passed into a script), and a value representing the smallest allowable
    number of arguments.  If the length of the list is smaller than the criterion, an error is returned
    and the script is killed.

   Inputs:

    x - LIST   : The list of arguments passed into the current working script, i.e. x=sys.argv.
    y - INTEGER: The minimum allowable number of arguments.  Caution; if called from bash, sys.argv
                 returns the function call as the x[0] element, so y is 1 greater than the number
                 of user inputs.

   Outputs:

    [none]

   -J.M.Court, 2014'''

   if len(x)<y:
      print 'Not enough arguments!'
      signoff()
      exit()


#-----Binify-----------------------------------------------------------------------------------------------------------

@mjit()
def binify(x,y,ye,binsize):                                               # Defining 'binify' subscript

   '''Binify

   Description:

    Takes a 2-dimenstional set of data which has already been evenly binned on the x-axis, and re-bins
    it into larger, evenly spaced bins on the x-axis.

   Inputs:

    x       -  LIST: The x-values of the two-dimensional data.
    y       -  LIST: The y-values of the two-dimensional data, must be the same length as x.
    ye      -  LIST: The errors associated with the y-values of the two-dimensional data, must be the
                     same length as x and y.
    binsize - FLOAT: the size of the new x-axis bins in which to re-bin the data.

   Outputs:

    xb      -  LIST: The re-binned x-values of the two-dimensional data, i.e. an array of the left-hand
                     edges of the new bins.
    yb      -  LIST: The re-binned y-values of the two-dimensional data.
    yeb     -  LIST: The errors associated with the rebinned y-values of the two-dimensional data.

   -J.M.Court, 2014'''

   binlx=binsize*np.floor(x[0]/binsize)                                   # Initialising 'bin lowest x', or the lowest x value of the current bin
   binct=0.0                                                              # Initialising 'bin count', or number of values sorted into the current bin

   xb=[x[0]]                                                              # Setting up arrays to append binned values into
   yb=[0]   
   yeb=[0]

   for xid in range(len(x)):

      if x[xid]-binlx < binsize:                                          ## If the difference between the current x and bin start x is less than bin width:

         binct+=1
         yb[-1]=yb[-1]+y[xid]                                                  #  Add y to current bin
         yeb[-1]=yeb[-1]+(ye[xid]**2)                                          #  Add y error in quadrature to current bin
 
      else:                                                               ## Otherwise:
         binlx+=binsize                                                        #  Create new bin with minimum x equal to current x
         xb.append(binlx)                                                      #  Append next x value into new array element
         yb[-1]=yb[-1]/binct                                                   #  Divide y in previous bin by bincount to get the average
         yeb[-1]=(np.sqrt(yeb[-1]))/binct                                      #  Sqrt error and divide by bincount
         yb.append(y[xid])                                                     #  Append current y value into new array element
         yeb.append((ye[xid])**2)

         binct=1                                                               #  Reset bin count to 1

   yb[-1]=yb[-1]/binct                                                    ## Clean up final bin
   yeb[-1]=(np.sqrt(yeb[-1]))/binct
   return np.array(xb),np.array(yb),np.array(yeb)


#-----BoolVal----------------------------------------------------------------------------------------------------------

def boolval(data,reverse=True):

   '''Boolean Evaluator

   Description:

    Given a list of Boolean values, interprets them as binary and returns the corresponding integer.
    By default reads lists as having the highest-value digit first.

   Inputs:

    data    - LIST: a list of lists Boolean values.
    reverse - BOOL: If set to False, then the Boolean strings will be interpreted as binaries with the
                    lowest value (1) first.

   Outputs:

    data    - LIST: the integer values represented by 'data' if its values are interpreted as binary.

   -J.M.Court, 2015'''

   keyr=range(len(data[0]))                                               # Set up the 'key range' to convert Bool list into int

   if reverse:
      keyr=keyr[::-1]                                                     # Reverse the key range

   mult=[1<<i for i in keyr]                                              # Mult is a list of all the powers of 2 from 2^0 to 2^(length of data)
   data=np.array(data)*mult
   data=np.sum(data,axis=1)                                               # Multiply Boolean list by mult, sum per row

   return data
 

#-----Circfold---------------------------------------------------------------------------------------------------------

def circfold(x,y,t,pcoords=True):

   '''Circular folder

   Description:

    Circularly folds data by converting the x array into an array of phase angles between 0 and 2pi
    and the y array into an array of radial distances.  The coordinates of these points are then
    returned in Cartesian or Polar co-ordinates, as specified by the user.

   Inputs:

    x      - ARRAY: The x (time) co-ordinates of the data.
    y      - ARRAY: The y co-ordinates of the data.
    t      - FLOAT: The period over which the data is to be folded.
    pcoord -  BOOL: [Optional: Default=True] If set to true, the output co-ordinates will be theta, r
                    as in polar co-ordinates.  If set to false, the output co-ordinates will be x, y
                    as in Cartesian co-ordinates.

   Outputs:

    See inputs: pcoord

   -J.M.Court, 2015'''

   mult=(2*pi)/float(t)
   a=x*mult

   if pcoords:
      return a%(2.0*pi),y
   else:
      s=y*np.sin(a)
      c=y*np.cos(a)
      return s,c


#-----CCor-------------------------------------------------------------------------------------------------------------

# Cross-Correlate

def ccor(data1,data2):
   assert len(data1)==len(data2)
   data1=np.array(data1)
   data2=np.array(data2)
   crosscor=[]
   crosscore=[]
   for i in range(len(data1)):
      r1=data1[(-i-1):]
      r2=data2[:i+1]
      mr1=np.mean(r1)
      mr2=np.mean(r2)
      sr1=np.std(r1)
      sr2=np.std(r2)
      ccorarray=(r1-mr1)*(r2-mr2)/(sr1*sr2)
      crosscor.append(np.mean(ccorarray))
      crosscore.append(np.std(ccorarray))
   for i in range(len(data1)-1):
      r1=data1[:(-i-1)]
      r2=data2[i+1:]
      mr1=np.mean(r1)
      mr2=np.mean(r2)
      sr1=np.std(r1)
      sr2=np.std(r2)
      ccorarray=(r1-mr1)*(r2-mr2)/(sr1*sr2)
      crosscor.append(np.mean(ccorarray))
      crosscore.append(np.std(ccorarray))

   times=np.array(range(len(crosscor)))
   times=times+1-len(data1)
   return times,crosscor,crosscore
   

#-----Eddington--------------------------------------------------------------------------------------------------------

def eddington(M):

   '''Eddington

   Description:

    Returns the Eddington Limit for a given black hole mass, assuming Hydrogen accreta.

   Inputs:

    M - FLOAT: The mass of the black hole in solar masses

   Outputs:

    L - FLOAT: The Eddington Luminosity in ergs/s

   -J.M.Court, 2016'''

   return 1.26E38*M


#-----EqRange----------------------------------------------------------------------------------------------------------

# Equal Length Range

def eqrange(array):

   '''Equal Length Array

   Description:

    Creates a range with an equal length to that of the array or list given as an argument

   Inputs:

    array - ARRAYLIKE: The array to be used as comparison.

   Outputs:

    A range with equal length to the input

   -J.M.Court, 2015'''

   return range(len(array))

#-----Eval_Burst-------------------------------------------------------------------------------------------------------

# Evaluate Burst

@mjit()
def eval_burst(t,y):

   '''Evaluate Burst

   Description:

    Takes a piece of data and its associated time array, and treats it as a single 'burst'-like pattern.
    Returns the peak flux, peak time, rise time, fall time of this burst.  Works best with Get_Bursts.

   Inputs:

    t         - ARRAY: The t (time) co-ordinates of the data.
    y         - ARRAY: The y co-ordinates of the data.

   Outputs:

    peak      - FLOAT: The highest y-value in the burst
    peak_time - FLOAT: The time at which the peak occurs
    rise_time - FLOAT: The time between the start of the burst and its peak
    fall_time - FLOAT: The time between the peak of the burst and its end

   '''

   peak=max(y)
   trough=min(y)
   p_ind=np.array(y).tolist().index(peak)
   rise_time=t[p_ind]-t[0]
   fall_time=t[-1]-t[p_ind]
   peak_time=t[p_ind]
   return peak,trough,peak_time,rise_time,fall_time


#-----FlnCheck---------------------------------------------------------------------------------------------------------

def filenamecheck(filename,validext,continu=False):

   '''Filename Checker

   Description:

    Takes a filename and a string representing the expected file extension.  If the extension of the
    file does not match expectations, either kill the script or return 'False'.

   Inputs:

    filename - STRING: The filename to be checked.
    validext - STRING: The extenstion expected for the file (WITHOUT the leading '.')
    continu  -   BOOL: [Optional: Default=False] If True, returns a value of False for an incorrect
                       file extension.  If False (default), kills the script upon finding an
                       incorrect file extension.

   Outputs:

    iscorr   -   BOOL: True if the file has the correct extension, False otherwise

   -J.M.Court, 2015'''

   flext=(filename.split('.')[-1])

   if flext != validext:
      if continu:
         return False
      else:
         print 'Invalid input file!  Must use .'+str(validext)+' file!'
         signoff()
         exit()
   else:
      return True


#-----Foldify----------------------------------------------------------------------------------------------------------

def foldify(t,y,ye,period,binsize,phres=None,name='',compr=False,verb=True):

   '''Foldify

   Description:

    Folds a two-dimensional set of data over a period in the first dimension using the folding script
    which is provided in PyAstronomy.  Also calculates how much the peak-trough difference of the
    data has been compressed by the fold, and tells the user.

   Inputs:

    t       -   LIST: The t- or x-values of the two-dimensional data.  The period to be folded over
                      is a period in this dimension.
    y       -   LIST: The y-values of the two-dimensional data, must be the same length as x.
    ye      -   LIST: The errors associated with the y-values of the two-dimensional data, must be
                      the same length as x and y.
    period  -  FLOAT: The period of the data, in the same units as the data's x-values.
    binsize -  FLOAT: The size of the new x-axis bins in which to re-bin the data.
    phres   -  FLOAT: [Optional: Default=None] The resolution, in frequency space, at which data will
                      be shown.
    name    - STRING: [Optional: Default=''] the name of the file to be folded.  This name will be
                      used in text outputs printed to screen.
    compr   -   BOOL: [Optional: Default=False] if set True, gives a fourth output which is the amount
                      by which the min-max difference of the data-set was compressed by folding.
    verb    -   BOOL: [Optional: Default=True] if set false, suppresses all non-error text output.

   Outputs:

    newt    -  LIST: A clipped version of the input t which now only corresponds to one period.
    newy    -  LIST: The folded y-values of the two-dimensional data for one period.
    newye   -  LIST: The errors associated with the folded y-values of the two-dimensional data for
                     one period.
    compr   - FLOAT: see 'compr' input option

   -J.M.Court, 2015'''

   try:
      from PyAstronomy.pyasl import foldAt
   except:
      print 'PyAstronomy module not found; aborting!'
      return t,y,ye

   tn=tnorm(t,binsize)
   if verb:
      print 'Folding File '+name+'...'
   ptdiff=max(y)-min(y)                                                   # Flux range before folding

   phases=foldAt(tn,period)
 
   if phres==None:
      try:
         phres=float(raw_input('Input Phase Resolution (0-1): '))
         assert phres<=1.0
      except:
         print 'Invalid phase resolution!  Aborting!'
         return t,y,ye

   npbins=int(1.0/phres)
   phasx =np.arange(0,1,phres)
   phasy =np.zeros(npbins)
   phasye=np.zeros(npbins)
   ny=np.zeros(npbins)

   for i in range(len(y)):
      k=int(phases[i]*npbins)
      phasy[k]+=y[i]
      phasye[k]+=(ye[i]**2)
      ny[k]+=1

   phasy=phasy/ny
   phasye=np.sqrt(phasye)/ny

   afdiff=max(phasy)-min(phasy)                                            # Flux range after folding
   if verb:
      print 'Flattened by '+str(100-afdiff/ptdiff*100)+'%'

   if compr:
      return phasx,phasy,phasye,(afdiff/ptdiff)

   else:
      return phasx,phasy,phasye


#-----Fold_Bursts------------------------------------------------------------------------------------------------------

@mjit()
def fold_bursts(times,data,q_lo=50,q_hi=90,do_smooth=False,alg='cubic spline',savgol=5):

   '''Return Phases Using Bursts as Reference Points

   Description:

    Takes a lightcurve and identifies 'bursts' in the data.  The peak of each burst is considered
    to be at zero phase, and all other points are assigned phases by linearly interpolating between them.

   Inputs:

    data       - ARRAY: The data in which bursts are sought.
    q_low      - FLOAT: [Optional: Default=30] The percentile value of the data which will be used as
                        the low-pass threshold.  This threshold determines the edges of already-located
                        bursts, and thus changing it will change the quality of bursts but not the
                        quantity.
    q_mid      - FLOAT: [Optional: Default=50] The percentile value of the data which will be used as
                        the med-pass threshold.  This threshold determines how deep an trough must be
                        before the regions either side of it are considered separate bursts-candidate
                        regions.  Must be larger than q_low.
    q_hi       - FLOAT: [Optional: Default=70] The percentile value of the data which will be used as
                        the hi-pass threshold.  This threshold determines how much peak flux a previously
                        defined burst-candidate region must be before it is considered a true burst.
                        larger than q_med.
    smooth     -  BOOL: [Optional: Default=False] Apply a Savitsky-Golay Filter to smooth the lightcurve
                        before fetching peaks.
    alg        -STRING: [Optional: Default='cubic spline'] The algorithm to use to obtain peaks.


   Outputs:

    burst_locs -  LIST: A list of tuples containing the start and end indices of each burst.

   -J.M.Court, 2015'''

   assert len(times)==len(data)
   peaks=get_bursts(data,q_lo,q_hi,just_peaks=True,smooth=do_smooth,alg=alg,times=times,savgol=savgol)
   peaks.sort()
   #if peaks[0]!=0:
   #   peaks=np.hstack([0,peaks])
   #if peaks[-1]!=len(times)-1:
   #   peaks=np.hstack([peaks,len(times)-1])

   p0=peaks[0]
   pe=peaks[-1]

   #data=data[p0:pe+1]
   #times=times[p0:pe+1]
   phases=np.zeros(len(times))
   peaks=np.array(peaks)

   npeaks=[0]   

   for i in range(len(peaks)-1):
      if peaks[i+1]-peaks[i]>0.1*(peaks[-1]/float(len(peaks))):           # Remove any peaks separated by less than 10% of average separation
         npeaks.append(peaks[i+1])
   peaks=npeaks

   numpeaks=len(peaks)

   phases=get_phases_intp(data,windows=1,q_lo=20,q_hi=90,peaks=peaks,givespline=False)

   return np.array(phases),numpeaks,(peaks[0],peaks[-1])  


#-----Gauss------------------------------------------------------------------------------------------------------------

@mjit()
def gauss(mean,standev,x):

   '''Gauss.  Returns a Gaussian'''

   return (1.0/(standev*(2*np.pi)**0.5))*np.exp(-(x-mean)**2/(2*(standev**2)))

#-----Get_Bursts-------------------------------------------------------------------------------------------------------

def get_bursts(data,q_lo=50,q_hi=90,just_peaks=False,smooth=False,savgol=5,alg='cubic spline',times=None):

   '''Get Bursts

   Description:

    Takes a lightcurve and identifies 'bursts' in the data; short, discrete regions of increased flux.
    Returns the locations of all peaks identified as a list of tuples, each of which consist of two
    integers which correspond to the indices of the start and end of a peak in the original data.

   Inputs:

    data       - ARRAY: The data in which bursts are sought.
    q_low      - FLOAT: [Optional: Default=50] The percentile value of the data which will be used as
                        the low-pass threshold.  This threshold determines the edges of already-located
                        bursts, and thus changing it will change the quality of bursts but not the
                        quantity.
    q_hi       - FLOAT: [Optional: Default=90] The percentile value of the data which will be used as
                        the hi-pass threshold.  This threshold determines how much peak flux a previously
                        defined burst-candidate region must be before it is considered a true burst.
                        larger than q_med.
    just_peaks -  BOOL: [Optional: Default=False] If set to true, the function will return a list of peak
                        indices instead of a list of peak datasets.
    smooth     -  BOOL: [Optional: Default=False] If set to true, applies a univariate spline to the data
                        to smooth it.
    savgol     -   INT: [Optional: Default=5] The window size for the Savitsky-Golay filter
    alg        -STRING: [Optional: Default='cubic spline'] The algorithm to use to obtain peaks.
    times      - ARRAY: [Optional: Default=None] If 'load' selected for algorithm, user must also give the
                        time array associated with the data.


   Outputs:

    burst_locs -  LIST: A list of tuples containing the start and end indices of each burst.

   -J.M.Court, 2015'''


   if smooth:                                                             # If user has requested smoothing...
      savgol=int(savgol)
      if savgol%2==0:
         savgol+=1
      data=sgnl.savgol_filter(data,savgol,3)                              #  Smooth it!

   high_thresh=np.percentile(data,q_hi)
   low_thresh=np.percentile(data,q_lo)
   over_thresh=data>low_thresh                                            # Create a Boolean array by testing whether the input array is above the mid
                                                                          #  threshold.  Each region of consecutive 'True' objects is considered a burst-
                                                                          #  -candidate region.
   peak_locs=[]
   burst_locs=[]
   if alg=='cubic spline':
      while True:                                                            
                                                                          
         masked=np.array(data)*over_thresh                                # Reduce all data outside of burst-candidate regions to zero
         if max(masked)<high_thresh:                                      # If highest peak in all remaining burst-candidate regions is below the high threshold,
                                                                          #  assume there are no more bursts to be found.
            break

         peak_loc=masked.tolist().index(max(masked))                      # Find peak in remaining data
         peak_locs.append(peak_loc)                                       # Construct list of peak location
         i=peak_loc
         while i<len(data) and over_thresh[i]:                            # Scrub the True objects in the Boolean array corresponding to that peak's candidate
                                                                          #  region, thus removing it
            over_thresh[i]=False
            i+=1
         i=peak_loc-1
         while i>=0 and over_thresh[i]:
            over_thresh[i]=False
            i-=1

   elif alg=='load':
      if times==None:
         raise Exception('Must provide data times when using "load" mode!')
      times=np.array(times)
      loadfilename=raw_input('Burst File Name: ')
      loadfile=open(loadfilename,'r')
      with open(loadfilename,'r') as f:
         first_line=f.readline()
      if len(first_line.split(','))>6:
         dataex=np.array(first_line.split(','))    
         for l in dataex:
            peak_locs.append(float(l))  
      else:
         for line in loadfile:
            l=line.split(',')
            peak_locs.append(float(l[0]))
      loadfile.close()
      peak_locs=(np.array(peak_locs)-times[0])/(times[1]-times[0])
      peak_locs=peak_locs.astype(int)
         

   if just_peaks: return peak_locs
   peak_locs.sort()                                                       # Sort the list so peaks can be returned in chronological order

   start_col=get_dip(data,0,peak_locs[0])

   for i in range(len(peak_locs)):
      if i==len(peak_locs)-1:
         n_peak=len(data)-1
      else:
         n_peak=peak_locs[i+1]
      end_col=get_dip(data,peak_locs[i],n_peak)
      if end_col>start_col:                                               # Force any null-length 'bursts' to be removed
         burst_locs.append((start_col,end_col))
         start_col=end_col

   return burst_locs


#-----Get_Bursts_Windowed----------------------------------------------------------------------------------------------

def get_bursts_windowed(data,windows,q_lo=50,q_hi=90,smooth=False):

   '''Get Bursts

   Description:

    Takes a lightcurve and identifies 'bursts' in the data; short, discrete regions of increased flux.
    Returns the locations of all peaks identified as a list of tuples, each of which consist of two
    integers which correspond to the indices of the start and end of a peak in the original data.
    Currently can only return the location of all burst peaks.

   Inputs:

    data       - ARRAY: The data in which bursts are sought.
    windows    -   INT: The number of windows into which to divide the data
    q_low      - FLOAT: [Optional: Default=50] The percentile value of the data which will be used as
                        the low-pass threshold.  This threshold determines the edges of already-located
                        bursts, and thus changing it will change the quality of bursts but not the
                        quantity.
    q_hi       - FLOAT: [Optional: Default=90] The percentile value of the data which will be used as
                        the hi-pass threshold.  This threshold determines how much peak flux a previously
                        defined burst-candidate region must be before it is considered a true burst.
                        larger than q_med.
    smooth     -  BOOL: [Optional: Default=False] If set to true, applies a univariate spline to the data
                        to smooth it.
    savgol     -   INT: [Optional: Default=1] The window size for the Savitsky-Golay filter

   Outputs:

    burst_locs -  LIST: A list of tuples containing the peak of each burst.

   -J.M.C.Court, 2016'''

   windows=int(windows)
   windowlength=len(data)/windows
   win_starts=[windowlength*i for i in range(windows)]
   win_ends=[windowlength*(i+1) for i in range(windows-1)]
   win_ends.append(len(data)+1)
   peak_locs=[]

   for i in range(windows):

      new_bursts=np.array(get_bursts(data[win_starts[i]:win_ends[i]],q_lo=q_lo,q_hi=q_hi,just_peaks=True,smooth=smooth,savgol=(win_ends[i]-win_starts[i])/4.0))
      new_bursts+=win_starts[i]
      peak_locs+=new_bursts.tolist()

   return peak_locs


#-----Get Dip----------------------------------------------------------------------------------------------------------

def get_dip(data,start,finish,smooth=False,savgol=1):

   '''Get Dips

   Description:

    Returns the index of the lowest value between two given points in a dataset

   Inputs:

    data    - LIST:  The dataset in which a trough is to be found
    start   -  INT:  The index of the startpoint of the user-defined sub-range
    finish  -  INT:  The index of the endpoint of the user-defined sub-range
    smooth  -  BOOL: [Optional: Default=False] If set to true, applies a univariate spline to the data
                     to smooth it.

   Outputs:

    key_col -  INT: The index of the lowest value in the user-defined sub-range

   -J.M.Court, 2015'''

   data=np.array(data)
   data_l=np.arange(len(data))
   if smooth:
      savgol=int(savgol)
      if savgol%2==0:
         savgol+=1
      data=sgnl.savgol_filter(data,savgol,3)
   data=data*(data_l>=start)*(data_l<finish)
   data[data==0]=max(data)      
   keycol_loc=data.tolist().index(min(data))
   return keycol_loc 


#-----Get Phase--------------------------------------------------------------------------------------------------------

@mjit()
def get_phases(data,windows=1,q_lo=20,q_hi=90):

   '''Get Phases

   Description:

    Given a set of data points evenly spaced in time, returns the phase coordinate of each point.

   Inputs:

    data   -  ARRAY: The data in which bursts are sought.
    windows -   INT: [Optional: Default=1] The number of windows into which to divide the data for analysis
    q_low   - FLOAT: [Optional: Default=20] The percentile value of the data which will be used as
                     the low-pass threshold.  This threshold determines the edges of already-located
                     bursts, and thus changing it will change the quality of bursts but not the
                     quantity.
    q_hi    - FLOAT: [Optional: Default=90] The percentile value of the data which will be used as
                     the hi-pass threshold.  This threshold determines how much peak flux a previously
                     defined burst-candidate region must be before it is considered a true burst.
                     larger than q_med.

   Outputs:

    phases - ARRAY: The phase coordinates of the input data

   -J.M.C.Court, 2016'''

   data_keys=range(len(data))                                             # Generate the data indices
   data=np.array(data)                                                    # Format the input data

   data_phas=np.zeros(len(data),dtype=float)                              # Create the phase array
   peak_keys=get_bursts_windowed(data,windows,q_lo=50,q_hi=90,smooth=False)
   peak_keys.sort()                                                       # Fetch and sort the indices corresponding to peaks
   dip_keys=[]                                                            # Create empty array for indices of dips
   if peak_keys[0]!=0:
      dip_keys.append(get_dip(data,0,peak_keys[0]))                       # If there is no peak in the first slot, prepend a false dip at element 0
   for i in range(len(peak_keys)-1):
      dip_keys.append(get_dip(data,peak_keys[i],peak_keys[i+1]))          # For every pair of peaks, fetch the index of the dip between them
   if peak_keys[-1]!=data_keys[-1]:
      dip_keys.append(get_dip(data,peak_keys[-1],data_keys[-1]))          # If there is no peak in the last slot, append a false dip at element -1
   data_keys=np.array(data_keys)
   peak_keys=np.array(peak_keys)
   dip_keys=np.array(dip_keys)

   if peak_keys[0]==0:                                                    # Prepend false dips/peaks at 0 or at a negative index to force all comparisons
                                                                          # to be valid.  This has the effect of making the fit messy near the beginning and
                                                                          # end of the time series, but this is not aviodable.
      dip_keys=np.hstack((-1,dip_keys))                                   #
   elif dip_keys[0]==0:                                                   #
      peak_keys=np.hstack((-1,peak_keys))                                 #                                       ^
   elif peak_keys[0]<dip_keys[0]:                                         #                                       |
      dip_keys=np.hstack((0,dip_keys))                                    #                                       |
      peak_keys=np.hstack((-1,peak_keys))                                 #                                       |
   else:                                                                  #                                       |
      peak_keys=np.hstack((0,peak_keys))                                  #                                       |
      dip_keys=np.hstack((-1,dip_keys))                                   #                                       |
                                                                          #                                       |
   if peak_keys[-1]==data_keys[-1]:                                       # Appending false dips/peaks, see above |
      dip_keys=np.hstack((dip_keys,data_keys[-1]+1))
   elif dip_keys[-1]==data_keys[-1]:
      peak_keys=np.hstack((peak_keys,data_keys[-1]+1))
   elif peak_keys[-1]>dip_keys[-1]:
      dip_keys =np.hstack((dip_keys ,data_keys[-1]))
      peak_keys=np.hstack((peak_keys,data_keys[-1]+1)) 
   else:
      peak_keys=np.hstack((peak_keys,data_keys[-1]))
      dip_keys =np.hstack((dip_keys ,data_keys[-1]+1)) 
   
   for key in data_keys:                                                  # For each key, extract phase
      if key in peak_keys:                                                # If key is a peak, set phase to 0.5
         data_phas[key]=0.5
         continue
      elif key in dip_keys:                                               # If key is a dip, set phase to 0.0
         continue
      prev_peak=peak_keys[peak_keys<key][-1]                              # Otherwise, collect the keys of the previous/next peaks and dips
      prev_dip =dip_keys[ dip_keys<key][-1]
      next_peak=peak_keys[peak_keys>key][0]
      next_dip =dip_keys[ dip_keys>key][0]
      is_in_fall=(prev_peak>prev_dip)                                     # If the considered element's most recent peak was more recent than
                                                                          #  the most recent dip, then it is in the 'fall' regime (phase>0.5)
      prev_marker=float(max(prev_peak,prev_dip))                          # Collect the indices of the previous/next significant points of either type
      next_marker=float(min(next_peak,next_dip))
      phase=0.5*(float(key)-prev_marker)/(next_marker-prev_marker)        # Define phase as half the considered point's fractional distance between them
      if is_in_fall:
         phase+=0.5                                                       # Add 0.5 in the fall regime
      data_phas[key]=phase

   return data_phas  

#-----Get Phase INTP---------------------------------------------------------------------------------------------------

def get_phases_intp(data,windows=1,q_lo=20,q_hi=90,peaks=None,givespline=False):

   if peaks==None:
      peak_keys=get_bursts(data,q_lo=q_lo,q_hi=q_hi,just_peaks=True,smooth=False,savgol=5)
   else:
      peak_keys=peaks
   peak_keys.sort()
   peak_keys=np.array(peak_keys)

   data=np.array(data)
   p_phases=np.arange(0,len(peak_keys),1.0)

   spline=intp.PchipInterpolator(peak_keys, p_phases, extrapolate=True)
   spline.firstpeak=peak_keys[0]                                          # Store the keys of the first and last peaks in the spline object
   spline.lastpeak=peak_keys[-1]
   if givespline:
       return spline
   phases=spline(range(len(data)))
   phases=np.remainder(phases,1)
   return phases

#-----GTIMask----------------------------------------------------------------------------------------------------------

@mjit()
def gtimask(times,gtis):

   '''GTI Mask

   Description:

    Takes a time series and a GTI file taken from FITS and creates a mask to place over the time series
    with 'True' values over timestamps within the GTI and 'False' values elsewhere.

   Inputs:

    times -       LIST: A list of times, the x-axis of the data to be considered.  Must be in same
                        physical units as GTI.
    gtis  - FITS TABLE: A FITS table object containing a list of 2-element lists, each of which is
                        the start and endpoint respectively of a good time index.

   Outputs:

    mask  -       LIST: A mask to put over data to hide values outside of the GTI.

   -J.M.Court, 2015'''

   times=np.array(times)
   mask=np.zeros(len(times),dtype=bool)                                   # Set up initial blank list of 'False'

   for gti in gtis:                                                       # For every GTI index:
      smask=(times>gti[0]) & (times<gti[1])                               # Create a submask which is the 'and'ed product of times>gti_start and times<gti_end
      mask=mask|smask                                                     # 'or' the submask with the main mask
   return mask
         

#-----LBinify----------------------------------------------------------------------------------------------------------

@mjit()
def lbinify(x,y,ye,logres):

   '''Logarithmic Binify

   Decription:

    Takes a 2-dimenstional set of data which has already been evenly binned on the x-axis, and re-bins
    it into larger bins on the x-axis which are evenly spaced in log-space.  The original linear
    binning is retained in regions of the data where the logarithmic binning would be smaller than
    the data resolution.

   Inputs:

    x      -  LIST: The x-values of the two-dimensional data.
    y      -  LIST: The y-values of the two-dimensional data, must be the same length as x.
    ye     -  LIST: The errors associated with the y-values of the two-dimensional data, must be the
                     same length as x and y.
    logres - FLOAT: the size of the new x-axis bins, in log10-space, in which to re-bin the data.

   Outputs:

    xb     -  LIST: The re-binned x-values of the two-dimensional data, i.e. an array of the left-hand
                     edges of the new bins.
    yb     -  LIST: The re-binned y-values of the two-dimensional data.
    yeb    -  LIST: The errors associated with the rebinned y-values of the two-dimensional data.

   Outputs:

   -J.M.Court, 2015'''

   hinge=((x[1]-x[0])*10**logres)/((10**logres)-1)                        # Find the 'hinge' point at which to switch between linear and logarithmic binning

   lbin=np.log10(x[0])
   xb =10**(np.arange(lbin,np.log10(x[-1]),logres))                       # Setting up arrays to append binned values into
   yb =np.zeros(len(xb))
   yeb=np.zeros(len(xb))

   hingel=sum((xb)<=hinge)                                                # Getting the ID of the hinge-point in the log

   xbl=len(xb)

   for i in range(hingel,xbl):

      lowid=int(((10**((i*logres)+lbin))-x[0])/(x[1]-x[0]))               # Calculate the ID of the lowest linear bin that corresponds to this log bin
      uppid=int(((10**(((i+1)*logres)+lbin))-x[0])/(x[1]-x[0]))           # Calculate the ID of the highest linear bin that corresponds to this log bin
      if uppid>lowid:
         yb[i]=np.mean(y[lowid:uppid])
         yeb[i]=(np.sqrt(sum(np.array(ye[lowid:uppid])**2)))/int(uppid-lowid)
      else:
         yb[i]=0                                                          # If no data found, error=power=0
         yeb[i]=0

   mask=x<hinge
   lmask=xb>hinge

   xf=np.append(x[mask],xb[lmask])
   yf=np.append(y[mask],yb[lmask])
   yef=np.append(ye[mask],yeb[lmask])

   return xf,yf,yef


#-----LeahyN-----------------------------------------------------------------------------------------------------------

@mjit()
def leahyn(data,counts,datres):

   '''Leahy Normaliser

   Description:

    Takes a raw FFT-algorithm output spectrum and normalises it by the Leahy convention.

   Inputs:

    data   - LIST: the spectrum to be normalised.
    counts -  INT: the number of photon events in the data sample for which the spectrum was created.
    datres -  INT: the number of time bins in the data sample for which the spectrum was created.

   Outputs:

    leahy  - LIST: the Leahy-normalised spectrum.

   -J.M.Court, 2015'''

   leahy=2*(abs(data[0:datres/2.0])**2)/counts
   return leahy


#-----Lh2RMS-----------------------------------------------------------------------------------------------------------

@mjit()
def lh2rms(leahy,rate,bg,const):

   '''Leahy 2 RMS Converter

   Description:

    Takes a Leahy-normalised spectrum and re-normalises it by the (RMS/Mean)^2 convention.

   Inputs:

    leahy -  LIST: the leahy-normalised spectrum to be re-normalised.
    rate  - FLOAT: the source + background count rate (per second) of the data sample from which the
                    spectrum was created.
    bg    - FLOAT: the background count rate (per second) of the data sample from which the spectrum
                    was created.
    const - FLOAT: the average power given in Leahy normalisation for pure white noise.  Theoretically
                    const=2, but in practice is slightly lower and varies between telescopes.

   Outputs:

    rms   -  LIST: the (RMS/Mean)^2-normalised spectrum.

   -J.M.Court, 2015'''

   denom=(rate-bg)**2
   if denom==0:
      mult=0.0
   else:
      mult=1.0/denom
   rms= (leahy-const)*(rate)*mult
   return rms


#-----LhConst----------------------------------------------------------------------------------------------------------

def lhconst(data):

   '''LHS Const

   Decription:

    Finds the normalisation of white noise in a given power Leahy-normalised power spectrum.  Assumes
    no spectral features in the last 20% of the spectrum.

   Inputs:

    data  -  LIST: the data to be converted.

   Outputs:

    const - FLOAT: the normalisation of white noise.

   -J.M.Court, 2015'''

   def leahynoise(x,a):
      return a

   olen=len(data)
   olen=int((4/5.0)*olen)
   datav=data[olen:]
   const=optm.curve_fit(leahynoise,range(len(datav)),datav)
   const=const[0][0]
   if const>2.5 or const<1.5:
      print "WARNING: Leahy constant of "+str(const)+" outside of accepted range!"
   return const


#-----Lomb_Scargle-----------------------------------------------------------------------------------------------------

@mjit()
def lomb_scargle(x,y,ye,freqs):

   '''Lomb Scargle

   Decription:

    Returns the Lomb-Scargle periodogram of data provided by the user, scanning over a set of frequencies
    also provided by the user.

   Inputs:

    x     -  LIST: the time-array for the data to be converted.
    y     -  LIST: the data associated with time array x.
    ye    -  LIST: the errors associated with each point in y.
    freqs -  LIST: the list of discrete frequencies over which the user wants the L-S periodogram to scan.

   Outputs:

    pgram - ARRAY: the Lomb-Scargle power associated with each frequency provided in freqs.

   -J.M.Court, 2015'''

   assert len(x)==len(y)
   x=np.array(x)
   y=np.array(y)
   ye=np.array(ye)
   x=x[y>0]                                                               # Weed out any negative values here before they break things
   ye=ye[y>0]
   y=y[y>0]
   w=safe_div(np.ones(len(ye)),ye**2)
   y=y-(sum(y*w)/sum(w))

   freqs=np.array(freqs)*2*pi

   wt=np.multiply.outer(x,freqs)

   sin2wt=np.sin(2*wt)
   cos2wt=cosfromsin(2*wt,sin2wt)

   tau=(np.arctan(np.sum(sin2wt,axis=0)/np.sum(cos2wt,axis=0)))/(2*freqs)

   wttau=wt-(freqs*tau)

   sinw=np.sin(wttau)
   cosw=cosfromsin(wttau,sinw)

   yT=np.vstack(y)

   ysin=yT*sinw
   ycos=yT*cosw

   norm=2*np.sum(w*(y**2)/(len(y)-2))

   w=np.vstack(w)
   pgram=((np.sum(w*ycos,axis=0)**2)/np.sum((w*cosw)**2,axis=0)+(np.sum(w*ysin,axis=0)**2)/np.sum((w*sinw)**2,axis=0)) / norm

   return pgram


#-----MXRebin----------------------------------------------------------------------------------------------------------

@mjit()
def mxrebin(spcdata,spcerrs,xaxis,good,bfac):

   '''Matrix X-Rebin

   Description:

    Takes 2-Dimensional arrays of data and errors which are linearly binned on the x-axis and rebins
    them by a factor of the user's choosing, returning new data array, error array, x-axis and Boolean
    'good' array.

   Inputs:

    spcdata   - ARRAY: some 2 dimensional array of values.
    spcerrs   - ARRAY: a 2 dimensional array containing the errors associated with the values in spcdata.
                       Must have the same dimensions as spcdata.
    xaxis     - ARRAY: an array of values corresponding to the x-values of columns in spcdata.  Must
                       be the same length as a row of spcerrs and must be linearly spaced.
    good      - ARRAY: A Boolean array with as many elements as spcdata has columns.  Its entries are
                       True unless the corresponding row in spcdata does not correspond to a valid time
                       within the GTIs of the photon count data.
    bfac      -   INT: the binning factor, i.e. the ratio between new bin size and old bin size.

   Outputs:

    b_spcdata - ARRAY: the rebinned 2-dimensional data array
    b_spcerrs - ARRAY: the rebinned 2-dimensional error array
    b_xaxis   - ARRAY: the rebinned x-axis
    b_good    - ARRAY: the rebinned 'good' array

   -J.M.Court, 2015'''

   spx=len(spcdata[:,0])                                                  # x-dimension of new matrix matches old matrix
   spy=int(len(spcdata[0,:])/bfac)                                        # y-dimension of new matrix equals the y dimension of the old matrix divided by the binning factor
   b_good=np.zeros(spy,dtype=bool)                                        # array to label good columns
   b_spcdata=np.zeros([spx,spy])                                          # Create the new matrix
   b_spcerrs=np.zeros([spx,spy])                                          # Create the new error matrix

   for i in range(spx):                                                   # For each freq row of fourgr:
      for j in range(spy):                                                # For each time col of fourgr:

         celltot=0                                                        # Create a running total to deposit into the matrix cell
         errtot=0                                                         # Create a running total to deposit into the error matrix cell

         for k in range(bfac):                                            # Calculate extent of current bin

            b_good[j]=(True)                                              # Assume column is good

            if not good[j*bfac+k]:                                        # If any of the columns being summed were flagged as bad, flag new column as bad

               celltot=0
               b_good[j]=(False)
               break                                                      # Don't attempt to find spectrum in column if any constituent columns are bad

            celltot+=spcdata[i,j*bfac+k]                                  # Sum all values that fall within the given bin
            errtot+=((spcerrs[i,j*bfac+k])**2)

         b_spcdata[i,j]=celltot/bfac                                         # Divide the cell total by the bin multiplier to convert to a mean
         b_spcerrs[i,j]=np.sqrt(errtot)/bfac

   b_xaxis=[]                                                             # Calculate new x-axis

   for i in range(int(len(xaxis)/bfac)):
      b_xaxis.append(xaxis[i*bfac])

   return b_spcdata,b_spcerrs,b_xaxis,b_good


#-----Nones------------------------------------------------------------------------------------------------------------

def nones(shape):
    
    '''Function to basically do what np.zeros and np.ones do, but creating an array of Nones.

    -J.Coxon, 2015'''

    makeNone = np.vectorize(lambda x: None)
    return makeNone(np.zeros(shape))


#-----PDColEx----------------------------------------------------------------------------------------------------------

@mjit()
def pdcolex2(y1,y2,ye1,ye2,gmask):

   '''Plot Demon Colour Extract (2D)

   Description:

    Takes two equally spaced flux series in the same set of time co-ordinates and constructs an array
    representing the file2 / file1 colour over the same time period.  Also returns the sum of the flux
    of the two series.

   Inputs:

    y1,y2   - ARRAYS: The first and second flux series respectively
    ye1,ye2 - ARRAYS: The errors of the first and second flux series respectively
    gmask   -  ARRAY: A mask of Boolean values with False at every point outside of a GTI in this
                      time series

   Outputs:

    flux    -  ARRAY: The sum total flux from the two bands
    fluxe   -  ARRAY: The errors of 'flux'
    col21   -  ARRAY: The file2/file1 colour
    col21e  -  ARRAY: The error on "col21"

   -J.M.Court, 2015'''

   warnings.filterwarnings("ignore")                                      # Div 0 errors are a real possibility.  This is me ignoring them...

   y=   {}                                                                # Prepare libraries for data
   ye=  {}
   col= {}
   cole={}

   y[1]=y1[gmask]                                                         # Mask data, store in library
   y[2]=y2[gmask]
   ye[1]=ye1[gmask]
   ye[2]=ye2[gmask]

   flux=y[1]+y[2]                                                         # Get total flux
   fluxe=np.sqrt(ye[1]**2+ye[2]**2)                                       # Get flux error

   for i in range(1,3):                                                   # For the ith possible numerator band
      for j in range(1,3):                                                # For the jth possible denominator band
         if j!=i:                                                         # Prevents taking x/x colour
            ld=int(str(i)+str(j))
            col[ld]=(y[i]/y[j])                                           # Fetch colour
            cole[ld]=col[ld]*np.sqrt(((ye[i]/y[i])**2)+((ye[j]/y[j])**2)) # Fetch colour error

   return flux,fluxe,y,ye,col,cole

@mjit()
def pdcolex3(y1,y2,y3,ye1,ye2,ye3,gmask):

   '''Plot Demon Colour Extract (3D)

   See help for Plot Demon Colour Extract (2D)

   -J.M.Court, 2015'''

   warnings.filterwarnings("ignore")                                      # Div 0 errors are a real possibility.  This is me ignoring them...

   y=   {}                                                                # Prepare libraries for data
   ye=  {}
   col= {}
   cole={}

   y[1]=y1[gmask]                                                         # Mask data, store in library
   y[2]=y2[gmask]
   y[3]=y3[gmask]
   ye[1]=ye1[gmask]
   ye[2]=ye2[gmask]
   ye[3]=ye3[gmask]

   flux=y[1]+y[2]+y[3]                                                    # Get total flux
   fluxe=np.sqrt(ye[1]**2+ye[2]**2+ye[3]**2)                              # Get flux error

   for i in range(1,4):                                                   # For the ith possible numerator band
      for j in range(1,4):                                                # For the jth possible denominator band
         if j!=i:                                                         # Prevents taking x/x colour
            ld=int(str(i)+str(j))
            col[ld]=(y[i]/y[j])                                           # Fetch colour
            cole[ld]=col[ld]*np.sqrt(((ye[i]/y[i])**2)+((ye[j]/y[j])**2)) # Fetch colour error

   return flux,fluxe,y,ye,col,cole

#-----PDload-----------------------------------------------------------------------------------------------------------

def pdload(filename,isplotd):
    
   '''PlotDemon loader

   Description:

    Calls plotdld or csvload depending on the type of file to be opened  .
    
   Inputs:
   
    filename - STRING: The absolute or relative path to the location of the file that will be opened.
    isplotd  - BOOL  : Whether the file to be opened has been identified as a .plotd file
    
   Outputs:

    times    -      ARRAY: An array, the elements of which are the left-hand edges of the time bins
                           into which counts have been binned.  Units of seconds.
    counts   -      ARRAY: The number of counts in each bin defined in 'times'.
    errors   -      ARRAY: The 1-sigma errors associated with each value in 'counts'
    binsize  -      FLOAT: The size of each bin in 'times'.  Saved for speed upon loading.
    gti      - FITS TABLE: The table of GTI values from the event data .fits file.
    mxpcus   -        INT: The maximum number of PCUs active at any one time during the observation.
    bgpcu    -      FLOAT: An estimate of the count rate of the background flux during the full 
                           observation, in counts per second per PCU, multiplied by the number of
                           PCUs.
    flavour  -     STRING: A useful bit of text to put on plots to help identify them later on.
    chanstr  -     STRING: A string containing the high and low channel numbers separated by a dash.
    mission  -     STRING: The name of the satellite
    obsdata  -      TUPLE: The first element is the name of the object, the second is the observation
                           ID.
    version  -     STRING: The Version of FITSGenie in which the file was created
    
   -J.M.Court, 2015'''
    
   if isplotd:
       return plotdld(filename)
   else:
       return csvload(filename)


#-----PlotdLd----------------------------------------------------------------------------------------------------------

def plotdld(filename):

   '''.Plotd Load

   Description:

    Opens a .plotd file and retrieves the useful data from it.

   Inputs:

    filename - STRING: The absolute or relative path to the location of the file that will be opened.

   Outputs:

    times    -      ARRAY: An array, the elements of which are the left-hand edges of the time bins
                           into which counts have been binned.  Units of seconds.
    counts   -      ARRAY: The number of counts in each bin defined in 'times'.
    errors   -      ARRAY: The 1-sigma errors associated with each value in 'counts'
    binsize  -      FLOAT: The size of each bin in 'times'.  Saved for speed upon loading.
    gti      - FITS TABLE: The table of GTI values from the event data .fits file.
    mxpcus   -        INT: The maximum number of PCUs active at any one time during the observation.
    bgpcu    -      FLOAT: An estimate of the count rate of the background flux during the full 
                           observation, in counts per second per PCU, multiplied by the number of
                           PCUs.
    flavour  -     STRING: A useful bit of text to put on plots to help identify them later on.
    chanstr  -     STRING: A string containing the high and low channel numbers separated by a dash.
    mission  -     STRING: The name of the satellite
    obsdata  -      TUPLE: The first element is the name of the object, the second is the observation
                           ID.
    version  -     STRING: The Version of FITSGenie in which the file was created

   -J.M.Court, 2015'''

   try:
      readfile=open(filename,'rb')                                        # Open the .plotd file
   except:
      print ''
      print 'File '+filename+' not found!  Aborting!'
      signoff()
      exit()
   data=pickle.load(readfile)                                            # Unpickle the .plotd file

   times=np.array(data['time'])                                          # Unleash the beast! [extract the file]
   rates=np.array(data['rate'])
   errors=np.array(data['errs'])
   tstart=data['tstr']
   binsize=data['bsiz']
   gti=data['gtis']
   mxpcus=data['pcus']
   bgest=data['bkgr']
   bgsub=data['bsub']
   bgdata=data['bdat']
   flavour=data['flav']
   chanstr=data['chan']
   mission=data['miss']
   obsdata=data['obsd']
   version=data['vers']

   bgpcu=bgest*mxpcus                                                     # Collect background * PCUs

   readfile.close()

   return times,rates,errors,tstart,binsize,gti,mxpcus,bgpcu,bgsub,bgdata,flavour,chanstr,mission,obsdata,version


#-----csvLoad----------------------------------------------------------------------------------------------------------

def linesplit(lx,iscomma):
   if iscomma:
      return lx.split(',')
   return lx.split()

def csvload(filename):

   '''.csv Loader
   
   Description:
   
    Extracts information from a csv file.  Column delimiters automatically set using python string object's
    .split() function.  Assumes the file is of format 'time,rate' if two columns, 'time,rate,rate_error' if
    three columns of 'time,time_error,rate,rate_error; if four or more columns.
   
   Inputs:

    filename - STRING: The absolute or relative path to the location of the file that will be opened.

   Outputs:

    times    -      ARRAY: An array, the elements of which are the left-hand edges of the time bins
                           into which counts have been binned.  Units of seconds.
    counts   -      ARRAY: The number of counts in each bin defined in 'times'.
    errors   -      ARRAY: The 1-sigma errors associated with each value in 'counts'
    binsize  -      FLOAT: The size of each bin in 'times'.  Saved for speed upon loading.
    gti      - FITS TABLE: The table of GTI values from the event data .fits file.
    mxpcus   -        INT: The maximum number of PCUs active at any one time during the observation.
    bgpcu    -      FLOAT: An estimate of the count rate of the background flux during the full 
                           observation, in counts per second per PCU, multiplied by the number of
                           PCUs.
    flavour  -     STRING: A useful bit of text to put on plots to help identify them later on.
    chanstr  -     STRING: A string containing the high and low channel numbers separated by a dash.
    mission  -     STRING: The name of the satellite
    obsdata  -      TUPLE: The first element is the name of the object, the second is the observation
                           ID.
    version  -     STRING: The Version of FITSGenie in which the file was created
    
   -J.M.Court, 2015'''
   
   f=open(filename,'r')

   times=[]
   rates=[]
   errors=[]
   got_firstline=False

   haserrs=False

   xind=0

   for line in f:
       if not got_firstline:
          if ((len(line.split())<2) and (len(line.split(','))<2)):        # Assume lines with <2 columns, or a non-number in col
             continue
          if ',' in line:
             iscomma=True
          else:
             iscomma=False
          l=linesplit(line,iscomma)
          try:                                                            #  zero are part of the header, skip 'em
             float(l[0])
          except:
             continue
          llen=len(l)
          if llen<3:                                                      # Interpret csv data depending on how many columns
             yind=1
             haserrs=False
          elif llen<4:
             yind=1
             haserrs=True
             eind=2
          else:
             yind=2
             haserrs=True
             eind=3
          got_firstline=True
       l=linesplit(line,iscomma)
       times.append(float(l[xind]))
       rates.append(float(l[yind]))
       if haserrs:
          errors.append(float(l[eind]))

   if not haserrs:
       errors=np.sqrt(np.array(rates))                                    # Assume root flux errors if none in file     

   times=np.array(times)
   rates=np.array(rates)
   errors=np.array(errors)

   tstart=times[0]
   binsize=times[1]-times[0]
   gti=None
   mxpcus=1
   bgpcu=1
   bgsub=0
   bgdata=None
   flavour='csv'
   chanstr='Unknown'
   mission='Unknown'
   obsdata=['Unknown','Unknown']
   version='From csv'
   return times,rates,errors,tstart,binsize,gti,mxpcus,bgpcu,bgsub,bgdata,flavour,chanstr,mission,obsdata,version

   
#-----PlotdSv----------------------------------------------------------------------------------------------------------

@mjit()
def plotdsv(filename,times,rates,errors,tstart,binsize,gti,mxpcus,bgest,bgsub,bgdata,flavour,chanstr,mission,obsdata,version):

   '''.Plotd Save

   Description:

    Takes the input of the data products required to create a .plotd file (to read with plotdemon)
    and creates a .plotd file at a location given as the first input.

   Inputs:

    filename -     STRING: The absolute or relative path to the location of the file that will be 
                           created.
    times    -      ARRAY: An array, the elements of which are the left-hand edges of the time bins
                           into which counts have been binned.  Units of seconds.
    counts   -      ARRAY: The number of counts in each bin defined in 'times'.
    binsize  -      FLOAT: The size of each bin in 'times'.  Saved for speed upon loading.
    gti      - FITS TABLE: The table of GTI values from the event data .fits file.
    mxpcus   -        INT: The maximum number of PCUs active at any one time during the observation.
    bgest    -      FLOAT: An estimate of the count rate of the background flux during the full 
                           observation, in counts per second per PCU.
    flavour  -     STRING: A useful bit of text to put on plots to help identify them later on.
    chanstr  -     STRING: A string containing the high and low channel numbers separated by a dash.
    mission  -     STRING: The name of the satellite
    obsdata  -      TUPLE: The first element is the name of the object, the second is the observation
                           ID.
    version  -     STRING: The Version of FITSGenie in which the file was created

   Outputs:

    [none]

   -J.M.Court, 2015'''

   savedata={}                                                            # Open library object to save in file

   savedata['time']=times                                                 # Dump each piece of data into an appropriate library element
   savedata['rate']=np.array(rates)
   savedata['errs']=np.array(errors)
   savedata['tstr']=tstart
   savedata['bsiz']=binsize
   savedata['gtis']=gti
   savedata['pcus']=mxpcus
   savedata['bkgr']=bgest
   savedata['flav']=flavour
   savedata['chan']=chanstr
   savedata['miss']=mission
   savedata['bsub']=bgsub
   savedata['bdat']=bgdata
   savedata['obsd']=obsdata
   savedata['vers']=version

   filename=uniqfname(filename,'plotd')                                   # Get the next available name of form filename(x).plotd
   wfile = open(filename, 'wb')                                           # Open file to write to

   pickle.dump(savedata,wfile)                                            # Pickle the data (convert into bitstream) and dump to file
   wfile.close()                                                          # Close file

   return filename


#-----RMS_N------------------------------------------------------------------------------------------------------------

@mjit()
def rms_n(data,counts,datres,rate,bg,const):

   '''RMS Normaliser

   Description:

    Takes a raw FFT-algorithm output spectrum and normalises it by the (RMS/Mean)^2 convention.

   Inputs:

    data   -  LIST: the spectrum to be normalised.
    counts -   INT: the number of photon events in the data sample for which the spectrum was created.
    datres -   INT: the number of time bins in the data sample for which the spectrum was created.
    rate   - FLOAT: the source + background count rate (per second) of the data sample from which the
                    spectrum was created.
    bg     - FLOAT: the background count rate (per second) of the data sample from which the spectrum
                    was created.
    const  - FLOAT: the average power given in Leahy normalisation for pure white noise.  Theoretically
                    const=2, but in practice is slightly lower and varies between telescopes.

   Outputs:

    rms    -  LIST: the Leahy-normalised spectrum.

   -J.M.Court, 2015'''

   leahy=leahyn(data,counts,datres)                                       # Leahy normalise the data
   rms=lh2rms(leahy,rate,bg,const)                                        # Convert to RMS
   return rms


#-----RMS--------------------------------------------------------------------------------------------------------------

@mjit()
def rms(data,data_err=[0]):

   '''RMS

   Description:

    Returns the fractional RMS of a 1-dimensional data set

   Inputs:

    data     - LIST: The data to find the RMS of.
    data_err - LIST: The errors on those data points

   Outputs:

    rms  - FLOAT: The rms of the data

   -J.M.Court, 2015'''

   if np.mean(data)==0:
      return 'div0'
   nr=1.0/len(data)
   mn=np.mean(data)
   vd=np.sum((data-mn)**2)-np.sum((data_err)**2)
   rms= ((vd*nr)**0.5)/abs(mn) 
   return rms


#-----Safe_Div---------------------------------------------------------------------------------------------------------

@mjit()
def safe_div(x,y):

   '''Safe Div

   Description:

    Divides the first inputs by the second inputs if the latter is nonzero.  If an element of the second input is zero,
    the corresponding element in the output is zero.

   Inputs:

    x - ARRAY: The numerator
    y - ARRAY: The denominator

   Outputs:

    r - ARRAY: The result of division, or zero

   -J.M.Court, 2015'''

   r=np.zeros(len(y))
   r[y!=0]=x[y!=0]/y[y!=0]
   return r


#-----SignOff----------------------------------------------------------------------------------------------------------

def signoff():

   '''Sign Off

   Description:

    Prints an underline with some spaces.  That's all.

   -J.M.Court, 2015'''

   print ''
   print '------------------------------------------------'
   print ''


#-----sinfromcos-------------------------------------------------------------------------------------------------------

@mjit()
def sinfromcos(x,cosx):

   '''Sine from Cosine

   Description:

    Returns the sine of an array of values when also given their cosines.  Using this function is
    faster than using sine if the cosine values are already stored.

   Inputs:

    x    - ARRAY: the array of values to calculate the sine of.
    cosx - ARRAY: the cosines array of values to calculate the sines of.

   Outputs:

    sinx - ARRAY: the sines of x.

   -J.M.Court, 2015'''

   sinx=np.absolute((1-cosx**2)**0.5)
   signx=np.sign(((x+pi)%(2*pi))-pi)
   return sinx*signx

@mjit()
def cosfromsin(x,sinx):

   '''Sine from Cosine

   Description:

    Returns the cosine of an array of values when also given their sines.  Using this function is
    faster than using cosine if the sine values are already stored.

   Inputs:

    x    - ARRAY: the array of values to calculate the cosine of.
    sinx - ARRAY: the sines array of values to calculate the cosines of.

   Outputs:

    cosx - ARRAY: the sines of x.

   -J.M.Court, 2015'''

   cosx=np.absolute((1-sinx**2)**0.5)
   signx=np.sign(((x-pi/2)%(2*pi))-pi)
   return cosx*signx


#-----SLPlot-----------------------------------------------------------------------------------------------------------

def slplot(x,y,ye,xlabel,ylabel,title,figid="",typ='both',errors=True):

   '''Standard/Log Plotter

   Description:

    Takes a two-dimensional array of data and plots it twice on the same figure; once on standard
    linear axes, and once on logarithmic axes.

   Inputs:

    x      -   LIST: The x-values of the two-dimensional data.
    y      -   LIST: The y-values of the two-dimensional data, must be the same length as x.
    xlabel - STRING: The title of the x-axis on the graphs.
    ylabel - STRING: The title of the y-axis on the graphs.
    title  - STRING: The title of the figure.
    figid  - STRING: Gives the Figure a name such that it can be easily manipulated or closed outside
                     of this function.
    typ    - STRING: User can input 'log' or 'lin' to only display plot of that type.
    errors -   BOOL: True or False: whether to display errorbars on plots

   Outputs:

    filename - STRING: The filename actually used when saving

   -J.M.Court, 2015'''

   pl.close(figid)                                                        # Close any previous plot of this type
   pl.figure(figid)                                                       # Create spectrum plot

   if typ in ('lin','both'):

      if typ=='lin':
         pl.subplot(111)                                                  # If 'lin' passed as typ word, only create one subplot
      else:
         pl.subplot(211)                                                  # If 'both' passed as typ word, make 1st of 2 subplots
      pl.grid(True,which="both")
      pl.xlabel(xlabel)
      pl.ylabel(ylabel)
      pl.title(title)
      if errors:
         pl.errorbar(x,y,ye)                                              # Plot data
      else:
         pl.plot(x,y)
      pl.plot([x[0],x[-1]],[0,0])

   if typ in ('log','both'):

      if typ=='log':
         ax=pl.subplot(111)                                               # If 'log' passed as typ word, only create one subplot
      else:
         ax=pl.subplot(212)                                               # If 'both' passed as typ word, make 2nd of 2 subplots
      pl.xlabel(xlabel)
      pl.ylabel(ylabel)
      pl.title(title)
      if errors:
         pl.errorbar(x,abs(y),ye,fmt='k')                                 # Plot data
      else:
         pl.plot(x,abs(y),'k')                                            # Plot log-log data
      ax.set_xscale('log')
      ax.set_yscale('log')
      pl.grid(True,which="both")

   #if typ in ('lin','log','both'):
   #   pl.show(block=False)                                               # Show both plots together
   else:
      print 'Invalid typ!  No plot shown.'                                # Complain if none of 'lin', 'log' or 'both are given as typ word


#-----Spliner----------------------------------------------------------------------------------------------------------

@mjit()
def spliner(data,errors=None):

   '''Spliner

   Description:

    Smooths evenly-sampled data using a first-order Univariate spline algorithm

   Inputs:

    data        - ARRAY: The data to be smoothed
    errors      - ARRAY: [Optional] The errors associated with the errors

   Outputs:

    smooth_data - ARRAY: The smoothed data values

   - J.M.C.Court,2016'''

   if type(errors)==type(None):
      errors=np.ones(len(data))                                           #  If no errors given, do not weight points
   else:
      assert len(errors)==len(data)
   spline=intp.UnivariatrSpline(range(len(data)),data,w=errors**-1,k=2)   #  Construct a univariate spline function from the data
   #pl.figure()                                                           #  Uncomment these 4 lines to allow for display of smoothed and pre-smoothed data
   #pl.plot(data_keys,data,'0.5')
   smooth_data=spline(range(len(data)))                                   #  Reconstruct the data using the spline
   #pl.plot(data_keys,data,'k')
   #pl.show(block=True)
   return smooth_data

#-----SpecaLd----------------------------------------------------------------------------------------------------------

def specald(filename):

   '''.Speca Load

   Description:

    Opens a .speca file and retrieves the useful data from it.

   Inputs:

    filename - STRING: The absolute or relative path to the location of the file that will be opened.

   Outputs:

    spcdata  -  ARRAY: An array, the elements of which are the non-normalised power spectra taken over
                       different equal-length time intervals of the photon count data being considered.
    good     -  ARRAY: A Boolean array with as many elements as spcdata has columns.  Its entries are
                       True unless the corresponding row in spcdata does not correspond to a valid time
                       within the GTIs of the photon count data.
    rates    -  ARRAY: The count rate per second of photons in the time interval represented by the
                       corresponding row of spcdata.
    prates   -  ARRAY: An array of floats with as many elements as spcdata has columns.  Denotes the
                       highest number of total counts in each timing window when the counts are binned
                       on a time of binsize * binfac
    trates   -  ARRAY: An array of floats with as many elements as spcdata has columns.  Denotes the
                       lowest number of total counts in each timing window when the counts are binned
                       on a time of binsize * binfac
    phcts    -    INT: The total number of photons detected, overall.
    bg       -  ARRAY: An estimate of the count rate of the background flux during the full observation,
                       in counts per second per PCU, multiplied by the number of PCUs active in the
                       corresponding row of specdata.
    binsize  -  FLOAT: The size in seconds of the time bins used when converting event data into time
                       binned photon count data.
    foures   -  FLOAT: The length of time corresponding to a single row of spcdata, in seconds.
    bgest    -  FLOAT: An estimate of the count rate of the background flux during the full observation,
                       in counts per second per PCU.
    flavour  - STRING: A useful bit of text to put on plots to help identify them later on.
    cs       - STRING: A string containing the high and low channel numbers separated by a dash.
    mission  - STRING: The name of the satellite
    obsdata  -  TUPLE: The first element is the name of the object, the second is the observation ID.
    wtype    - STRING: A string denoting the functional geometry of the windows used for Fourier analysis.
    slide    -  FLOAT: The time separation of the start times of windows represented by two consecutive
                       rows of spcdata, in seconds
    binfac   -    INT: See trates or prates
    version  - STRING: The Version of FITSGenie in which the file was created

   -J.M.Court, 2015'''

   try:
      readfile=open(filename,'rb')                                        # Open the .speca file
   except:
      print ''
      print 'File not found!  Aborting!'
      signoff()
      exit()
   data=pickle.load(readfile)                                             # Unpickle the .speca file

   spcdata=data['data']                                                   # Unleash the beast! [extract the file]
   good=np.array(data['good'])
   rates=np.array(data['rate'])
   prates=np.array(data['prts'])
   trates=np.array(data['trts'])
   phcts=data['phct']
   n_pcus=np.array(data['npcu'])
   binsize=data['bsiz']
   bgest=data['bkgr']
   foures=data['fres']
   flavour=data['flav']
   cs=data['chan']
   mission=data['miss']
   obsdata=data['obsd']
   wtype=data['wndw']
   slide=data['slid']
   binfac=data['sbin']
   version=data['vers']

   readfile.close()

   print ''

   print 'Power spectra taken over '+str(foures)+'s of data each'
   print str(sum(good))+'/'+str(len(spcdata))+' power spectra are good'
   if flavour!='':
      print "Flavour is '"+str(flavour)+"'"                               # Only print flavour if there is a flavour to print

   bg=n_pcus*bgest

   return spcdata,good,rates,prates,trates,phcts,bg,binsize,foures,bgest,flavour,cs,mission,obsdata,wtype,slide,binfac,version


#-----SpecaSv----------------------------------------------------------------------------------------------------------

@mjit()
def specasv(filename,spcdata,good,rates,prates,trates,phcts,npcus,binsize,bgest,foures,flavour,cs,mission,obsdata,wtype,slide,spcbinfac,version):

   '''.Speca Save

   "Should've gone to SpecaSaver..."

   Description:

    Takes the input of the data products required to create a .speca file (to read with specangel)
    and creates a .speca file at a location given as the first input.

   Inputs:

    filename  - STRING: The absolute or relative path to the location of the file that will be created.
    spcdata   -  ARRAY: An array, the elements of which are the non-normalised power spectra taken over
                        different equal-length time intervals of the photon count data being considered.
    good      -  ARRAY: A Boolean array with as many elements as spcdata has columns.  Its entries are
                        True unless the corresponding row in spcdata does not correspond to a valid time
                        within the GTIs of the photon count data.
    rates     -  ARRAY: An array of floats with as many elements as spcdata has columns.  Denotes the
                        total number of photon counts in the time interval represented by the corresponding
                        row in spcdata, divided by the width of a column in seconds.
    prates    -  ARRAY: An array of floats with as many elements as spcdata has columns.  Denotes the
                        highest number of total counts in each timing window when the counts are binned
                        on a time of binsize * spcbinfac
    trates    -  ARRAY: An array of floats with as many elements as spcdata has columns.  Denotes the
                        lowest number of total counts in each timing window when the counts are binned
                        on a time of binsize * spcbinfac
    phcts     -    INT: The total number of photons detected, overall.
    npcus     -  ARRAY: An array of ints with as many elements as spcdata has columns.  States the number
                        of detectors that were active when the data represented by the corresponding row
                        of spcdata was recorded.
    binsize   -  FLOAT: The size in seconds of the time bins used when converting event data into time
                        binned photon count data.
    bgest     -  FLOAT: An estimate of the count rate of the background flux during the full observation,
                        in counts per second per PCU.
    foures    -  FLOAT: The length of time corresponding to a single row of spcdata, in seconds.
    flavour   - STRING: A useful bit of text to put on plots to help identify them later on.
    cs        - STRING: A string containing the high and low channel numbers separated by a dash.
    mission   - STRING: The name of the satellite
    obsdata   -  TUPLE: The first element is the name of the object, the second is the observation ID.
    wtype     - STRING: A string denoting the functional geometry of the windows used for Fourier analysis.
    slide     -  FLOAT: The time separation of the start times of windows represented by two consecutive
                        rows of spcdata, in seconds
    spcbinfac -    INT: See trates or prates
    version   - STRING: The Version of FITSGenie in which the file was created

   Outputs:

    filename - STRING: The filename actually used when saving

   -J.M.Court, 2015'''

   savedata={}                                                            # Open library object to save in file

   savedata['data']=spcdata                                               # Dump each piece of data into an appropriate library element
   savedata['good']=good
   savedata['rate']=rates
   savedata['prts']=prates
   savedata['trts']=trates
   savedata['phct']=phcts
   savedata['npcu']=npcus
   savedata['bsiz']=binsize
   savedata['bkgr']=bgest
   savedata['fres']=foures
   savedata['flav']=flavour
   savedata['chan']=cs
   savedata['miss']=mission
   savedata['obsd']=obsdata
   savedata['wndw']=wtype
   savedata['slid']=slide
   savedata['sbin']=spcbinfac
   savedata['vers']=version

   filename=uniqfname(filename,'speca')                                   # Get the next available name of form filename(x).speca
   wfile = open(filename, 'wb')                                           # Open file to write to

   pickle.dump(savedata,wfile)                                            # Pickle the data (convert into bitstream) and dump to file
   wfile.close()                                                          # Close file

   return filename


#-----SRinR------------------------------------------------------------------------------------------------------------

def srinr(t,binning,domain,minv=None,maxv=None):

   '''Subrange in Range: Subrange Creator

   Description:

    Takes an array of ordered values and asks the user to select a subrange.  Converts the
    users entries (in the same units as the data) into data IDs and checks that they are valid and
    within the array.

   Inputs:

    t       -   LIST: An evenly spaced list of values.
    binning -   LIST: The spacing of the list of values.
    domain  - STRING: The name of the physical quantity represented by the data in t.

   Outputs:

    new_mn  -    INT: The array index of the subrange minimum.
    new_mx  -    INT: The array index of the subrange maximum.
    boolv   -   BOOL: A flag to denote whether range was clipped.

   -J.M.Court, 2015'''

   old_mn=0
   old_mx=len(t)

   try:
      if minv==None:
         new_mn=float(raw_input('Minimum '+domain+': '))                  # Fetch new min value from user
      else:
         new_mn=minv
      new_mn=len(t[t<new_mn]) 
   except:
      new_mn=old_mn                                                       # Treat garbage input as 'no change'

   try:
      if maxv==None:
         new_mx=float(raw_input('Maximum '+domain+': '))                  # Fetch new max value from user
      else:
         new_mx=maxv
      new_mx=len(t[t<=new_mx])-1      
   except:
      new_mx=old_mx                                                       # Treat garbage input as 'no change'

   if new_mn>=new_mx:
      print 'Invalid clipping!  Resetting to full.'
      boolv=False                                                         # Set the flag to propagate the fact no clipping took place
      new_mn=old_mn
      new_mx=old_mx
   else:
      boolv=True
   return new_mn,new_mx,boolv


#-----TNorm------------------------------------------------------------------------------------------------------------

def tnorm(t,res):

   '''Time-Normer

   Description:

    Takes an array of evenly-spaced time values and normalises them by subtracting the lowest value
    from each.  Then forces each time to equal an integer multiple of the pre-normalising resolution
    to deal with computational inaccuracies in the subtraction process.

   Inputs:

    t   -  LIST: The list of evenly spaced time-values.  They need not be ordered, but the lowest
                 value must be first.
    res - FLOAT: The resolution of the list t; entered by hand again to negate computational
                 error.

   Outputs:

    nt  -  LIST: The normalised time values.

   -J.M.Court, 2014'''

   t=np.array(t)
   sttime=t[0]
   for i in range(len(t)):
      t[i]=res*np.floor((t[i]-sttime)/float(res))                         # Rescaling time to start at 0, rounding to deal with errors incurred 
                                                                          # by subtraction of large numbers
   return t


#-----TokenLoc---------------------------------------------------------------------------------------------------------

def tokenloc(enigma,token):

   envals=[]
   entoks=[]

   for i in range(len(enigma)):

      if enigma[i][0]==',':

         et=enigma[i][1:]                                           # Extract the tokens and how many bytes are allocated to each in the datawords

      else:

         et=enigma[i]

      entoks.append(et[0])

      if et[0]=='C':

         tk=et.split('[')[1].split(']')[0].split(',')
         envals.append(int(enigma[i].split('{')[1]))

   if not (token in entoks):

      return None,None

   cloc=entoks.index(token)                                         # Fetch the 'channels' token
   r1=sum(envals[:cloc])
   r2=sum(envals[cloc+1:])

   return tk,(r1,r2)

#-----UniqFName--------------------------------------------------------------------------------------------------------

def uniqfname(filename,extension):

   '''Unique Filename Generator

   Description:

    When given the path to a file, ascertains whether a file already exists in that location.  If it
    does, this function finds the lowest value of n such that 'filename(n).extension' is unique and
    returns this new name.  If the file does not exist, returns 'filename.extension'.

   Inputs:

    filename  - STRING: The full path to the proposed file location, excluding the extension.
    extension - STRING: The extension of the proposed file location, excluding the leading '.'.

   Outputs:

    uniqname  - STRING: A string containing the best available unique filename.

   -J.M.Court, 2015'''

   filenamex=filename

   n=1

   while os.path.isfile(filenamex+'.'+extension):
      # print filename+'.'+extension+' already exists!                    # Can uncomment this for more verbosity.
      filenamex=filename+'('+str(n)+')'
      n+=1

   uniqname=filenamex+'.'+extension

   return uniqname


#-----VCRebin----------------------------------------------------------------------------------------------------------

@mjit()
def vcrebin(vecdata,bfac):

   '''Vector X-Rebin

   Description:

    Takes 1-Dimensional array of data and which is linearly binned and rebins it by a factor of the
    user's choosing, returning new data array.

   Inputs:

    vecdata   - ARRAY: some 2 dimensional array of values.

   Outputs:

    b_vecdata - ARRAY: the rebinned 2-dimensional data array

   -J.M.Court, 2015'''

   bfac=int(bfac)

   spx=len(vecdata)                                                       # x-dimension of new vector
   b_vecdata=np.zeros(spx)                                                # Create the new vector

   for i in range(spx):                                                   # For each freq row of fourgr:

      celltot=0                                                           # Create a running total to deposit into the vector cell

      for k in range(bfac):                                               # Calculate extent of current bin

         celltot+=vecdata[i*bfac+k]                                       # Sum all values that fall within the given bin

      b_vecdata[i]=celltot/bfac                                           # Divide the cell total by the bin multiplier to convert to a mean

   return b_vecdata



#-----XtrFilLoc--------------------------------------------------------------------------------------------------------

@mjit()
def xtrfilloc(filepath):

   '''Extract File Location

   Description:

    Given the relative path to a file, extracts the file name and it's location.

   Inputs:

    filepath - STRING: The path to a file

   Outputs:

    filename - STRING: The name of the file
    fileloca - STRING: The location of the file

   -J.M.Court, 2015'''

   filename=(filepath.split('/')[-1])                                     # Identify filename without directory

   if filename!=filepath:
      fileloca=filepath[:-len(filename)]                                  # Fetch directory
   else:
      fileloca=os.getcwd()

   if fileloca[0]!='/':
      fileloca=os.getcwd()+fileloca[1:]                                   # Dump current directory onto the front to make this absolute

   return filename,fileloca


