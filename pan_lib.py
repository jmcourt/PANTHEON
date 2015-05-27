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
#  FLNCHECK  - checks to see whether a proposed input file has the correct file extension.
#
#  FOLDIFY   - takes a time series with its associated y-axis data and y-axis errors.  Folds this data
#              over a time period of the user's choosing, and returns them as the tuple x,y,y_error.
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
#  PDCOLEX   - extracts colours from a set of 2 or 3 lightcurves
#
#  PLOTDLD   - load and unpickle a .plotd file and extract its data.
#
#  PLOTDSV   - collect a selection of data products as a library, pickle it and save as a .plotd file.
#
#  RMS_N     - takes the raw power spectrum output from the scipy FFT algorithm and normalises it using
#              (RMS/Mean)^2 normalisation.
#
#  SIGNOFF   - prints an dividing line with some space.  That's all it does.
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
from matplotlib.ticker import ScalarFormatter
from numba import jit
from numpy import pi, sin, cos, array, arange, ceil, exp, floor, log10, mean, ones, sqrt, zeros
from numpy import append as npappend
from numpy import sum as npsum


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

@jit
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

   binlx=binsize*floor(x[0]/binsize)                                      # Initialising 'bin lowest x', or the lowest x value of the current bin
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
         yeb[-1]=(sqrt(yeb[-1]))/binct                                         #  Sqrt error and divide by bincount
         yb.append(y[xid])                                                     #  Append current y value into new array element
         yeb.append((ye[xid])**2)

         binct=1                                                               #  Reset bin count to 1

   yb[-1]=yb[-1]/binct                                                    ## Clean up final bin
   yeb[-1]=(sqrt(yeb[-1]))/binct
   return array(xb),array(yb),array(yeb)


#-----BoolVal----------------------------------------------------------------------------------------------------------

def boolval(data,reverse=True):

   '''Boolean Evaluator

   Description:

    Given a list of Boolean values, interprets them as binary and returns the corresponding integer.
    By default reads lists as having the highest first.

   Inputs:

    data    - LIST: a list of lists Boolean values.
    reverse - BOOL: [optional] defaults to True.  If set to False, then the Boolean strings will be
                    interpreted as binaries with the lowest value (1) first.

   Outputs:

    data    - LIST: the integer values represented by 'data' if its values are interpreted as binary.

   -J.M.Court, 2015'''

   keyr=range(len(data[0]))                                               # Set up the 'key range' to convert Bool list into int

   if reverse:
      keyr=keyr[::-1]                                                     # Reverse the key range

   mult=[1<<i for i in keyr]                                              # Mult is a list of all the powers of 2 from 2^0 to 2^(length of data)
   data=array(data)*mult
   data=npsum(data,axis=1)                                                # Multiply Boolean list by mult, sum per row

   return data
 

#-----Circfold---------------------------------------------------------------------------------------------------------

def circfold(x,y,t,pcoords=True):

   mult=(2*pi)/float(t)
   a=x*mult
   s=y*sin(a)
   c=y*cos(a)

   if pcoords:
      return a,y
   else:
      return s,c


#-----FlnCheck---------------------------------------------------------------------------------------------------------

def flncheck(filename,validext,cont=False):

   '''Filename Checker

   Description:

    Takes a filename and a string representing the expected file extension.  If the extension of the
    file does not match expectations, either kill the script or return 'False'.

   Inputs:

    filename - STRING: The filename to be checked.
    validext - STRING: The extenstion expected for the file (WITHOUT the leading '.')
    cont     -   BOOL: [Optional] if True, returns a value of False for an incorrect file extension.
                       If False (default), kills the script upon finding an incorrect file extension.

   Outputs:

    iscorr   -   BOOL: True if the file has the correct extension, False otherwise

   -J.M.Court, 2015'''

   flext=(filename.split('.')[-1])

   if flext != validext:
      if cont:
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
    phres   -  FLOAT: [optional] The resolution, in frequency space, at which data will be shown.
    name    - STRING: [optional] the name of the file to be folded.  This name will be used in text
                      outputs printed to screen.
    compr   -   BOOL: [optional] if set True, gives a fourth output which is the amount by which the
                      min-max difference of the data-set was compressed by folding.
    verb    -   BOOL: [optional] if set false, suppresses all non-error text output.

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
   phasx =arange(0,1,phres)
   phasy =zeros(npbins)
   phasye=zeros(npbins)
   ny=zeros(npbins)

   a=[];ax=[]
   b=[];bx=[]

   for i in range(len(y)):
      k=int(phases[i]*npbins)
      phasy[k]+=y[i]
      phasye[k]+=(ye[i]**2)
      ny[k]+=1

   phasy=phasy/ny
   phasye=sqrt(phasye)/ny

   afdiff=max(phasy)-min(phasy)                                            # Flux range after folding
   if verb:
      print 'Flattened by '+str(100-afdiff/ptdiff*100)+'%'

   if compr:
      return phasx,phasy,phasye,(afdiff/ptdiff)

   else:
      return phasx,phasy,phasye


#-----GTIMask----------------------------------------------------------------------------------------------------------

@jit
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

   times=array(times)
   mask=zeros(len(times),dtype=bool)                                      # Set up initial blank list of 'False'

   for gti in gtis:                                                       # For every GTI index:
      smask=(times>gti[0]) & (times<gti[1])                               # Create a submask which is the 'and'ed product of times>gti_start and times<gti_end
      mask=mask|smask                                                     # 'or' the submask with the main mask
   return mask
         

#-----LBinify----------------------------------------------------------------------------------------------------------

@jit
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
   hingep=int((hinge-x[0])/(x[1]-x[0]))

   lbin=log10(x[0])
   xb =10**(arange(lbin,log10(x[-1]),logres))                             # Setting up arrays to append binned values into
   yb =zeros(len(xb))
   yeb=zeros(len(xb))
   ct =zeros(len(xb))

   hingel=sum((xb)<=hinge)                                                # Getting the ID of the hinge-point in the log

   xbl=len(xb)

   lx=(log10(x)-lbin)/logres

   for i in range(hingel,xbl):

      lowid=int(((10**((i*logres)+lbin))-x[0])/(x[1]-x[0]))               # Calculate the ID of the lowest linear bin that corresponds to this log bin
      uppid=int(((10**(((i+1)*logres)+lbin))-x[0])/(x[1]-x[0]))           # Calculate the ID of the highest linear bin that corresponds to this log bin
      if uppid>lowid:
         yb[i]=mean(y[lowid:uppid])
         yeb[i]=(sqrt(sum(array(ye[lowid:uppid])**2)))/int(uppid-lowid)
      else:
         yb[i]=0                                                          # If no data found, error=power=0
         yeb[i]=0

   mask=x<hinge
   lmask=xb>hinge

   xf=npappend(x[mask],xb[lmask])
   yf=npappend(y[mask],yb[lmask])
   yef=npappend(ye[mask],yeb[lmask])

   return xf,yf,yef


#-----LeahyN-----------------------------------------------------------------------------------------------------------

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


#-----LHImprov---------------------------------------------------------------------------------------------------------

def lhimprov(spec):

   def rmsnoise(x,a,b):
      return x*a+b

   const_corr=optm.curve_fit(rmsnoise,arange(len(spec)),spec*arange(len(spec)))
   return const_corr[0][0]


#-----MXRebin----------------------------------------------------------------------------------------------------------

@jit
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
   b_good=zeros(spy,dtype=bool)                                           # array to label good columns
   b_spcdata=zeros([spx,spy])                                             # Create the new matrix
   b_spcerrs=zeros([spx,spy])                                             # Create the new error matrix

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
         b_spcerrs[i,j]=sqrt(errtot)/bfac

   b_xaxis=[]                                                             # Calculate new x-axis

   for i in range(int(len(xaxis)/bfac)):
      b_xaxis.append(xaxis[i*bfac])

   return b_spcdata,b_spcerrs,b_xaxis,b_good


#-----PDColEx----------------------------------------------------------------------------------------------------------

@jit
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
   fluxe=sqrt(ye[1]**2+ye[2]**2)                                          # Get flux error

   for i in range(1,3):                                                   # For the ith possible numerator band
      for j in range(1,3):                                                # For the jth possible denominator band
         if j!=i:                                                         # Prevents taking x/x colour
            ld=int(str(i)+str(j))
            col[ld]=(y[i]/y[j])                                           # Fetch colour
            cole[ld]=col[ld]*sqrt(((ye[i]/y[i])**2)+((ye[j]/y[j])**2))    # Fetch colour error

   return flux,fluxe,col,cole

@jit
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
   fluxe=sqrt(ye[1]**2+ye[2]**2+ye[3]**2)                                 # Get flux error

   for i in range(1,4):                                                   # For the ith possible numerator band
      for j in range(1,4):                                                # For the jth possible denominator band
         if j!=i:                                                         # Prevents taking x/x colour
            ld=int(str(i)+str(j))
            col[ld]=(y[i]/y[j])                                           # Fetch colour
            cole[ld]=col[ld]*sqrt(((ye[i]/y[i])**2)+((ye[j]/y[j])**2))    # Fetch colour error

   return flux,fluxe,col,cole


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
      print 'File not found!  Aborting!'
      signoff()
      exit()
   data=pickle.load(readfile)                                            # Unpickle the .plotd file

   times=array(data['time'])                                              # Unleash the beast! [extract the file]
   rates=array(data['flux'])
   errors=array(data['errs'])
   tstart=data['tstr']
   binsize=data['bsiz']
   gti=data['gtis']
   mxpcus=data['pcus']
   bgest=data['bkgr']
   flavour=data['flav']
   chanstr=data['chan']
   mission=data['miss']
   obsdata=data['obsd']
   version=data['vers']

   bgpcu=bgest*mxpcus                                                     # Collect background * PCUs

   readfile.close()

   return times,rates,errors,tstart,binsize,gti,mxpcus,bgpcu,flavour,chanstr,mission,obsdata,version


#-----PlotdSv----------------------------------------------------------------------------------------------------------

@jit
def plotdsv(filename,times,counts,errors,tstart,binsize,gti,mxpcus,bgest,flavour,chanstr,mission,obsdata,version):

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
   savedata['flux']=array(counts)/float(binsize)
   savedata['errs']=array(errors)/float(binsize)
   savedata['tstr']=tstart
   savedata['bsiz']=binsize
   savedata['gtis']=gti
   savedata['pcus']=mxpcus
   savedata['bkgr']=bgest
   savedata['flav']=flavour
   savedata['chan']=chanstr
   savedata['miss']=mission
   savedata['obsd']=obsdata
   savedata['vers']=version

   filename=uniqfname(filename,'plotd')                                   # Get the next available name of form filename(x).plotd
   wfile = open(filename, 'wb')                                           # Open file to write to

   pickle.dump(savedata,wfile)                                            # Pickle the data (convert into bitstream) and dump to file
   wfile.close()                                                          # Close file

   return filename


#-----RMS_N------------------------------------------------------------------------------------------------------------

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

#-----SignOff----------------------------------------------------------------------------------------------------------

def signoff():

   '''Sign Off

   Description:

    Prints an underline with some spaces.  That's all.

   -J.M.Court, 2015'''

   print ''
   print '------------------------------------------------'
   print ''


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
      pl.title('Log '+str(title))
      if errors:
         pl.errorbar(x,abs(y),ye,fmt='k')                                 # Plot data
      else:
         pl.plot(x,abs(y),'k')                                            # Plot log-log data
      ax.set_xscale('log')
      ax.set_yscale('log')
      pl.grid(True,which="both")

   if typ in ('lin','log','both'):
      pl.show(block=False)                                                # Show both plots together
   else:
      print 'Invalid typ!  No plot shown.'                                # Complain if none of 'lin', 'log' or 'both are given as typ word


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
   good=array(data['good'])
   rates=array(data['rate'])
   prates=array(data['prts'])
   trates=array(data['trts'])
   phcts=data['phct']
   n_pcus=array(data['npcu'])
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

@jit
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
      new_mn=old_mn
      new_mx=old_mx
   return new_mn,new_mx


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

   t=array(t)
   sttime=t[0]
   for i in range(len(t)):
      t[i]=res*floor((t[i]-sttime)/float(res))                            # Rescaling time to start at 0, rounding to deal with errors incurred by subtraction of large numbers
   return t


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

@jit
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
   b_vecdata=zeros(spx)                                                   # Create the new vector

   for i in range(spx):                                                   # For each freq row of fourgr:

      celltot=0                                                           # Create a running total to deposit into the vector cell

      for k in range(bfac):                                               # Calculate extent of current bin

         celltot+=vecdata[i*bfac+k]                                       # Sum all values that fall within the given bin

      b_vecdata[i]=celltot/bfac                                           # Divide the cell total by the bin multiplier to convert to a mean

   return b_vecdata



#-----XtrFilLoc--------------------------------------------------------------------------------------------------------

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


