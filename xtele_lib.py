#!/usr/bin/python

# |----------------------------------------------------------------------|
# |------------------------------XTELE_LIB-------------------------------|
# |----------------------------------------------------------------------|

# A selection of useful functions which are placed here to reduce clutter in the other files of
# xte-l-extract.
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
#  CHRANGE   - clips a set of event data given to it to screen out photons outside of some energy range
#              given by the user.
#
#  EVMCHAN   - converts an RXTE channel ID into an E_125us_64M_0_1s DATAMODE channel range ID.
#
#  FOLDIFY   - takes a time series with its associated y-axis data and y-axis errors.  Folds this data
#              over a time period of the user's choosing, and returns them as the tuple x,y,y_error.
#
#  GETPCUS   - return the total number of PCUs which contributed photons to a given set of event words.
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
#  RMS_N     - takes the raw power spectrum output from the scipy FFT algorithm and normalises it using
#              (RMS/Mean)^2 normalisation.
#
#  SLPLOT    - plots an x-y line plot of two sets of data, and then below plots the same data on another
#              set of axes in log-log space.
#
#  SMFOLD    - uses foldify to fold data, and prints the factor by which the data's amplitude was flattened.
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
#

#-----Importing Modules------------------------------------------------------------------------------------------------

import os,cPickle
import pylab as pl
from numpy import array, arange, ceil, exp, floor, log10, mean, sqrt, zeros
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
      exit()


#-----Binify-----------------------------------------------------------------------------------------------------------

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
   return xb,yb,yeb


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
 

#-----ChRange----------------------------------------------------------------------------------------------------------

def chrange(data,low,high,datamode):

   '''Channel Ranger

   Description:

    Takes raw event or GoodXenon data from PCA on RXTE and removes all photons outside of an energy
    channel range given by the user.

   Inputs:

    data     - FITS DATA: the event data; raw photon arrival times and data words.
    low      -       INT: the lowest channel the user wants to consider.
    high     -       INT: the highest channel the user wants to consider.
    datamode -    STRING: the datamode in which the data was taken.

   Outputs:

    ch_data  - FITS DATA: the same as data, but with all photons outside of the given channel range
                          removed.

   -J.M.Court, 2015'''

   words=data.field(1)

   if low<=0 and high>=255:
      return data                                                         # Don't bother searching through if the user wants full range

   if datamode=='E_125us_64M_0_1s':
      low=evmchan(low)                                                    # Convert the channels given into range IDs
      high=evmchan(high)
      r=5,11                                                              # Identify where in the E_125 data word the channel is hidden
   elif datamode=='GoodXenon_2s':
      r=17,25                                                             # GoodXenon data contains the channels as written, no need to convert
   else:
      print datamode,'not yet supported, using full range!'
      return data

   words=array(boolval((words[:,r[0]:r[1]]).tolist()))

   mask1=(words>=low)
   mask2=(words<=high)

   ch_data=data[mask1&mask2]

   return ch_data


#-----EvMChan----------------------------------------------------------------------------------------------------------

def evmchan(chan):

   '''Event Mode Channel-Get

   Description:

    Converts a channel number into the ID of the range which contains that channel in E_125us_64M_0_1s
    data from PCA on RXTE.  This is the least interesting function I have ever written.

   Inputs:

    chan   - INT: the real channel number for PCA data from RXTE.

   Outputs:

    n_chan - INT: the ID of the range containing the relevant channel.

   -J.M.Court, 2015'''

   chan=int(chan)

   if chan<5:     n_chan=0                                                # This is just a list of ifs.  It checks if the value falls into each and, if not, carries on.
   elif chan<7:   n_chan=1
   elif chan<8:   n_chan=2
   elif chan<9:   n_chan=3
   elif chan<10:  n_chan=4
   elif chan<11:  n_chan=5
   elif chan<12:  n_chan=6
   elif chan<13:  n_chan=7
   elif chan<14:  n_chan=8
   elif chan<15:  n_chan=9
   elif chan<16:  n_chan=10
   elif chan<18:  n_chan=11
   elif chan<20:  n_chan=12
   elif chan<22:  n_chan=13
   elif chan<24:  n_chan=14
   elif chan<26:  n_chan=15
   elif chan<28:  n_chan=16
   elif chan<30:  n_chan=17
   elif chan<32:  n_chan=18
   elif chan<34:  n_chan=19
   elif chan<36:  n_chan=20
   elif chan<38:  n_chan=21
   elif chan<40:  n_chan=22
   elif chan<42:  n_chan=23
   elif chan<44:  n_chan=24
   elif chan<47:  n_chan=25
   elif chan<50:  n_chan=26
   elif chan<53:  n_chan=27
   elif chan<56:  n_chan=28
   elif chan<59:  n_chan=29
   elif chan<62:  n_chan=30
   elif chan<65:  n_chan=31
   elif chan<68:  n_chan=32
   elif chan<72:  n_chan=33
   elif chan<76:  n_chan=34
   elif chan<80:  n_chan=35
   elif chan<84:  n_chan=36
   elif chan<88:  n_chan=37
   elif chan<92:  n_chan=38
   elif chan<97:  n_chan=39
   elif chan<102: n_chan=40
   elif chan<107: n_chan=41
   elif chan<112: n_chan=42
   elif chan<117: n_chan=43
   elif chan<122: n_chan=44
   elif chan<127: n_chan=45
   elif chan<132: n_chan=46
   elif chan<138: n_chan=47
   elif chan<144: n_chan=48
   elif chan<150: n_chan=49
   elif chan<156: n_chan=50
   elif chan<162: n_chan=51
   elif chan<168: n_chan=52
   elif chan<175: n_chan=53
   elif chan<182: n_chan=54
   elif chan<189: n_chan=55
   elif chan<196: n_chan=56
   elif chan<203: n_chan=57
   elif chan<210: n_chan=58
   elif chan<218: n_chan=59
   elif chan<226: n_chan=60
   elif chan<234: n_chan=61
   elif chan<244: n_chan=62
   else:          n_chan=63

   return n_chan

#-----Foldify----------------------------------------------------------------------------------------------------------

def foldify(t,y,ye,period,binsize):                                       # Defining 'foldify' subscript

   '''Foldify

   Description:

    Takes a 2-dimensional set of data; given a user-input period in the t data, this function
    additively folds the y data over that period.

   Inputs:

    t       -  LIST: The t- or x-values of the two-dimensional data.  The period to be folded over
                     is a period in this dimension.
    y       -  LIST: The y-values of the two-dimensional data, must be the same length as x.
    ye      -  LIST: The errors associated with the y-values of the two-dimensional data, must be
                     the same length as x and y.
    period  - FLOAT: The period of the data, in the same units as the data's x-values.
    binsize - FLOAT: the size of the new x-axis bins in which to re-bin the data.

   Outputs:

    newt    -  LIST: A clipped version of the input t which now only corresponds to one period.
    newy    -  LIST: The folded y-values of the two-dimensional data for one period.
    newye   -  LIST: The errors associated with the folded y-values of the two-dimensional data for
                     one period.

   -J.M.Court, 2014'''

   ndat=len(t)                                                            # Calculating the number of data points to use

   rperiod=binsize*floor(period/binsize)                                  # Round the period to the nearest multiple of binsize

   newt=arange(0,period,binsize)                                          # Create new time array with same binning as initial data and range of 1 period
   newy=zeros(len(newt))                                                  # Set up equal length blank arrays for flux and flux error
   newye=zeros(len(newt))
   ndpy=zeros(len(newt))                                                  # Set up array to keep track of how many times data points were placed into each bin

   for i in range(ndat):                                                  # For all data points...
      fx=t[i]                                                             # Get the time
      fy=y[i]                                                             # Get the flux
      fye=(ye[i]**2)                                                      # Get the square of the flux error
      fx-=floor(fx/period)*period                                         # Subtract an integer number of periods from the timestamp so it is in the range [0,period]
      while fx<0:
         fx+=period                                                       # Bump up any values which, through computational error, get placed at a negative location
      fx=floor(fx/binsize)                                                # Convert the timestamp into an index within the folded arrays
      newy[fx]+=fy                                                        # Add flux and data into their new folded bins
      newye[fx]+=fye
      ndpy[fx]+=1                                                         # Add one to the total number of data points contributing to that bin

   for i in range(len(newt)):
      if ndpy[i]!=0:
         newy[i]=newy[i]/ndpy[i]                                          # Divide through by the number of data points thrown into each bin
         newye[i]=sqrt(newye[i])/ndpy[i]

   return newt, newy, newye


#-----Get PCUs---------------------------------------------------------------------------------------------------------

def getpcus(words,datamode):

   '''Get PCUs

   Description:

    If given a list of RXTE event words, and the data-mode in which their associated data is stored,
    returns the number of PCUs that contributed photons to the dataset.  Currently works for two
    datamodes: GoodXenon_2s and E_125us_64M_0_1s with the intention to eventually generalise it to
    any event data.

   Inputs:

    words    - 2D ARRAY: The array of event words, each of which is an array of True/False statements.
    datamode -   STRING: The DATAMODE in which the data to be analysed is stored.

   Outputs:

    pcus     -      INT: The number of PCUs active when the data was taken.  If a non-recognised event
                         word is given, 1 is returned instead along with an error message.

   -J.M.Court, 2015'''

   if datamode=='E_125us_64M_0_1s':
      r=1,4
   elif datamode=='GoodXenon_2s':
      r=7,10
   else:
      print datamode,'not yet supported!'
      return 1

   pcus=0

   words=(words[:,r[0]:r[1]]).tolist()

   if [False,False,False] in words: pcus+=1
   if [False,False,True] in words: pcus+=1
   if [False,True,False] in words: pcus+=1
   if [False,True,True] in words: pcus+=1
   if [True,False,False] in words: pcus+=1

   if pcus==0:
      print 'Error!  0 PCUs present in data!'
      return 1
   else:
      return pcus


#-----LBinify----------------------------------------------------------------------------------------------------------

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
   if denom==0: denom=1
   rms= (leahy-const)*(rate)/denom
   return rms


#-----LHConst----------------------------------------------------------------------------------------------------------

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

   olen=len(data)
   olen=int(4.0/5.0*olen)
   data=data[olen:]
   const=mean(data)
   if const>2.5 or const<1.5:
      print "WARNING: Could not find Leahy constant!"
      const=2
   return const


#-----MXRebin----------------------------------------------------------------------------------------------------------

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

   -J.M.Court,2015'''

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


#-----PlotdLd----------------------------------------------------------------------------------------------------------

def plotdld(filename,):

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

   -J.M.Court, 2015'''

   print 'Opening '+str(filename)

   readfile=open(filename,'rb')
   data=cPickle.load(readfile)                                            # Unpickle the .speca file

   times=array(data['time'])                                              # Unleash the beast! [open the file]
   counts=array(data['flux'])
   errors=array(data['errs'])
   binsize=data['bsiz']
   gti=data['gtis']
   mxpcus=data['pcus']
   bgest=data['bkgr']
   flavour=data['flav']

   bgpcu=bgest*mxpcus                                                     # Collect background * PCUs

   readfile.close()

   print ''

   print 'Power spectra taken over '+str(foures)+'s of data each'
   print str(sum(good))+'/'+str(len(spcdata))+' power spectra are good'
   if flavour!='':
      print "Flavour is '"+str(flavour)+"'"                               # Only print flavour if there is a flavour to print

   return times,counts,errors,binsize,gti,mxpcus,bgpcu,flavour


#-----PlotdSv----------------------------------------------------------------------------------------------------------

def plotdsv(filename,times,counts,binsize,gti,mxpcus,bgest,flavour):

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

   Outputs:

    [none]

   -J.M.Court, 2015'''

   savedata={}                                                            # Open library object to save in file

   savedata['time']=times                                                 # Dump each piece of data into an appropriate library element
   savedata['flux']=counts
   savedata['errs']=sqrt(abs(array(counts)))
   savedata['bsiz']=binsize
   savedata['gtis']=gti
   savedata['pcus']=mxpcus
   savedata['bkgr']=bgest
   savedata['flav']=flavour

   filename=uniqfname(filename,'plotd')                                   # Get the next available name of form filename(x).plotd
   wfile = open(filename, 'wb')                                           # Open file to write to

   cPickle.dump(savedata,wfile)                                           # Pickle the data (convert into bitstream) and dump to file

   print "PlotDemon file saved to "+filename

   wfile.close()                                                          # Close file


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


#-----SLPlot-----------------------------------------------------------------------------------------------------------

def slplot(x,y,ye,xlabel,ylabel,title,figid="",typ='both'):

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
    figid  - STRING: [optional] gives the Figure a name such that it can be easily manipulated or
                     closed outside of this function.
    typ    - STRING: [optional] user can input 'log' or 'lin' to only display plot of that type.

   Outputs:

    [none]

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
      pl.errorbar(x,y,ye)                                                 # Plot data
      pl.plot([x[0],x[-1]],[0,0])

   if typ in ('log','both'):

      if typ=='log':
         ax=pl.subplot(111)                                               # If 'log' passed as typ word, only create one subplot
      else:
         ax=pl.subplot(212)                                               # If 'both' passed as typ word, make 2nd of 2 subplots
      pl.xlabel(xlabel)
      pl.ylabel(ylabel)
      pl.title('Log '+str(title))
      pl.errorbar(x,abs(y),ye)                                            # Plot log-log data
      ax.set_xscale('log')
      ax.set_yscale('log')
      pl.grid(True,which="both")

   if typ in ('lin','log','both'):
      pl.show(block=False)                                                # Show both plots together
   else:
      'Invalid typ!  No plot shown.'                                      # Complain if none of 'lin', 'log' or 'both are given as typ word


#-----SmFold-----------------------------------------------------------------------------------------------------------

def smfold(t,y,ye,period,binsize,name=''):

   '''Smart Folder

   Description:

    Folds a two-dimensional set of data over a period in the first dimension using the foldify script
    which is also provided in xtele_lib.  Also calculates how much the peak-trough difference of the
    data has been compressed by the fold, and tells the user.

   Inputs:

    t       -   LIST: The t- or x-values of the two-dimensional data.  The period to be folded over
                      is a period in this dimension.
    y       -   LIST: The y-values of the two-dimensional data, must be the same length as x.
    ye      -   LIST: The errors associated with the y-values of the two-dimensional data, must be
                      the same length as x and y.
    period  -  FLOAT: The period of the data, in the same units as the data's x-values.
    binsize -  FLOAT: the size of the new x-axis bins in which to re-bin the data.
    name    - STRING: [optional] the name of the file to be folded.  This name will be used in text
                      outputs printed to screen.

   Outputs:

    newt    -  LIST: A clipped version of the input t which now only corresponds to one period.
    newy    -  LIST: The folded y-values of the two-dimensional data for one period.
    newye   -  LIST: The errors associated with the folded y-values of the two-dimensional data for
                     one period.

   -J.M.Court, 2015'''

   tn=tnorm(t,binsize)
   print 'Folding File '+name+'...'
   ptdiff=max(y)-min(y)                                                   # Flux range before folding
   newt,newy,newye=foldify(tn,y,ye,period,binsize)                        # Fold file 3 using 'foldify' in xtel_lib
   afdiff=max(newy)-min(newy)                                             # Flux range after folding
   print 'Flattened by '+str(100-afdiff/ptdiff*100)+'%'
   return newt,newy,newye


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
    phcts    -    INT: The total number of photons detected, overall.
    bg       -  ARRAY: An estimate of the count rate of the background flux during the full observation,
                       in counts per second per PCU, multiplied by the number of PCUs active in the
                       corresponding row of specdata.
    binsize  -  FLOAT: The size in seconds of the time bins used when converting event data into time
                       binned photon count data.
    foures   -  FLOAT: The length of time corresponding to a single row of spcdata, in seconds.
    flavour  - STRING: A useful bit of text to put on plots to help identify them later on.

   -J.M.Court, 2015'''

   print 'Opening '+str(filename)

   readfile=open(filename,'rb')
   data=cPickle.load(readfile)                                            # Unpickle the .speca file

   spcdata=data['data']                                                   # Unleash the beast! [open the file]
   good=array(data['good'])
   rates=array(data['rate'])
   phcts=data['phct']
   n_pcus=array(data['npcu'])
   binsize=data['bsiz']
   bgest=data['bkgr']
   foures=data['fres']
   flavour=data['flav']

   readfile.close()

   print ''

   print 'Power spectra taken over '+str(foures)+'s of data each'
   print str(sum(good))+'/'+str(len(spcdata))+' power spectra are good'
   if flavour!='':
      print "Flavour is '"+str(flavour)+"'"                               # Only print flavour if there is a flavour to print

   bg=n_pcus*bgest

   return spcdata,good,rates,phcts,bg,binsize,foures,flavour


#-----SpecaSv----------------------------------------------------------------------------------------------------------

def specasv(filename,spcdata,good,rates,phcts,npcus,binsize,bgest,foures,flavour):

   '''.Speca Save

   "Should've gone to SpecaSaver..."

   Description:

    Takes the input of the data products required to create a .speca file (to read with specangel)
    and creates a .speca file at a location given as the first input.

   Inputs:

    filename - STRING: The absolute or relative path to the location of the file that will be created.
    spcdata  -  ARRAY: An array, the elements of which are the non-normalised power spectra taken over
                       different equal-length time intervals of the photon count data being considered.
    good     -  ARRAY: A Boolean array with as many elements as spcdata has columns.  Its entries are
                       True unless the corresponding row in spcdata does not correspond to a valid time
                       within the GTIs of the photon count data.
    rates    -  ARRAY: An array of floats with as many elements as spcdata has columns.  Denotes the
                       total number of photon counts in the time interval represented by the corresponding
                       row in spcdata, divided by the width of a column in seconds.
    phcts    -    INT: The total number of photons detected, overall.
    npcus    -  ARRAY: An array of ints with as many elements as spcdata has columns.  States the number
                       of detectors that were active when the data represented by the corresponding row
                       of spcdata was recorded.
    binsize  -  FLOAT: The size in seconds of the time bins used when converting event data into time
                       binned photon count data.
    bgest    -  FLOAT: An estimate of the count rate of the background flux during the full observation,
                       in counts per second per PCU.
    foures   -  FLOAT: The length of time corresponding to a single row of spcdata, in seconds.
    flavour  - STRING: A useful bit of text to put on plots to help identify them later on.

   Outputs:

    [none]

   -J.M.Court, 2015'''

   savedata={}                                                            # Open library object to save in file

   savedata['data']=spcdata                                               # Dump each piece of data into an appropriate library element
   savedata['good']=good
   savedata['rate']=rates
   savedata['phct']=phcts
   savedata['npcu']=npcus
   savedata['bsiz']=binsize
   savedata['bkgr']=bgest
   savedata['fres']=foures
   savedata['flav']=flavour

   filename=uniqfname(filename,'speca')                                   # Get the next available name of form filename(x).speca
   wfile = open(filename, 'wb')                                           # Open file to write to

   cPickle.dump(savedata,wfile)                                           # Pickle the data (convert into bitstream) and dump to file

   print "SpecAngel file saved to "+filename

   wfile.close()                                                          # Close file


#-----SRinR------------------------------------------------------------------------------------------------------------

def srinr(t,binning,domain):

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
      new_mn=float(raw_input('Minimum '+domain+': '))                     # Fetch new min value from user
      new_mn=len(t[t<new_mn])  
   except:
      new_mn=old_mn                                                       # Treat garbage input as 'no change'

   try:
      new_mx=float(raw_input('Maximum '+domain+': '))                     # Fetch new max value from user
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

   -J.M.Court,2015'''

   filenamex=filename

   n=1

   while os.path.isfile(filenamex+'.'+extension):
      # print filename+'.'+extension+' already exists!                    # Can uncomment this for more verobisty.
      filenamex=filename+'('+str(n)+')'
      n+=1

   uniqname=filenamex+'.'+extension

   return uniqname


