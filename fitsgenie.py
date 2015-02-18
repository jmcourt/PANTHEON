#! /usr/bin/env python

# |----------------------------------------------------------------------|
# |-----------------------------FITS GENIE-------------------------------|
# |----------------------------------------------------------------------|

# Call as ./fitsgenie.py FILE1 [LCHAN] [HCHAN] [BINNING] [FOURIER RES] [BGEST] [FLAVOUR]

# Takes 1 FITS Event file and produces .speca and .plotd formatted products to be analysed by plotdemon
# and specangel.
#
# Arguments:
#
#  FILE1
#   The absolute path to the file to be used.
#
#  [LCHAN]
#   Optional: The lowest channel on the PCA instrument on RXTE which will be used to populate the data.
#   Default of 0 (minimum).
#
#  [HCHAN]
#   Optional: The highest channel on the PCA instrument on RXTE which will be used to populate the data.
#   Default of 255 (maximum).
#
#  [BINNING]
#   Optional: The size, in seconds, of bins into which data will be sorted.  Takes the value of the
#   time resolution of the data if not specified by the user.  Default of 2^-15s
#
#  [FOURIER RES]
#   Optional: The size of the individual time windows in which the data is to be split.  Fourier
#   spectra will be made of each of these windows.  Default of 128s.
#
#  [BGEST]
#   Optional: The approximate average background count rate during the observation in cts/s.  Default
#   of 30cts/s.
#
#  [FLAVOUR]
#   Optional: A useful bit of text to put on plots to help identify them later on.
#
#

#-----Importing Modules------------------------------------------------------------------------------------------------

import sys
import pan_lib as pan
from astropy.io import fits
from math import log
from numpy import arange, array, histogram, sqrt, zeros
from scipy.fftpack import fft
import pylab as pl


#-----User-set Parameters----------------------------------------------------------------------------------------------

ptdbinfac=16                                                              # To save space and time, the time bins for saved plotdemon data will be greater than the time
                                                                          # bins for the not-saved specangel data by this factor.  Must be power of 2.
usrmin=-13                                                                # The smallest time resolution to consider is 2^usrmin seconds


#-----Welcoming Header-------------------------------------------------------------------------------------------------

print ''
print '-------Running XTEGet: J.M.Court, 2015------'
print ''


#-----Checking Validity of Filename------------------------------------------------------------------------------------

args=sys.argv
pan.argcheck(args,2)                                                      # Must give at least 2 args (Filename and the function call)

filename=args[1]                                                          # Fetch file name from arguments


#-----Opening FITS file, identifying mission---------------------------------------------------------------------------

event=fits.open(filename)                                                 # Unleash the beast! [open the file]
mission=event[1].header['TELESCOP']                                       # Fetch the name of the telescope

if mission in ['XTE','SUZAKU']:
   print mission,'data detected!'
else:
   print mission,'data not yet supported!'
   pan.signoff()
   exit()

if mission == 'XTE' :
   etype='channel'                                                        # XTE requires an input of channel IDs
   escale=''
   escaleb=''
   import xtepan_lib as inst                                              # Import XTE extraction functions
elif mission == 'SUZAKU':
   etype='energy'                                                         # SUZAKU requires an input of raw energies
   escale='eV'
   escaleb=' (eV)'
   import szkpan_lib as inst                                              # Import SUZAKU extraction functions
else:
   print "This error shouldn't happen...  sorry 'bout that!"
   pan.signoff()
   exit()

obsdata=inst.getobs(event,event[1].header['DATAMODE'],filename)           # Fetch object and Obs_ID


#-----Checking validity of remaining inputs----------------------------------------------------------------------------

print 'Object =',obsdata[0]
print 'Obs_ID =',obsdata[1]
print ''

maxen=inst.maxen()                                                        # Get the value of the highest energy or channel for the instrument

print 'Inputs:'
print ''

if len(args)>2:
   lowc=int(args[2])                                                      # Collect minimum channel label from user
   print 'Min Channel =',lowc
else:
   try:
      lowc=int(raw_input("Minimum "+etype+escaleb+": "))
   except:
      lowc=0
      print "Using min "+etype+" of 0"+escale+"!"

if len(args)>3:
   highc=int(args[3])                                                     # Collect maximum channel label from user
   print 'Max Channel =',highc
else:
   try:
      highc=int(raw_input("Maximum "+etype+escaleb+": "))
   except:
      highc=maxen
      print "Using max "+etype+" of "+str(maxen)+escale+"!"

if lowc<0:    lowc=0                                                      # Force channels to be in range 0,255
if highc>maxen: highc=maxen

if lowc>highc:
   print 'Invalid '+etype+'"!  Aborting!'                                 # Abort if user gives lowc>highc
   pan.signoff()
   exit()

cs=str(int(lowc))+'-'+str(int(highc))

if len(args)>4:
   bszt=float(args[4])                                                    # Collect binsize from inputs if given, else ask user, else use resolution encoded in .fits file
   print 'Bin size (s)=',bszt
else:
   try:
      bszt=float(raw_input("Photon count bin-size (s): "))
   except:
      bszt=0
      print "Using max time resolution..."

if len(args)>5:
   foures=float(args[5])                                                  # Collect Fourier resolution from inputs if given, else ask user, else use 128s
   print 'Fourier Res.=',foures
else:
   try:
      foures=float(raw_input("Length of time per Fourier spectrum (s): "))
   except:
      foures=128
      print "Using 128s per spectrum..."

if len(args)>6:
   bgest=float(args[6])                                                   # Collect background estimate from inputs if given, else ask user, else use 30c/s
   print 'Background  =',bgest
else:
   try:
      bgest=float(raw_input("Estimate of background (c/s): "))
   except:
      bgest=30
      print "Using 30c/s background..."
   print ''

if len(args)>7:
   flavour=args[7]                                                        # Collect flavour if given, else flavourless
   print 'Flavour     =',flavour
else:
   flavour=''

print ''

#-----Masking data---------------------------------------------------------------------------------------

gti=inst.getgti(event)                                                    # Extract GTI indices
datas=inst.getdat(event)                                                  # Extract event data

print 'Discarding photons outside of '+etype+' range '+str(lowc)+escale+'-'+str(highc)+escale+'...'

datas=inst.discnev(datas)
olen=str(len(datas))
datas=inst.chrange(datas,lowc,highc,event[1].header['DATAMODE'])
tstart=inst.getini(event)

phcts=len(datas)
pcg=str(int(100*phcts/float(olen)))+'%'

print str(phcts)+'/'+olen+' photons fall within '+etype+' range ('+pcg+')!'

if phcts==0:
   print 'Aborting!'
   pan.signoff()
   exit()

print ''

times=inst.gettim(datas)                                                  # Extracting list of photon incident times as a separate object
pcwrds=inst.getwrd(datas)
sttim=times[0]
times=times-sttim


#-----Fetching Bin Size------------------------------------------------------------------------------------------------

bsz=inst.getbin(event)                                                    # Fetch 'Binning' as the time resolution of the data

if bszt>bsz:                                                              # If user enters a lower binning resolution than maximum, use that instead
   bsz=bszt

n=usrmin                                                                  # Rounding bsz to the nearest (greater) power of 2
while (2**n)<bsz:
   n+=1
bsz=2**n

print 'SpecAngel binsize rounded to 2^'+str(n)+'s ('+str(bsz)+'s)!'
print 'PlotDemon binsize rounded to 2^'+str(n+int(log(ptdbinfac,2)))+'s ('+str(bsz*ptdbinfac)+'s)!'


#-----Fetching Fourier Range Size--------------------------------------------------------------------------------------

if foures>max(times):
   foures=128

n=0                                                                       # Rounding foures to the nearest (greater) power of 2
while (2**n)<foures:
   n+=1
foures=2**n

print 'Fourier timeframe rounded to 2^'+str(n)+'s ('+str(foures)+'s)!'

print ''


#-----Rescaling GTI----------------------------------------------------------------------------------------------------

for j in range(len(gti)):
   gti[j]=gti[j][0]-sttim,gti[j][1]-sttim


#-----Setting up power spectra-----------------------------------------------------------------------------------------

ndat=int(max(times)/bsz)
datres=int(foures/bsz)                                                    # Work out how many data points corresponds to the user given time interval 'foures'
numstep=(ndat/datres)                                                     # Calculate how many intervals of 'datres' can be divided into the data length

print 'Analysing data...'
print ''

fourgrlin=[]                                                              # Set up matrix
bad=0                                                                     # Counter to count ranges which fall out of the GTIs
good=[]                                                                   # Array to keep track of which ranges were good
rates=[]                                                                  # Array of count rates to be populated
npcus=[]                                                                 
t=arange(0,foures+bsz,bsz)                                                # Setting up SpecAngel resolution time series per Fourier bin
tc=arange(0,foures+bsz*ptdbinfac,bsz*ptdbinfac)                           # Setting up PlotDemon coarse resolution time series per Fourier bin
ta=arange(0,(foures*numstep),bsz*ptdbinfac)                               # Setting up PlotDemon resolution full time series


#-----Populating power spectra-----------------------------------------------------------------------------------------

fullhist=[]                                                               # Create empty flux array to pass to plotdemon
fullerrs=[]
tcounts=0                                                                 # Initiate photon counter

for step in range(numstep):                                               ## For every [foures]s interval in the data:
   stpoint=step*foures                                                         #  Calculate the startpoint of the interval
   edpoint=stpoint+foures                                                      #  Calculate the endpoint of the interval

   in_gti=False                                                                #  Assume the subrange is not in the GTI

   for j in range(len(gti)):
      if gti[j][0]<=step*foures<(step+1)*foures<=gti[j][1]: in_gti=True        #  Change in_gti flag if this range is wholly within one GTI

   mask=times>=stpoint
   datrow=times[mask]                                                          #  Take all photons in the event data which occurred after the start point
   wrdrow=pcwrds[mask]
   mask=datrow<edpoint
   datrow=datrow[mask]                                                         #  Remove all photons which occurred after the end point
   wrdrow=wrdrow[mask]

   fc,null=histogram(datrow,tc+step*foures)                                    #  Coarsely bin this subrange of event data
   del null
   fullhist=fullhist+list(fc)
   fullerrs=fullerrs+list(sqrt(array(fc)))

   if in_gti:

      f,txis=histogram(datrow,t+step*foures)                                   #  Bin well this subrange of event data

      pcus=inst.getpcu(wrdrow,event[1].header['DATAMODE'])                     #  Count active PCUs by assuming any that recorded 0 events in the time period were inactive
      npcus.append(pcus)

      counts=sum(f)
      tcounts+=counts
      rates.append(float(counts)/foures)

      tsfdata=fft(f)                                                           #  Fourier transform the interval

      tsfdata=pan.leahyn(tsfdata,counts,datres)                                #  Normalise to Leahy Power
      good.append(True)                                                        #  Flag this column as good

   else:

      tsfdata=zeros(datres/2)
      npcus.append(0)
      rates.append(0)
      good.append(False)                                                       #  Flag this column as bad

   fourgrlin.append(tsfdata)                                                   #  Append the FT'd data to the matrix

   prog=step+1
   if (prog % 5)==0 or prog==numstep:
      print str(prog)+'/'+str(numstep)+' series analysed...'                   # Display progress every 5 series

pcg=str(int(100*tcounts/float(phcts)))+'%'

print ''
print str(tcounts)+'/'+str(phcts)+' photons in GTI ('+pcg+')!'

if tcounts==0:
   print 'Aborting!'
   pan.signoff()
   exit()


#-----Save .speca and .plotd files-------------------------------------------------------------------------------------

print ''
print 'Saving...'
print ''

filext=(filename.split('.')[-1])                                          # Identify file extension from the original filename
if filext!=filename:
   filename=filename[:-len(filext)-1]                                     # Remove file extension, if present

filename=filename+'_'+cs

pfilename=pan.plotdsv(filename,ta,fullhist,fullerrs,tstart,bsz*ptdbinfac,gti,max(npcus),bgest,flavour,cs,mission,obsdata)
print "PlotDemon file saved to "+pfilename
sfilename=pan.specasv(filename,fourgrlin,good,rates,tcounts,npcus,bsz,bgest,foures,flavour,cs,mission,obsdata)
print "SpecAngel file saved to "+sfilename

#-----Footer-----------------------------------------------------------------------------------------------------------

pan.signoff()


