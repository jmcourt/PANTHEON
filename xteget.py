#!/usr/bin/python

# |----------------------------------------------------------------------|
# |-----------------------------XTE GET SPEC-----------------------------|
# |----------------------------------------------------------------------|

# Call as ./xtegetspec.py FILE1 [LCHAN] [HCHAN] [BINNING] [FOURIER RES] [BGEST] [FLAVOUR]

# Takes 1 RXTE FITS Event file and produces an interactive spectrogram
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

#-----Importing Modules------------------------------------------------------------------------------------------------

import sys
import xtele_lib as xtl
from astropy.io import fits
from math import log
from numpy import arange, histogram, zeros
from scipy.fftpack import fft
import pylab as pl


#-----User-set Parameters----------------------------------------------------------------------------------------------

ptdbinfac=16                                                              # To save space and time, the time bins for saved plotdemon data will be greater than the time
                                                                          # bins for the not-saved specangel data by this factor.  Must be power of 2.


#-----Welcoming Header-------------------------------------------------------------------------------------------------

print ''
print '-------Running XTEGetSpec: J.M.Court, 2015------'
print ''


#-----Checking Validity of Arguments-----------------------------------------------------------------------------------

args=sys.argv
xtl.argcheck(args,2)                                                      # Must give at least 2 args (Filename and the function call)

filename=args[1]                                                          # Fetch file name from arguments

if len(args)>2:
   lowc=int(args[2])                                                      # Collect minimum channel label from user
else:
   try:
      lowc=int(raw_input("Minimum Channel: "))
   except:
      lowc=0
      print "Using min channel of 0!"

if len(args)>3:
   highc=int(args[3])                                                     # Collect maximum channel label from user
else:
   try:
      highc=int(raw_input("Maximum Channel: "))
   except:
      highc=0
      print "Using max channel of 255!"

if lowc<0:    lowc=0                                                      # Force channels to be in range 0,255
if highc>255: highc=255

if lowc>highc:
   print 'Invalid channels!  Aborting!'                                   # Abort if user gives lowc>highc
   print ''
   print '------------------------------------------------'
   exit()

if len(args)>4:
   bszt=float(args[4])                                                    # Collect binsize from inputs if given, else ask user, else use resolution encoded in .fits file
else:
   try:
      bszt=float(raw_input("Photon count bin-size (s): "))
   except:
      bszt=0
      print "Using max time resolution..."

if len(args)>5:
   foures=float(args[5])                                                  # Collect Fourier resolution from inputs if given, else ask user, else use 128s
else:
   try:
      foures=float(raw_input("Length of time per Fourier spectrum (s): "))
   except:
      foures=128
      print "Using 128s per spectrum..."

if len(args)>6:
   bgest=float(args[6])                                                   # Collect background estimate from inputs if given, else ask user, else use 30c/s
else:
   try:
      bgest=float(raw_input("Estimate of background (c/s): "))
   except:
      bgest=30
      print "Using 30c/s background..."
   print ''

if len(args)>7:
   flavour=args[7]                                                        # Collect flavour if given, else flavourless
else:
   flavour=''


#-----Opening FITS file, masking---------------------------------------------------------------------------------------

event=fits.open(filename)                                                 # Unleash the beast! [open the file]
gti=event[2].data                                                         # Extract GTI indices
datas=event[1].data                                                       # Extract event data

print 'Discarding photons outside of channel range '+str(lowc)+'-'+str(highc)+'...'

mask=datas['Event'][:,0]==True                                            # Creating a mask to obscure any data not labelled as photons
datas=datas[mask]                                                         # Applying the mask
olen=str(len(datas))
datas=xtl.chrange(datas,lowc,highc,event[1].header['DATAMODE'])

print str(len(datas))+'/'+olen+' photons fall within channel range!'
print ''

times=datas.field(0)                                                      # Extracting list of photon incident times as a separate object
words=datas.field(1)
sttim=times[0]
times=times-sttim


#-----Fetching Bin Size------------------------------------------------------------------------------------------------

bsz=event[1].header['TIMEDEL']                                            # Fetch 'Binning' as the time resolution of the data

if bszt>bsz:                                                              # If user enters a lower binning resolution than maximum, use that instead
   bsz=bszt

n=-15                                                                     # Rounding bsz to the nearest (greater) power of 2
while (2**n)<bsz:
   n+=1
bsz=2**n

print 'SpecAngel binsize rounded to 2^'+str(n)+'s ('+str(bsz)+'s)!'
print 'PlotDemon binsize rounded to 2^'+str(n+int(log(ptdbinfac,2)))+'s ('+str(bsz*ptdbinfac)+'s)!'


#-----Fetching Fourier Range Size--------------------------------------------------------------------------------------

if foures>max(times):
   foures=128

n=0                                                                       # Rounding fourez to the nearest (greater) power of 2
while (2**n)<foures:
   n+=1
foures=2**n

print 'Fourier timeframe rounded to 2^'+str(n)+'s ('+str(foures)+'s)!'

print ''


#-----Rescaling GTI----------------------------------------------------------------------------------------------------

for j in range(len(gti)):
   gti[j]=gti[j][0]-sttim,gti[j][1]-sttim


#-----Setting up power spectra-----------------------------------------------------------------------------------------

#print event[1].header['DATAMODE'],event[1].header['TEVTB2']              # Leaving this here for diagnostics; if something looks odd, check the event word

ndat=int(max(times)/bsz)
datres=int(foures/bsz)                                                    # Work out how many data points corresponds to the user given time interval 'foures'
numstep=(ndat/datres)                                                     # Calculate how many intervals of 'datres' can be divided into the data length

print 'Analysing data...'
print ''

fourgrlin=[]                                                              # Set up matrix
bad=0                                                                     # Counter to count ranges which fall out of the GTIs
good=[]                                                                   # Array to keep track of which ranges were good
phcts=[]                                                                  # Array of count rates to be populated
npcus=[]                                                                 
t=arange(0,foures+bsz,bsz)                                                # Setting up SpecAngel resolution time series per Fourier bin
tc=arange(0,foures+bsz*ptdbinfac,bsz*ptdbinfac)                           # Setting up PlotDemon coarse resolution time series per Fourier bin
ta=arange(0,(foures*numstep),bsz*ptdbinfac)                               # Setting up PlotDemon resolution full time series


#-----Populating power spectra-----------------------------------------------------------------------------------------

fullhist=[]
fulltxis=[]

for step in range(numstep):                                               ## For every [foures]s interval in the data:
   stpoint=step*foures                                                         #  Calculate the startpoint of the interval
   edpoint=stpoint+foures                                                      #  Calculate the endpoint of the interval

   in_gti=False                                                                #  Assume the subrange is not in the GTI

   for j in range(len(gti)):
      if gti[j][0]<=step*foures<(step+1)*foures<=gti[j][1]: in_gti=True        #  Change in_gti flag if this range is wholly within one GTI

   mask=times>=stpoint
   datrow=times[mask]                                                          #  Take all photons in the event data which occurred after the start point
   wrdrow=words[mask]
   mask=datrow<edpoint
   datrow=datrow[mask]                                                         #  Remove all photons which occurred after the end point
   wrdrow=wrdrow[mask]

   fc,null=histogram(datrow,tc+step*foures)                                    #  Coarsely bin this subrange of event data
   fullhist=fullhist+list(fc)

   if in_gti:

      f,txis=histogram(datrow,t+step*foures)                                   #  Bin well this subrange of event data

      pcus=xtl.getpcus(wrdrow,event[1].header['DATAMODE'])                     #  Count active PCUs by assuming any that recorded 0 events in the time period were inactive
      npcus.append(pcus)

      counts=sum(f)
      phcts.append(float(counts)/foures)

      tsfdata=fft(f)                                                           #  Fourier transform the interval

      tsfdata=xtl.leahyn(tsfdata,counts,datres)                                #  Normalise to Leahy Power
      good.append(True)                                                        #  Flag this column as good

   else:

      tsfdata=zeros(datres/2)
      npcus.append(0)
      phcts.append(0)
      good.append(False)                                                       #  Flag this column as bad

   fourgrlin.append(tsfdata)                                                   #  Append the FT'd data to the matrix

   prog=step+1
   if (prog % 5)==0 or prog==numstep:
      print str(prog)+'/'+str(numstep)+' series analysed...'                   # Display progress every 5 series


#-----Save .speca and .plotd files-------------------------------------------------------------------------------------

print ''
print 'Saving...'
print ''

filext=(filename.split('.')[-1])                                          # Identify file extension from the original filename
if filext!=filename:
   filename=filename[:-len(filext)-1]                                     # Remove file extension, if present

filename=filename+'_'+str(lowc)+'-'+str(highc)

xtl.plotdsv(filename,ta,fullhist,bsz,gti,max(npcus),bgest,flavour)
xtl.specasv(filename,fourgrlin,good,phcts,npcus,bsz,bgest,foures,flavour) # Save data, good array, counts array, #pcus array, background estimate, Fourier resolution,
                                                                          #  and flavour as a pickled binary library object; see specasv in xtele_lib

#-----Footer-----------------------------------------------------------------------------------------------------------

print ''
print '------------------------------------------------'
print ''


