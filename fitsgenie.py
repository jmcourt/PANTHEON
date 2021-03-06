#! /usr/bin/env python

# |----------------------------------------------------------------------|
# |-----------------------------FITS GENIE-------------------------------|
# |----------------------------------------------------------------------|

# Call as ./fitsgenie.py FILE1 PROD_REQ [LCHAN] [HCHAN] [BINNING] [FOURIER RES] [FOURIER SEP] [BGEST] [FLAVOUR]

# Takes 1 FITS Event file and produces .speca and .plotd formatted products to be analysed by plotdemon
# and specangel.
#
# Arguments:
#
#  FILE1
#   The absolute path to the file to be used.
#
#  PROD_REQ
#   The products requested by the user.  The following inputs are valid:
#      'spec','speca','s' will cause FITSGenie to produce only a .speca file as output
#      'plot','plotd','p' will cause FITSGenie to produce only a .plotd file as output
#      'both','all','b','a','sp','ps' will cause both files to be output
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
#  [FOURIER SEP]
#   Optional: The separation of the startpoints of individual time windows in which the data is to be
#   split.  Fourier spectra will be made of each of these windows.  Default of 128s.
#
#  [BGEST]
#   Optional: The approximate average background count rate during the observation in cts/s.  Default
#   of 30cts/s.
#
#  [FLAVOUR]
#   Optional: A useful bit of text to put on plots to help identify them later on.
#
#

#-----User-set Parameters----------------------------------------------------------------------------------------------

ptdbinfac=1                                                               # To save space and time, the time bins for saved plotdemon data will be greater than the time
                                                                          # bins for the not-saved specangel data by this factor.  Must be power of 2.
spcbinfac=4096                                                            # The binning factor for SpecAngel data to use when searching for data peaks and troughs
usrmin=-13                                                                # The smallest time resolution to consider is 2^usrmin seconds
version=6.1

#-----Welcoming Header-------------------------------------------------------------------------------------------------

print ''
print '-------Running FITSGenie: J.M.Court, 2015-------'
print ''


#-----Importing Modules------------------------------------------------------------------------------------------------

try:

   import sys
   import pan_lib as pan
   from astropy.io import fits


except ImportError:

   print 'Modules missing!  Aborting!'
   print ''
   print '------------------------------------------------'
   print ''
   exit()


#-----Checking Validity of Filename------------------------------------------------------------------------------------

args=sys.argv
pan.argcheck(args,1)                                                      # Must give at least 1 args (the function call)

if len(args)<2:
   filename=raw_input('Filename: ')
   print ''
else:
   filename=args[1]                                                       # Fetch file name from arguments

if len(args)<3:
   print ''
   print 'FitsGenie can produce [P]lotDemon files, [S]pecangel files or [B]oth.'
   print ''
   prod_req=raw_input('Select product(s): ')
   print ''
else:
   prod_req=args[2]                                                       # Fetch products request from arguments


#-----Identifying Products---------------------------------------------------------------------------------------------

if prod_req.lower() in ['spec','speca','s','both','all','a','b','ps','sp']:
   spec_on=True
   print '.speca File will be created!'
else:
   spec_on=False
   print '.speca File will NOT be created!'

if prod_req.lower() in ['plot','plotd','p','both','all','a','b','ps','sp']:
   plot_on=True
   print '.plotd File will be created!'
else:
   plot_on=False
   print '.plotd File will NOT be created!'

del prod_req
print ''

#-----Opening FITS file, identifying mission---------------------------------------------------------------------------

try:
   assert filename[-6:] not in ('.speca','.plotd')                        # Don't try to open plotd/speca files please...
   event=fits.open(filename)                                              # Unleash the beast! [open the file]
except:
   print 'Could not open file "'+filename+'"!'
   print 'Aborting!'
   pan.signoff()
   exit()

from math import log                                                      # Import remaining modules.  This is a slight speedup when running
                                                                          # a script to attempt to FITSgenie some valid files and some invalid
                                                                          # files
from numpy import arange, array, histogram, zeros
from numpy import sum as npsum
from scipy.fftpack import fft
import pylab as pl

try:
   mission=event[1].header['TELESCOP']                                    # Fetch the name of the telescope
except:
   print 'Could not identify mission!'
   print 'Aborting!'
   pan.signoff()
   exit()

if mission in ['XTE','SUZAKU','SWIFT']:
   print mission,'data detected!'
else:
   print mission,'data not yet supported!'
   pan.signoff()
   exit()

if mission == 'XTE' :
   etype='channel'                                                        # XTE requires an input of channel IDs
   escale=''
   escaleb=''

   #try:
   import xtepan_lib as inst                                           # Import XTE extraction functions
   #except:
   #   print 'XTE PANTHEON Library not found!  Aborting!'
   #   pan.signoff()
   #   exit()

elif mission == 'SUZAKU':
   etype='energy'                                                         # SUZAKU requires an input of raw energies
   escale='eV'
   escaleb=' (eV)'

   try:
      import szkpan_lib as inst                                           # Import SUZAKU extraction functions
   except:
      print 'Suzaku PANTHEON Library not found!  Aborting!'
      pan.signoff()
      exit()

elif mission == 'SWIFT':
   etype='energy'                                                         # SUZAKU requires an input of raw energies
   escale='eV'
   escaleb=' (eV)'

   try:
      import swfpan_lib as inst                                           # Import SUZAKU extraction functions
   except:
      print 'Swift PANTHEON Library not found!  Aborting!'
      pan.signoff()
      exit()

else:
   print "This error shouldn't happen...  sorry 'bout that!"
   pan.signoff()
   exit()

try:
   obsdata=inst.getobs(event,event[1].header['DATAMODE'],filename)        # Fetch object and Obs_ID
   print event[1].header['DATAMODE'],'format detected.'
except:
   print 'Could not identify DATAMODE!'
   print 'Aborting!'
   pan.signoff()
   exit()

if event[1].header['DATAMODE'][:2] in ['B_','SB']:
   spec_on=False
   bin_dat=True
   print 'No .speca file can be produced!'
   if not plot_on:
      print 'Aborting!'
      pan.signoff()
      exit()
else:
   bin_dat=False

#-----Checking validity of remaining inputs----------------------------------------------------------------------------

print 'Object =',obsdata[0]
print 'Obs_ID =',obsdata[1]
print ''

maxen=inst.maxen(event[1].header['DATAMODE'])                             # Get the value of the highest energy or channel for the instrument

print 'Inputs:'
print ''

if len(args)>3:
   lowc=int(args[3])                                                      # Collect minimum channel label from user
   print 'Min Channel =',lowc
else:
   try:
      lowc=int(raw_input("Minimum "+etype+escaleb+": "))
   except:
      lowc=0
      print "Using min "+etype+" of 0"+escale+"!"

if len(args)>4:
   highc=int(args[4])                                                     # Collect maximum channel label from user
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
   print 'Invalid '+etype+'!  Aborting!'                                  # Abort if user gives lowc>highc
   pan.signoff()
   exit()

cs=str(int(lowc))+'-'+str(int(highc))

if len(args)>5:
   bszt=float(args[5])                                                    # Collect binsize from inputs if given, else ask user, else use resolution encoded in .fits file
   print 'Bin size (s)=',bszt
else:
   try:
      bszt=float(raw_input("Photon count bin-size (s): "))
   except:
      bszt=0
      print "Using max time resolution..."

if len(args)>6:
   foures=float(args[6])                                                  # Collect Fourier resolution from inputs if given, else ask user, else use 128s
   print 'Fourier Res.=',foures
elif not spec_on:
   foures=16
else:
   try:
      foures=float(raw_input("Length of time per Fourier Window (s): "))
   except:
      foures=128
      print "Using 128s per spectrum..."

if len(args)>7:
   slide=float(args[7])                                                   # Collect Fourier resolution from inputs if given, else ask user, else use 128s
   print 'Fourier Sep.=',slide
elif not spec_on:
   slide=16
else:
   try:
      slide=float(raw_input("Separation of Fourier Windows (s): "))
   except:
      slide=foures
      print "Using "+str(slide)+"s per spectrum..."

if len(args)>8:
   bgest=float(args[8])                                                   # Collect background estimate from inputs if given, else ask user, else use 30c/s
   print 'Background  =',bgest
elif not spec_on:
   bgest=0
else:   
   try:
      bgest=float(raw_input("Estimate of background (c/s): "))
   except:
      bgest=30
      print "Using 30c/s background..."
   print ''

if len(args)>9:
   flavour=args[9]                                                        # Collect flavour if given, else flavourless
   print 'Flavour     =',flavour
else:
   flavour=''

print ''

wtype='Boxcar'                                                            # Setting all windows to BoxCar; will make transition to SpecAngel 4.0 smoother if this takes a value


#-----Masking data---------------------------------------------------------------------------------------

gti=inst.getgti(event)                                                    # Extract GTI indices
datas=inst.getdat(event)                                                  # Extract event data

print 'Discarding photons outside of '+etype+' range '+str(lowc)+escale+'-'+str(highc)+escale+'...'

#try:
datas=inst.discnev(datas,event[1].header['DATAMODE'])                     # Discarding non-events / reformatting XTE Binned data into a less awful structure
#except:
#   print 'Could not filter data!'
#   print 'Aborting!'
#   pan.signoff()
#   exit()

if event[1].header['DATAMODE'][:2] == 'B_':
   olen=str(npsum(array(datas)))
else:
   olen=str(len(datas))
datas=inst.chrange(datas,lowc,highc,event[1].header)
tstart=inst.getini(event)

if event[1].header['DATAMODE'][:2] == 'B_':
   phcts=npsum(datas)
else:
   phcts=len(datas)

if float(olen)==0:
   print 'No photons!  Aborting!'
   pan.signoff()
   exit()

pcg=str(int(100*phcts/float(olen)))+'%'

print str(phcts)+'/'+olen+' photons fall within '+etype+' range ('+pcg+')!'

if phcts==0:
   print 'Aborting!'
   pan.signoff()
   exit()

print ''


#-----Fetching Bin Size------------------------------------------------------------------------------------------------

bsz=inst.getbin(event,event[1].header['DATAMODE'])                        # Fetch 'Binning' as the time resolution of the data
ores=bsz

if bszt>bsz:                                                              # If user enters a lower binning resolution than maximum, use that instead
   bsz=bszt

n=usrmin                                                                  # Rounding bsz to the nearest (greater) power of 2
while (2**n)<bsz:
   n+=1
bsz=2**n

bsz=float(bsz)


#-----Fetching Time Axis-----------------------------------------------------------------------------------------------

print 'SpecAngel binsize rounded to 2^'+str(n)+'s ('+str(bsz)+'s)!'
print 'PlotDemon binsize rounded to 2^'+str(n+int(log(ptdbinfac,2)))+'s ('+str(bsz*ptdbinfac)+'s)!'
times,datas=inst.gettim(datas,event[1].data,tstart,ores,event[1].header['DATAMODE']) # Extracting list of photon incident times as a separate object

pcwrds=inst.getwrd(datas,event[1].header['DATAMODE'])
sttim=times[0]
times=times-sttim


#-----Fetching Fourier Range Size--------------------------------------------------------------------------------------

if foures>max(times):
   foures=128

if slide>max(times):
   slide=foures

slidelock=slide==foures                                                   # If foures=slide, lock them together

n=0                                                                       # Rounding foures to the nearest (greater) power of 2
while (2**n)<foures:
   n+=1
foures=2**n

if slidelock:
   slide=foures
else:
   plot_on=False
   print 'No .plotd file can be produced!'
   if not spec_on:
      print 'Aborting!'
      pan.signoff()
      exit()

print 'Fourier window length rounded to 2^'+str(n)+'s ('+str(foures)+'s)!'
print 'Fourier window separation rounded to',str(slide)+'s!'

print ''


#-----Rescaling GTI----------------------------------------------------------------------------------------------------

for j in pan.eqrange(gti):
   gti[j]=gti[j][0]-sttim,gti[j][1]-sttim


#-----Setting up power spectra-----------------------------------------------------------------------------------------

ndat=int(max(times)/bsz)
datres=int(foures/bsz)                                                    # Work out how many data points corresponds to the user given time interval 'foures'
stpres=int(slide/bsz)                                                     # Work out how many data points corresponds to the user given time separation 'slide'
numstep=(ndat/stpres)                                                     # Calculate how many intervals of 'datres' can be divided into the data length

while ((numstep-1)*slide)+foures>max(times):
   numstep-=1                                                             # Awful, lazy way to make sure the final bin doesn't exceed the data time limit

print 'Analysing data...'
print ''

fourgrlin=[]                                                              # Set up matrix
bad=0                                                                     # Counter to count ranges which fall out of the GTIs
good=[]                                                                   # Array to keep track of which ranges were good
prates=[]                                                                 # 'Peak Rates'
trates=[]                                                                 # 'Trough Rates'
rates=[]                                                                  # Array of count rates to be populated
npcus=[]                                                                 
t=arange(0,foures+bsz,bsz)                                                # Setting up SpecAngel resolution time series per Fourier bin
tp=arange(0,foures+bsz*spcbinfac,bsz*spcbinfac)                           # Setting up SpecAngel coarse resolution time series
tc=arange(0,foures+bsz*ptdbinfac,bsz*ptdbinfac)                           # Setting up PlotDemon coarse resolution time series per Fourier bin
ta=arange(0,(foures*numstep),bsz*ptdbinfac)                               # Setting up PlotDemon resolution full time series


#-----Populating power spectra-----------------------------------------------------------------------------------------

fullhist=[]                                                               # Create empty flux array to pass to plotdemon
fullerrs=[]
tcounts=0                                                                 # Initiate photon counter

pcus=None

if not bin_dat:

   for step in range(numstep):                                           ## For every [foures]s interval in the data:
      stpoint=step*slide                                                    #  Calculate the startpoint of the interval
      edpoint=stpoint+foures                                                #  Calculate the endpoint of the interval

      in_gti=False                                                          #  Assume the subrange is not in the GTI

      for j in pan.eqrange(gti):
         if gti[j][0]<=stpoint<edpoint<=gti[j][1]: in_gti=True              #  Change in_gti flag if this range is wholly within one GTI

      mask=times>=stpoint
      datrow=times[mask]                                                    #  Take all photons in the event data which occurred after the start point
      wrdrow=inst.getwrdrow(pcwrds,mask,event[1].header['DATAMODE'])
      mask=datrow<edpoint
      datrow=datrow[mask]                                                   #  Remove all photons which occurred after the end point
      wrdrow=inst.getwrdrow(wrdrow,mask,event[1].header['DATAMODE'])

      fc,null=histogram(datrow,tc+stpoint)                                  #  Coarsely bin this subrange of event data
      fp,null=histogram(datrow,tp+stpoint)                                  #  Very Coarsely bin this subrange of event data

      del null

      fullhist=fullhist+list(fc)
      fullerrs=fullerrs+list((array(fc)**0.5))

      if in_gti:

         f,txis=histogram(datrow,t+stpoint)                                 #  Bin well this subrange of event data

         pcus=inst.getpcu(wrdrow,event[1].header,t_pcus=pcus)               #  Count active PCUs by assuming any that recorded 0 events in the time period were inactive
         npcus.append(pcus)

         counts=sum(f)
         peak=max(fp)
         trough=min(fp)

         rates.append(float(counts)/foures)
         prates.append(float(peak)*datres/(foures*spcbinfac))
         trates.append(float(trough)*datres/(foures*spcbinfac))
         tcounts+=counts

         tsfdata=fft(f)                                                     #  Fourier transform the interval

         tsfdata=pan.leahyn(tsfdata,counts,datres)                          #  Normalise to Leahy Power
         good.append(True)                                                  #  Flag this column as good

      else:

         tsfdata=zeros(datres/2)
         npcus.append(0)
         rates.append(0)
         prates.append(0)
         trates.append(0)
         good.append(False)                                                 #  Flag this column as bad
  
      fourgrlin.append(tsfdata)                                             #  Append the FT'd data to the matrix

      prog=step+1
      if (prog % 5)==0 or prog==numstep:
         print str(prog)+'/'+str(numstep)+' series analysed...'             # Display progress every 5 series

   pcg=str(int(100*tcounts/float(phcts)))+'%'

   print ''
   print str(tcounts)+'/'+str(phcts)+' ('+pcg+') photons in GTI '+str((gti[0][0],gti[-1][1]))+'!'

   if tcounts==0:
      print 'Aborting!'
      pan.signoff()
      exit()

else:                                                                        # Not doing Spectra for Binned data just yet...

   print 'Number of PCUs unknown!'
   npcus=[int(raw_input('Number of Active PCUS: '))]                         # Ask the user how many there are
   ta,fullhist,fullerrs=pan.binify(times,datas/ores*bsz*ptdbinfac,(datas**0.5)/ores*bsz*ptdbinfac,bsz)


#-----Save .speca and .plotd files-------------------------------------------------------------------------------------

print ''
print 'Saving...'
print ''

filext=(filename.split('.')[-1])                                          # Identify file extension from the original filename
if filext!=filename:
   print filename
   tfilename=filename[:-len(filext)-1]                                    # Remove file extension, if present
   if tfilename[-1]!='.':                                                 # Saving extensionless files with .. in the path name breaks without this
      filename=tfilename

filename=filename+'_'+cs+'_'+str(bsz)+'s'

if plot_on:
   pfilename=pan.plotdsv(filename,ta,array(fullhist)/bsz,array(fullerrs)/bsz,tstart,bsz*ptdbinfac,gti,max(npcus),bgest,'False',None,flavour,cs,mission,obsdata,version)
   print "PlotDemon file saved to "+pfilename
else:
   print "PlotDemon file not saved."

if spec_on:
   sfilename=pan.specasv(filename,fourgrlin,good,rates,prates,trates,tcounts,max(npcus),bsz,bgest,foures,flavour,cs,mission,obsdata,wtype,slide,spcbinfac,version)
   print "SpecAngel file saved to "+sfilename
else:
   print "SpecAngel file not saved."


#-----Footer-----------------------------------------------------------------------------------------------------------

pan.signoff()


