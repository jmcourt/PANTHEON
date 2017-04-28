#! /usr/bin/env python

# |----------------------------------------------------------------------|
# |------------------------------PLOT DEMON------------------------------|
# |----------------------------------------------------------------------|

# Call as ./plotdemon.py FILE1 [FILE2] [FILE3] BINNING

# Takes 1-3 .plotd files and plots relevant astrometric plots
#
# Arguments:
#
#  FILE1
#   The absolute path to the first file to be used (generally the lowest energy band)
#
#  [FILE2]
#   The absolute path to the second file to be used
#
#  [FILE3]
#   The absolute path to the third file to be used (generally the highest energy band)
#
#  [BINNING]
#   Optional: the size, in seconds, of bins into which data will be sorted.
#
#

#-----User-set Parameters----------------------------------------------------------------------------------------------

minbin=0.0078125                                                          # The minimum bin size the code is allowed to attempt to use.  This can prevent long hang-ups
version=4.3                                                               # The version of PlotDemon
cbin=32.0                                                                 # The number of bins to use when calculating inhomonogeneity in circfold

#-----Welcoming Header-------------------------------------------------------------------------------------------------

print ''
print '-------Running Plot Demon: J.M.Court, 2014------'
print ''


#-----Importing Modules------------------------------------------------------------------------------------------------

try:

   import sys,os,imp
   import pylab as pl
   import pan_lib as pan
   import numpy as np
   import scipy.signal as sig
   import scipy.interpolate as ipo
   from math import pi


except ImportError:

   print 'Modules missing!  Aborting!'
   print ''
   print '------------------------------------------------'
   print ''
   exit()

try:
   imp.find_module('PyAstronomy')                                         # Check if PyAstronomy exists
   module_pyastro=True
except ImportError:
   module_pyastro=False


#-----Opening Files----------------------------------------------------------------------------------------------------

args=sys.argv                                                             # Fetching arguments; softest energy band first please
pan.argcheck(args,2)

try:
   float(args[-1])                                                        # If the final argument can be converted to integer, assume user intends it as a binning 
   isbininp=True                                                          # "IS BINsize given as an INPut?"
except:
   isbininp=False

nfiles=len(args)-isbininp-1                                               # Fetch number of infiles (total args minus one or two iff binsize given)
if nfiles>3: nfiles=3

file1=args[1]
pan.filenamecheck(file1,'plotd')

ch={}                                                                     # Save channel info in a library

bg1=0
bg2=0
bg3=0

print 'Opening',file1                                                     # Opening file 1
x1r,y1r,ye1r,tst1,bsz1,gti,pcus1,bg1,bsub1,bdata1,flv1,ch[1],mis1,obsd1,v1=pan.plotdld(file1)
y1r=y1r/float(pcus1)                                                      # Normalising flux by dividing by the number of active PCUs and the binsize
ye1r=ye1r/float(pcus1)

flavour=flv1
if flavour=='':
   qflav=''
else:
   qflav=' "'+flavour+'"'

if nfiles>1:
   file2=args[2]
   if pan.filenamecheck(file2,'plotd'):
      print 'Opening',file2                                               # Opening file 2
      x2r,y2r,ye2r,tst2,bsz2,gti2,pcus2,bg2,bsub2,bdata2,flv2,ch[2],mis2,obsd2,v2=pan.plotdld(file2)
      y2r=y2r/float(pcus2)                                                # Normalising flux by dividing by the number of active PCUs and the binsize
      ye2r=ye2r/float(pcus2)
   else:
      nfiles=1
      print 'Warning: File 2 of incorrect file type!'
      print 'Only loading File 1.'
else: x2r=y2r=ye2r=tst2=bsz2=None

if nfiles>2:
   file3=args[3]
   if pan.filenamecheck(file2,'plotd'):
      print 'Opening',file3                                               # Opening file 3
      x3r,y3r,ye3r,tst3,bsz3,gti3,pcus3,bg3,bsub3,bdata3,flv3,ch[3],mis3,obsd3,v3=pan.plotdld(file3)
      y3r=y3r/float(pcus3)                                                # Normalising flux by dividing by the number of active PCUs and the binsize
      ye3r=ye3r/float(pcus3)
   else:
      nfiles=2
      print 'Warning: File 3 of incorrect file type!'
      print 'Only loading Files 1 & 2.'
else: x3r=y3r=ye3r=tst3=bsz3=None

bg=(bg1+bg2+bg3)/float(pcus1)

xit1=x1r[-1]
if nfiles>1:
   xit2=x2r[-1]
else:
   xit2=None
if nfiles==3:
   xit3=x3r[-1]
else:
   xit3=None
oet=max(xit1,xit2,xit3)                                                   # Fetch the observation end time

mint=0                                                                    # Save original start and endpoints for use in clipping
maxt=oet

if nfiles>1:                                                              # Checking that start-times of files 1 & 2 match
   if tst1!=tst2:
      if tst1>tst2:
         while x1r[0]+tst1>x2r[0]+tst2:                                   # Hack data off of the start of file 2 until its startpoint matches file 1
            if len(x2r)==0:
               print 'Times domains for files 1 & 2 do not overlap!  Aborting!'
               pan.signoff()
               exit()
            x2r=np.delete(x2r,0)
            y2r=np.delete(y2r,0)
            ye2r=np.delete(ye2r,0)
         if tst1+x1r[0]!=tst2+x2r[0]:
            print 'Starting times for files 1 & 2 do not match!  Aborting!'# If this butchering overshoots, give up
            pan.signoff()
            exit()
         else:
            tst2+=x2r[0]                                                  # Amend new start time
            x2r=x2r-x2r[0]
      else:
         while x2r[0]+tst2>x1r[0]+tst1:                                   # Or Hack data off of the start of file 1 until its startpoint matches file 2
            if len(x1r)==0:
               print 'Times domains for files 1 & 2 do not overlap!  Aborting!'
               pan.signoff()
               exit()
            x1r=np.delete(x1r,0)
            y1r=np.delete(y1r,0)
            ye1r=np.delete(ye1r,0)
         if tst1+x1r[0]!=tst2+x2r[0]:
            print 'Starting times for files 1 & 2 do not match!  Aborting!'
            pan.signoff()
            exit()
         else:
            tst1+=x1r[0]
            x1r=x1r-x1r[0]


if nfiles>2:                                                              # Checking that start-times of files 1 & 3 match (and thus 2 & 3 also match)
   if tst1!=tst3:
      if tst1>tst3:
         while x1r[0]+tst1>x3r[0]+tst3:                                   # Hack data off of the start of file 3 until its startpoint matches file 1
            if len(x3r)==0:
               print 'Times domains for files 1 & 3 do not overlap!  Aborting!'
               pan.signoff()
               exit()
            x3r=np.delete(x3r,0)
            y3r=np.delete(y3r,0)
            ye3r=np.delete(ye3r,0)
         if tst1+x1r[0]!=tst3+x3r[0]:
            print 'Starting times for files 1 & 3 do not match!  Aborting!'# If this butchering overshoots, give up
            pan.signoff()
            exit()
         else:
            tst3+=x3r[0]                                                  # Amend new start time
            x3r=x3r-x3r[0]
      else:
         while x3r[0]+tst3>x1r[0]+tst1:                                   # Or Hack data off of the start of files 1 & 2 until their startpointa matches file 3
            if len(x1r)==0:
               print 'Times domains for files 1 & 3 do not overlap!  Aborting!'
               pan.signoff()
               exit()
            x1r=np.delete(x1r,0)
            y1r=np.delete(y1r,0)
            ye1r=np.delete(ye1r,0)
            x2r=np.delete(x2r,0)
            y2r=np.delete(y2r,0)
            ye2r=np.delete(ye2r,0)
         if tst1+x1r[0]!=tst3+x3r[0]:
            print 'Starting times for files 1 & 3 do not match!  Aborting!'
            pan.signoff()
            exit()
         else:
            tst1+=x1r[0]
            x1r=x1r-x1r[0]


#-----Binning----------------------------------------------------------------------------------------------------------

if isbininp:
   binning=float(args[-1])                                                # Collect binsize input if given
else:
   while True:                                                            # Keep asking until a good response is given
      try:
         binning=float(raw_input("Enter bin size (s): "))                 # Ask for binsize in dialogue box
         break
      except:
         print 'Invalid bin size input!'

if minbin>max(binning,bsz1,bsz2,bsz3):
   print 'Warning!  User-entered bin is smaller than the minimum!'
   print 'Minimum can be changed in the user-input section of this code'

binning=max(binning,bsz1,bsz2,bsz3,minbin)                                # Prevent overbinning by setting minimum binning to the maximum of the binnings of the files

print ''
print 'Bin size='+str(binning)+'s'  
print 'Binning File 1...'

x1,y1,ye1=pan.binify(x1r,y1r,ye1r,binning)                                # Bin File 1 using 'binify' in pan_lib
if nfiles>1:
   print 'Binning File 2...'
   x2,y2,ye2=pan.binify(x2r,y2r,ye2r,binning)                             # Bin File 2 using 'binify' in pan_lib
   if nfiles>2:
      print 'Binning File 3...'
      x3,y3,ye3=pan.binify(x3r,y3r,ye3r,binning)                          # Bin File 3 using 'binify' in pan_lib

print 'Binning complete!'
print ''

wrongsize=False

x3l=len(x1)                                                               # Sloppy fix to make this work for 2 or 3 mismatched files

#-----Force file lengths to match--------------------------------------------------------------------------------------

if nfiles>1:                                                              # Checking file lengths match
   if len(x1)!=len(x2):
      print 'Warning!  Files 1&2 of different lengths!'
      wrongsize=True

if nfiles==3:
   if len(x1)!=len(x3):
      print 'Warning!  Files 1&3 of different lengths!'
      wrongsize=True
      x3l=len(x3)

if wrongsize:                                                             # Forcing file lengths to match if possible
   print 'Attempting to crop files...'
   mindex=min(len(x1),len(x2),x3l)-1
   if x1[mindex]!=x2[mindex]:
      print 'Cannot crop, aborting!'
      pan.signoff()
      exit()
   if nfiles==3:
      if x1[mindex]!=x3[mindex]:
         print 'Cannot crop, aborting!'
         pan.signoff()
         exit()
   mindex+=1
   x1=x1[:mindex]
   y1=y1[:mindex]
   ye1=ye1[:mindex]
   x2=x2[:mindex]
   y2=y2[:mindex]
   ye2=ye2[:mindex]
   if nfiles==3:
      x3=x3[:mindex]
      y3=y3[:mindex]
      ye3=ye3[:mindex]
   print 'Cropped succesfully!'
   print ''
   
#-----Fetch GTI Mask---------------------------------------------------------------------------------------------------

print 'Fetching GTI mask...'

gmask=pan.gtimask(x1,gti)                                                 # A mask to blank values that fall outside of the GTIs
if nfiles>1:
   gmask2=pan.gtimask(x2,gti2)
   gmask=gmask&gmask2
if nfiles>2:
   gmask3=pan.gtimask(x3,gti3)
   gmask=gmask&gmask3
print str(int(100*sum(gmask)/len(gmask)))+'% of data within GTI!'
print ''


#-----Fetch Colours----------------------------------------------------------------------------------------------------

def colorget(verbose=True):                                               # Define colorget to easily re-obtain colours if base data is modified
   if verbose:
      print 'Analysing Data...'
   times=x1[gmask]
   timese=np.zeros(len(times))
   ys={}
   col={}
   cole={}
   if nfiles==1:                                                          # If only one file given, flux and flux_error are just the flux and error of this one file
      flux=y1[gmask]                                                      # Use gmask to clip out the areas outside of GTI
      fluxe=ye1[gmask]
   elif nfiles==2:
      flux,fluxe,ys,yes,col,cole=pan.pdcolex2(y1,y2,ye1,ye2,gmask)        # Get 2/1 and 1/2 colour information using PDColEx in pan_lib
   elif nfiles==3:
      flux,fluxe,ys,yes,col,cole=pan.pdcolex3(y1,y2,y3,ye1,ye2,ye3,gmask) # Get ALL colour values with 3D PDColEx
   else:
      print 'Error!  Too much data somehow.'                              # This warning should never come up...
      pan.signoff()
      exit()
   return times,timese,flux,fluxe,ys,yes,col,cole

times,timese,flux,fluxe,ys,yes,col,cole=colorget()                        # Use colorget

print 'Done!'
print ''


#-----Setting up plot environment--------------------------------------------------------------------------------------

show_block=False                                                          # Do not force plots to stay open by default
plotopt=''
es=True                                                                   # Options to keep track of what form the data is in.  'es': with error bars.
cs=False                                                                  # 'cs' with colour key
ls=False                                                                  # 'ls' with delineation
folded=False                                                              # 'folded' has been folded over some period
saveplots=False

def doplot(x,xe,y,ye,ovr=False,ft='-k',per2=False):                       # Defining short function to determine whether errorbars are needed on the fly
                                                                          # 'ovr' allows to override colour and line options, so lightcurves can be made differently

   if ovr: formst=ft                                                      # If override given, accept input format; if none given, just plot lines
   elif ls: formst='-ok'                                                  # If deLineate mode on, connect points with lines and mark points
   else: formst='ok'                                                      # If neither deLineate nor override on, just plot points.  No lines here, buddy

   if ls and not ovr:
      plotx=np.append(x,x[0])
      ploty=np.append(y,y[0])
      plotxe=np.append(xe,xe[0])
      plotye=np.append(ye,ye[0])

   elif per2 and ovr:
      plotx=np.append(x,x+1.0)
      ploty=np.append(y,y)
      plotxe=np.append(xe,xe)
      plotye=np.append(ye,ye)      

   else:
      plotx=x
      ploty=y
      plotxe=xe
      plotye=ye
   if cs and not ovr:                                                     # If coloured mode on, colour first 5 data points unless override given
      if len(x)<5:                                                        # Abort if less than 5 data points present
         print 'Not enough data to colour!'
      else:
         pl.plot(x[0],y[0],'or',zorder=1)                                 # Plot a round marker over each of the first five points with colour ascending red->blue
         pl.plot(x[1],y[1],'oy',zorder=2)
         pl.plot(x[2],y[2],'og',zorder=3)
         pl.plot(x[3],y[3],'oc',zorder=4)
         pl.plot(x[4],y[4],'ob',zorder=5)
   if es:
      pl.errorbar(plotx,ploty,xerr=plotxe,yerr=plotye,fmt=formst,zorder=0)# Plot errorbar plot if errors turned on
   else:
      pl.plot(plotx,ploty,formst,zorder=0)                                # Else plot regular graph


def plot_save(saveplots,show_block):                                      # Add a function to redirect all show calls to savefigs if toggled
   if saveplots:
      pl.savefig(raw_input('Save plot as: '))
      print 'Plot saved!'
   else:
      pl.show(block=show_block)

def burstplot(key,text,units):

   if bursts==None:
      print 'No burst data to plot!  Run "burst get" first!'
      return

   print 'Plotting Histogram of Burst '+text+'...'

   pl.figure()
   pl.hist(bursts[key],bins=np.arange(min(bursts[key]),max(bursts[key]),(max(bursts[key])-min(bursts[key]))/21.0))
   pl.xlabel(text,'('+units+')')
   pl.ylabel('Frequency')
   pl.title('Histogram of Burst '+text)
   plot_save(saveplots,show_block)

fldtxt=''
bursts=None
burst_alg='cubic spline'

flux_axis=r'Flux (cts s$^{-1}$ PCU$^{-1}$)'
time_n=0                                                                  # Normalise dump time


#-----User Menu--------------------------------------------------------------------------------------------------------

def give_inst():                                                          # Define printing this list of instructions as a function
   print 'COMMANDS: Enter a command to manipulate data.'
   print ''
   print 'DATA:'
   print '* "rebin" to reset the data and load it with a different binning.'
   print '* "clip" to clip the data.'
   print '* "mask" to remove a range of data.'
   print '* "rms" to return the fractional rms of the data.'
   print '* "fold" to fold data over a period of your choosing'+(' (requires PyAstronomy module!)' if not module_pyastro else '')+'.'
   print '* "autofold" to automatically seek a period over which to fold data'+(' (requires PyAstronomy module!)' if not module_pyastro else '')+'.'
   print '* "circfold" to circularly fold data over a period of your choosing.'
   print '* "varifold" to fold over a non-constant period using an algorithm optimised for high-amplitude quasi-periodic flares.'
   print ''
   print '1+ DATASET PLOTS:'
   print '* "lc" to plot a simple graph of flux over time.'
   print '* "bg" to plot background over time, if background has been estimated for these files.'
   print '* "animate" to create an animation of the lightcurve as the binning is increased.'
   print '* "circanim" to create an animation of the lightcurve circularly folded as the period is increased.'
   print '* "lombscargle" to create a Lomb-Scargle periodogram of the lightcurve.'
   print '* "autocor" to plot the auto-correlation function.'
   if nfiles>1:                                                           # Only display 2-data-set instructions if 2+ datasets given
      print ''
      print '2+ DATASET PLOTS:'
      print '* "hardness21" to plot a hardness/time diagram of file2/file1 colour over time.'
      print '* "hardness12" to plot a hardness/time diagram of file1/file2 colour over time.'
      print '* "hid21" to plot a hardness-intensity diagram of file2/file1 colour against total flux.'
      print '* "hid12" to plot a hardness-intensity diagram of file1/file2 colour against total flux.'
      print '* "col21" to plot file2/file1 colour against time.'
      print '* "col12" to plot file1/file2 colour against time.'
      print '* "band" to plot the lightcurve of a single energy band.'
      print '* "bands" to plot lightcurves of all bands on adjacent axes.'
      print '* "xbands" to plot lightcurves of all bands on the same axes.'
      print '* "compbands21" to plot lightcurves of bands 2 and 1 against each other.'
      print '* "crosscor21" to plot the cross-correlation function of band 1 with band 2.'
      print '* "timeres crosscor21" to plot the time-resolved cross-correlation function of band 1 with band 2'
      print '* "all" to plot all available data products.'
   if nfiles==3:                                                           # Only display 3-data-set instructions if 3 datasets given
      print ''
      print '3 DATASET PLOTS:'
      print '* "hardness32" to plot a hardness/time diagram of file3/file2 colour over time.'
      print '* "hardness23" to plot a hardness/time diagram of file2/file3 colour over time.'
      print '* "hardness31" to plot a hardness/time diagram of file3/file1 colour over time.'
      print '* "hardness13" to plot a hardness/time diagram of file1/file3 colour over time.'
      print '* "hid32" to plot a hardness-intensity diagram of file3/file2 colour against total flux.'
      print '* "hid23" to plot a hardness-intensity diagram of file2/file3 colour against total flux.'
      print '* "hid31" to plot a hardness-intensity diagram of file3/file1 colour against total flux.'
      print '* "hid13" to plot a hardness-intensity diagram of file1/file3 colour against total flux.'
      print '* "col32" to plot file3/file2 colour against time.'
      print '* "col23" to plot file2/file3 colour against time.'
      print '* "col31" to plot file3/file1 colour against time.'
      print '* "col13" to plot file1/file3 colour against time.'
      print '* "compbands31" to plot lightcurves of bands 3 and 1 against each other.'
      print '* "compbands32" to plot lightcurves of bands 3 and 2 against each other.'
      print '* "ccd" to plot a colour-colour diagram (3/1 colour against 2/1 colour).'
      print '* "timeres crosscor31" to plot the time-resolved cross-correlation function of band 3 with band 1'
      print '* "timeres crosscor32" to plot the time-resolved cross-correlation function of band 3 with band 2'
      print '* "crosscor31" to plot the cross-correlation function of band 3 with band 1.'
      print '* "crosscor32" to plot the cross-correlation function of band 3 with band 2.'
   print ''
   print 'BURST ANALYSIS:'
   print '* "burst get" to interactively extract burst data for analysis.'
   print '* "burst peaks" for a histogram of peak heights of extracted bursts.'
   print '* "burst risetimes" for a histogram of rise times of extracted bursts.'
   print '* "burst falltimes" for a histogram of fall times of extracted bursts.'
   print '* "burst lengths" for a histogram of durations of extracted bursts.'
   print '* "burst help" for further information on burst analysis.'
   print ''
   print 'SAVING DATA TO ASCII:'
   print '* "export" to dump the lightcurve and colour data into an ASCII file.'
   print '* "bgdump" to export background lightcurve to an ASCII file.'
   print '* "timenorm" to toggle absolute or relative time values on x-axis.'
   print ''
   print 'TOGGLE OPTIONS:'
   print '* "errors" to toggle whether to display errors in plots.'
   print '* "lines" to toggle lines joining points in graphs.'
   print '* "ckey" to toggle colour key (red-blue) for the first five points in all plots.'
   print '* "save" to save to disk any plots which would otherwise be shown.'
   print ''
   print 'ADVANCED OPTIONS:'
   print '* "burstalg" to select algorithm for finding pulse peaks in lightcurve.'
   print ''
   print 'OTHER COMMANDS:'
   print '* "info" to display a list of facts and figures about the current PlotDemon session.'
   print '* "reflav" to rewrite the flavour text used for graph titles.'
   print '* "help" or "?" to display this list of instructions again.'
   print '* "quit" to quit.'

#give_inst()                                                              # Print the list of instructions
print ''
print ' --------------------'


#-----Entering Interactive Mode----------------------------------------------------------------------------------------

while plotopt not in ['quit','exit']:                                     # If the previous command given was not quit, continue

   print ''
   plotopt=raw_input('Give command [? for help]: ').lower()               # Fetch command from user
   print ''


   #-----Aliasing options----------------------------------------------------------------------------------------------

   if plotopt=='shid':                                                    # Let 'shid' refer to the soft 2/1 HID
      plotopt='hid21'

   elif plotopt=='hhid':                                                  # Let 'hhid' refer to the hard 3/1 HID
      plotopt='hid31'

   elif plotopt=='hid' and nfiles==2:
      plotopt='hid21'                                                     # Let 'hid' refer to the 2/1 HID if that is the only HID available


   #-----Hidden 'stick' option-----------------------------------------------------------------------------------------

   if plotopt=='stick':                                                   # For use when scripting with Plotdemon.  If turned on, this causes
                                                                          #  all plots to block when shown.
      show_block=not show_block
      if show_block:
         print 'Sticky Plots on!'
      else:
         print 'Sticky Plots off!'


   #-----'save' option-------------------------------------------------------------------------------------------------

   elif plotopt=='save':                                                  # Causes a plot to be saved when it would otherwise have been shown
      saveplots=not saveplots
      if saveplots:
         print 'Plot saving on!'
      else:
         print 'Plot saving off!'


   #-----'rebin' option------------------------------------------------------------------------------------------------

   elif plotopt=='rebin':                                                 # Rebin data

      bursts=None                                                         # Remove burst data
      fldtxt=''
      while True:                                                         # Keep asking until a good response is given
         try:
            binning=float(raw_input("Enter bin size (s): "))              # Ask for binsize in dialogue box
            assert binning>=minbin
            break
         except:
            print 'Invalid bin size input!'

      print 'Binning File 1...'
      x1,y1,ye1=pan.binify(x1r,y1r,ye1r,binning)                          # Bin File 1 using 'binify' in pan_lib
      if nfiles>1:
         print 'Binning File 2...'
         x2,y2,ye2=pan.binify(x2r,y2r,ye2r,binning)                       # Bin File 2 using 'binify' in pan_lib
         if nfiles>2:
            print 'Binning File 3...'
            x3,y3,ye3=pan.binify(x3r,y3r,ye3r,binning)                    # Bin File 3 using 'binify' in pan_lib

      if nfiles>1:                                                        # Checking file lengths match... need to migrate this to a function one day
         if len(x1)!=len(x2):
            wrongsize=True

      if nfiles==3:
         if len(x1)!=len(x3):
            wrongsize=True
            x3l=len(x3)

      if wrongsize:                                                       # Forcing file lengths to match if possible
         mindex=min(len(x1),len(x2),x3l)-1
         if x1[mindex]!=x2[mindex]:
            print 'Cannot crop, aborting!'
            pan.signoff()
            exit()
         if nfiles==3:
            if x1[mindex]!=x3[mindex]:
               print 'Cannot crop, aborting!'
               pan.signoff()
               exit()
         mindex+=1
         x1=x1[:mindex]
         y1=y1[:mindex]
         ye1=ye1[:mindex]
         x2=x2[:mindex]
         y2=y2[:mindex]
         ye2=ye2[:mindex]
         if nfiles==3:
            x3=x3[:mindex]
            y3=y3[:mindex]
            ye3=ye3[:mindex]

      gmask=pan.gtimask(x1,gti)                                           # Re-establish gmask

      print 'Binning complete!'
      print ''

      times,timese,flux,fluxe,ys,yes,col,cole=colorget()                  # Re-get colours
      folded=False                                                        # Re-allow clipping
      print 'Done!'
      print ''


   #-----'fold' Option-------------------------------------------------------------------------------------------------

   elif plotopt=='fold':                                                  # Fold lightcurve

      bursts=None                                                         # Remove burst data
      if folded:
         print 'Data already folded!  Rebin before re-folding.'
         continue

      if not module_pyastro:                                              # Only attempt to fold if pyastro is present
         print 'PyAstronomy Module not found!  Cannot perform fold!'      # Warn user they cannot fold as module is missing
         continue

      while True:                                                         # Keep asking user until they give a sensible period
         try:
            period=float(raw_input('Input period to fold over (s): '))    # Fetch period from user
            break
         except:
            print "Invalid period!"                                       # Keep trying until they give a sensible input

      while True:                                                         # Keep asking user until they give a sensible phase resolution
         try:
            phres=float(raw_input('Input phase resolution (0-1): '))      # Fetch phase resolution from user
            assert phres<1.0
            assert phres>0.0
            break
         except:
            print "Invalid phase resolution!"                             # Keep trying until they give a sensible input

      x1=x1[gmask];y1=y1[gmask];ye1=ye1[gmask]                            # Zeroing all data points outside of GTI
      x1,y1,ye1=pan.foldify(x1,y1,ye1,period,binning,phres=phres,name='ch. '+ch[1])    # Fold using foldify function from pan_lib

      fldtxt='Folded '

      if nfiles>1:
         x2=x2[gmask];y2=y2[gmask];ye2=ye2[gmask]                         # Zeroing all data points outside of GTI
         x2,y2,ye2=pan.foldify(x2,y2,ye2,period,binning,phres=phres,name='ch. '+ch[2]) # Fold data of file 2 if present

      if nfiles==3:
         x3=x3[gmask];y3=y3[gmask];ye3=ye3[gmask]                         # Zeroing all data points outside of GTI
         x3,y3,ye3=pan.foldify(x3,y3,ye3,period,binning,phres=phres,name='ch. '+ch[3]) # Fold data of file 3 if present

      gmask=np.ones(len(x1),dtype=bool)                                   # Re-establish gmask
      times,timese,flux,fluxe,ys,yes,col,cole=colorget()                  # Re-get colours
      folded=True

      print 'Folding Complete!'
      print ''


   #-----'autofold' Option-------------------------------------------------------------------------------------------------

   elif plotopt=='autofold':                                              # Autofold data lightcurve

      bursts=None                                                         # Remove burst data
      if folded:
         print 'Data already folded!  Rebin before re-folding.'
         continue

      if not module_pyastro:                                              # Only attempt to fold if pyastro is present
         print 'PyAstronomy Module not found!  Cannot perform fold!'      # Warn user they cannot fold as module is missing
         continue

      ls_st=max(4.0/(times[-1]-times[0]),0.005)
      ls_end=0.5/binning

      lsx=np.arange(ls_st,ls_end,(ls_end-ls_st)/2500.0)                   # Perform Lomb-Scargle Analysis on the data to seek best period
      lsy=pan.lomb_scargle(times,flux,fluxe,lsx)
      period=1.0/(lsx[lsy.tolist().index(max(lsy))])

      while True:                                                         # Keep asking user until they give a sensible phase resolution
         try:
            phres=float(raw_input('Input phase resolution (0-1): '))      # Fetch phase resolution from user
            assert phres<1.0
            assert phres>0.0
            break
         except:
            print "Invalid phase resolution!"                             # Keep trying until they give a sensible input

      print ''
      print 'Using period of '+str(period)+'!'

      x1=x1[gmask];y1=y1[gmask];ye1=ye1[gmask]                            # Zeroing all data points outside of GTI
      x1,y1,ye1=pan.foldify(x1,y1,ye1,period,binning,phres=phres,name='ch. '+ch[1])    # Fold using foldify function from pan_lib

      fldtxt='Folded '

      if nfiles>1:
         x2=x2[gmask];y2=y2[gmask];ye2=ye2[gmask]                         # Zeroing all data points outside of GTI
         x2,y2,ye2=pan.foldify(x2,y2,ye2,period,binning,phres=phres,name='ch. '+ch[2]) # Fold data of file 2 if present

      if nfiles==3:
         x3=x3[gmask];y3=y3[gmask];ye3=ye3[gmask]                         # Zeroing all data points outside of GTI
         x3,y3,ye3=pan.foldify(x3,y3,ye3,period,binning,phres=phres,name='ch. '+ch[3]) # Fold data of file 3 if present

      gmask=np.ones(len(x1),dtype=bool)                                   # Re-establish gmask
      times,timese,flux,fluxe,ys,yes,col,cole=colorget()                  # Re-get colours
      folded=True

      print 'Folding Complete!'
      print ''


   #-----'varifold' Option---------------------------------------------------------------------------------------------

   elif plotopt=='varifold':

      if folded:
         print 'Cannot perform burst analysis on folded data!'
         continue

      if burst_alg=='cubic spline':
         while True:
            try:
               iq_lo=float(raw_input('Low Threshold:  '))
               iq_hi=float(raw_input('High Threshold: '))
               assert iq_hi>iq_lo
               assert iq_hi<=100
               assert iq_lo>=0
               break
            except AssertionError:
               print 'Invalid Entry!  Valid entry is of the form High>Low.'
      else:
         iq_lo=0
         iq_hi=100

      while True:
         try:
            phase_res=float(raw_input('Input phase resolution (0-1): '))
            assert phase_res<1.0
            assert phase_res>0.0
            break
         except AssertionError:
            print 'Invalid Phase Resolution!'

      phases=pan.fold_bursts(times,flux,iq_hi,iq_lo,do_smooth=True,alg=burst_alg)

      ymask=(phases!=np.inf)

      nbins=int(1.0/phase_res)
      phases=(nbins*phases[ymask]).astype(int)

      x1=x1[gmask][ymask];y1=y1[gmask][ymask];ye1=ye1[gmask][ymask]       # Removing all data points outside of GTI

      newx1=[]
      newy1=[]
      newye1=[]

      for i in range(nbins):
         newx1.append(float(i)/float(nbins))
         newy1.append(np.mean(y1[phases==i]))
         newye1.append((np.sum(ye1[phases==i]**2))**0.5/len(ye1[phases==i]))

      x1=np.array(newx1)
      #y1=sig.savgol(np.array(newy1),5,3)
      #ye1=sig.savgol(np.array(newye1),5,3)
      y1=np.array(newy1)
      ye1=np.array(newye1)

      if nfiles>1:

         x2=x2[gmask][ymask];y2=y2[gmask][ymask];ye2=ye2[gmask][ymask]    # Removing all data points outside of GTI

         newx2=[]
         newy2=[]
         newye2=[]

         for i in range(nbins):
            newx2.append(float(i)/float(nbins))
            newy2.append(np.mean(y2[phases==i]))
            newye2.append((np.sum(ye2[phases==i]**2))**0.5/len(ye2[phases==i]))

         x2=np.array(newx2)
         y2=np.array(newy2)
         ye2=np.array(newye2)

      if nfiles==3:
 
         x3=x3[gmask][ymask];y3=y3[gmask][ymask];ye3=ye3[gmask][ymask]    # Removing all data points outside of GTI

         newx3=[]
         newy3=[]
         newye3=[]

         for i in range(nbins):
            newx3.append(float(i)/float(nbins))
            newy3.append(np.mean(y3[phases==i]))
            newye3.append((np.sum(ye3[phases==i]**2))**0.5/len(ye3[phases==i]))

         x3=np.array(newx3)
         y3=np.array(newy3)
         ye3=np.array(newye3)

      gmask=np.ones(len(x1),dtype=bool)                                   # Re-establish gmask
      times,timese,flux,fluxe,ys,yes,col,cole=colorget()                  # Re-get colours

      folded=True

      print 'Folding Complete!'
      print ''
      period='N/A'


   #-----'clip' Option-------------------------------------------------------------------------------------------------

   elif plotopt=='clip':                                                  # Clipping data

      bursts=None                                                         # Remove burst data
      if folded:
         print 'Cannot clip folded data!'

      else:

         print 'Clipping data'
         print ''
         print 'Time range is '+str(x1[0])+'s - '+str(x1[-1])+'s'

         print 'Please choose new range of data:'
         mint,maxt,srbool=pan.srinr(x1,binning,'time')                    # Fetch new time domain endpoints using srinr function from pan_lib

         if srbool:

            print 'Clipping...'

            x1=x1[mint:maxt]                                              # Clip file 1
            y1=y1[mint:maxt]
            ye1=ye1[mint:maxt]

            if nfiles>1:
               x2=x2[mint:maxt]                                           # Clip file 2
               y2=y2[mint:maxt]
               ye2=ye2[mint:maxt]

            if nfiles==3:
               x3=x3[mint:maxt]                                           # Clip file 3
               y3=y3[mint:maxt]
               ye3=ye3[mint:maxt]

            gmask=pan.gtimask(x1,gti)                                     # Re-establish gmask
            times,timese,flux,fluxe,ys,yes,col,cole=colorget()            # Re-get colours

            print 'Data clipped!'


   #-----'mask' Option-------------------------------------------------------------------------------------------------

   elif plotopt=='mask':

      bursts=None                                                         # Remove burst data
      if folded:
         print 'Cannot mask folded data!'

      else:

         print 'Masking data'
         print ''
         print 'Select time range to mask: '
         mint,maxt,srbool=pan.srinr(x1,binning,'time')                    # Fetch time domain endpoints of bad window using srinr function from pan_lib

         if srbool:
            print 'Masking...'

            gmask[mint:maxt]=False                                        # Force all values inside the bad window to appear as outside of GTIs

            times,timese,flux,fluxe,ys,yes,col,cole=colorget()            # Re-get colours

            print 'Data masked!'


   #-----'rms' Option--------------------------------------------------------------------------------------------------

   elif plotopt=='rms':

      if nfiles>1:                                                        # If more than one file loaded, prompt user to select one
         if nfiles==3:
            is_band_3=', 3'
         else:
            is_band_3=''
         selected_band=raw_input('Select Energy Band [1, 2'+is_band_3+', All]: ').lower()
         if selected_band not in ['1','2','3','all']:
            print 'Invalid band!'
            continue
      else:
         selected_band='all'

      if selected_band=='1':                                              # Fetch the rms
         rms=pan.rms(y1)
      elif selected_band=='2':
         rms=pan.rms(y2)
      elif selected_band=='3':
         rms=pan.rms(y3)
      else:
         rms=pan.rms(flux)

      if rms=='div0':                                                     # If the mean of the data is 0, fail safely
         print 'Error!  Div 0 Encountered!  Aborting!'
      else:
         print 'rms =',str(rms*100)+'%'                                   # Otherwise, print RMS


   #-----'lc' Option---------------------------------------------------------------------------------------------------

   elif plotopt=='lc':                                                    # Plot lightcurve

      taxis='Phase' if folded else 'Time (s)'

      pl.figure()
      doplot(times,timese,flux,fluxe,ovr=True,per2=folded)                # Plot flux/time using doplot from pan_lib
      pl.xlabel(taxis)
      pl.ylabel(flux_axis)
      pl.ylim(ymin=0)
      pl.title(fldtxt+'Lightcurve'+qflav)
      plot_save(saveplots,show_block)
      print 'Lightcurve plotted!'


   #-----'bg' Option---------------------------------------------------------------------------------------------------

   elif plotopt=='bg':                                                    # Plot background lightcurve

      is_bdata=(bdata1!=None)
      if nfiles>=2:
         is_bdata=is_bdata and (bdata2!=None)
      if nfiles==3:
         is_bdata=is_bdata and (bdata3!=None)

      if is_bdata:                                                        # Only proceed if all infiles have background data

         try:
            bdata_x=bdata1[0]
            bdata_y=bdata1[1]
            if nfiles>=2:
               bdata_y+=bdata2[1]
            if nfiles==3:
               bdata_y+=bdata3[1]                                         # Sum counts of all <=3 backgrounds
            
            pl.figure()
            pl.plot(bdata_x,bdata_y)
            pl.xlabel('Time (s)')
            pl.ylabel(flux_axis)
            pl.ylim(ymin=0)
            pl.title(fldtxt+'Lightcurve'+qflav)
            plot_save(saveplots,show_block)
            print 'Background plotted!'
         except:
            print 'Backgrounds inconsistenly formatted!'                  # Abort if backgrounds are somehow of different lengths'

      else:
         print 'Not all files have background data!'


   #-----'bgdump' Option-----------------------------------------------------------------------------------------------

   elif plotopt=='bgdump':                                                # Export background lightcurve

      is_bdata=(bdata1!=None)
      if nfiles>=2:
         is_bdata=is_bdata and (bdata2!=None)
      if nfiles==3:
         is_bdata=is_bdata and (bdata3!=None)

      if is_bdata:                                                        # Only proceed if all infiles have background data

         try:
            bdata_x=bdata1[0]
            bdata_y=bdata1[1]
            if nfiles>=2:
               bdata_y+=bdata2[1]
            if nfiles==3:
               bdata_y+=bdata3[1]                                         # Sum counts of all <=3 backgrounds

         except:
            print 'Backgrounds inconsistenly formatted!'                  # Abort if backgrounds are somehow of different lengths'
            continue

         ofilename=raw_input('Save textfile as: ')                     # Fetch filename from user

         ofil = open(ofilename, 'w')                                   # Open file

         for i in range(len(bdata_x)):

            row=['0.0']*3                                              # Create a row of strings reading 0.0, append data into it
            row[0]=str(bdata_x[i]+time_n)+' '                          # Column 01: Time
            row[1]=str(bdata_y[i])+' '                                 # Column 02: Rate
            row[2]='\n'
            ofil.writelines(row) 

         print 'Background saved as '+ofilename+'!'


      else:
         print 'Not all files have background data!'


   #-----'export' Option-----------------------------------------------------------------------------------------------

   elif plotopt=='export':                                                # Export lightcurve to ASCII file

      ofilename=raw_input('Save textfile as: ')                           # Fetch filename from user

      ofil = open(ofilename, 'w')                                         # Open file

      for i in range(len(times)):

         row=['0.0 ']*15                                                  # Create a row of strings reading 0.0, append data into it
         row[0]=str(times[i]+time_n)+' '                                         # Column 00: Time
         row[1]=str(timese[i])+' '                                        # Column 01: Time Error
         row[2]=str(flux[i])+' '                                          # Column 02: Total Flux
         row[3]=str(fluxe[i])+' '                                         # Column 03: Total Flux Error
         row[4]=str(y1[gmask][i])+' '                                     # Column 04: Band 1 Flux
         row[5]=str(ye1[gmask][i])+' '                                    # Column 05: Band 1 Flux Error
         row[14]='\n'                                                     # Column 14: Return (so further data will be appended to a new line)

         if nfiles>1:                                                     # If 2+ bands are given:

            row[6]=str(y2[gmask][i])+' '                                  # Column 06: Band 2 Flux
            row[7]=str(ye2[gmask][i])+' '                                 # Column 07: Band 2 Flux Error
            row[10]=str(col[21][i])+' '                                   # Column 10: [2/1] Colour
            row[11]=str(cole[21][i])+' '                                  # Column 11: [2/1] Colour Error

         if nfiles==3:                                                    # If 3 bands are given:

            row[8]=str(y3[gmask][i])+' '                                  # Column 08: Band 3 Flux
            row[9]=str(ye3[gmask][i])+' '                                 # Column 09: Band 3 Flux Error
            row[12]=str(col[31][i])+' '                                   # Column 12: [3/1] Colour
            row[13]=str(cole[31][i])+' '                                  # Column 13: [3/1] Colour Error

         ofil.writelines(row)                                             # Append row of data into open file

      ofil.close()                                                        # Close file

      print 'Data saved!'


   #-----'timenorm' Option---------------------------------------------------------------------------------------------
   
   elif plotopt=='timenorm':

      print ''
      if time_n==0:
         time_n=tst1
         print 'Using absolute time!'
      else:
         time_n=0
         print 'Using relative time!'
      print ''

   #-----'animate' Option----------------------------------------------------------------------------------------------

   elif plotopt=='animate':

      animsloc=raw_input('Folder to save images: ')
      print ''

      if os.path.exists(animsloc):                                        # Create the folder
         print 'Folder "'+animsloc+'" already exists...'
      else:
         print 'Creating folder "'+animsloc+'"...'
         os.makedirs(animsloc)

      here=os.getcwd()                                                    # Get current working directory (to move back to later)
      os.chdir(animsloc)                                                  # Change working directory to animation location

      animbin=0.0025                                                      # Start with an arbitrarily low binsize

      while animbin<max(bsz1,bsz2,bsz3,minbin):                           # Find lowest allowable binsize of the form 0.01*2^N
         animbin=animbin*2

      anstep=1                                                            # Track the number of steps taken
      dst=times[0]                                                        # Fetch largest and smallest times to use to force same scale on all graphs
      det=times[-1]

      while animbin<(det-dst)/4.0:                                        # Set the maximum binsize at one quarter of the observation length

         print "Creating",str(animbin)+"s binned lightcurve"

         x1,y1,ye1=pan.binify(x1r,y1r,ye1r,animbin)                       # Bin File 1 using 'binify' in pan_lib
         if nfiles>1:
            x2,y2,ye2=pan.binify(x2r,y2r,ye2r,animbin)                    # Bin File 2 using 'binify' in pan_lib
            if nfiles>2:
               x3,y3,ye3=pan.binify(x3r,y3r,ye3r,animbin)                 # Bin File 3 using 'binify' in pan_lib

         mina,maxa,srbool=pan.srinr(x1,binning,'time',minv=dst,maxv=det)  # Clip each individual lightcurve

         if srbool:

            x1=x1[mina:maxa]                                              # Clip file 1
            y1=y1[mina:maxa]
            ye1=ye1[mina:maxa]

            if nfiles>1:

               x2=x2[mina:maxa]                                           # Clip file 2
               y2=y2[mina:maxa]
               ye2=ye2[mina:maxa]

            if nfiles==3:

               x3=x3[mint:maxt]                                           # Clip file 3
               y3=y3[mint:maxt]
               ye3=ye3[mint:maxt]

         gmask=pan.gtimask(x1,gti)                                        # Re-establish gmask

         times,timese,flux,fluxe,ys,yes,col,cole=colorget(verbose=False)  # Re-get colours

         if anstep==1:
            if es:
              merr=max(fluxe)
            else:
              merr=0
            maxany=max(flux)+merr                                         # Calculate the scale of all plots based on the range of the first plot
            minany=min(flux)-merr
            if minany<0: minany=0

         taxis='Phase' if folded else 'Time (s)'

         pl.figure()
         doplot(times,timese,flux,fluxe,ovr=True)                         # Plot the graph using doplot from pan_lib
         pl.xlabel(taxis)
         pl.ylabel(flux_axis)
         pl.title('Lightcurve ('+str(animbin)+'s binning)')
         pl.xlim(dst,det)
         pl.ylim(max(minany,0),maxany)
         pl.savefig(str("%04d" % anstep)+'.png')                          # Save the figure with leading zeroes to preserve order when int convereted to string
         pl.close()

         anstep+=1                                                        # Increment the step tracker
         animbin=animbin*2                                                # Double the binsize

      print 'Cleaning up...'

      os.system ("convert -delay 10 -loop 0 *.png animation.gif")         # Use the bash command 'convert' to create the animated gif

      x1,y1,ye1=pan.binify(x1r,y1r,ye1r,binning)                          # Reset binning of File 1 using 'binify' in pan_lib
      if nfiles>1:
         x2,y2,ye2=pan.binify(x2r,y2r,ye2r,binning)                       # Reset binning of File 2 using 'binify' in pan_lib
         if nfiles>2:
            x3,y3,ye3=pan.binify(x3r,y3r,ye3r,binning)                    # Reset binning File 3 using 'binify' in pan_lib

      gmask=pan.gtimask(x1,gti)                                           # Re-establish gmask
      times,timese,flux,fluxe,ys,yes,col,cole=colorget(verbose=False)     # Re-get colours

      print ''
      print "Animation saved to",animsloc+'/animation.gif!'
      os.chdir(here)        


   #-----'circfold' Option----------------------------------------------------------------------------------------------

   elif plotopt=='circfold':

      goodfp=False

      while True:                                                         # Keep trying until user gives a sensible fold-time

         try:
            cftime=float(raw_input('Fold Period (s): '))
            break
         except:
            print 'Invalid entry!'

      theta,rad=pan.circfold(times,flux,cftime)
      limany=max(rad)                                                     # Fetch radius of plot

      pl.figure()
      ax = pl.subplot(111, polar=True)

      # The binned curve #

      bintheta=[]
      binrad=[]
      avetheta=[]
      averad=[]

      for i in range(0,int(cbin)):                                        # For each angular segment...
         mask=((theta>=(2*pi*i/cbin)) & (theta<(2*pi*(i+1)/cbin)))        # Filter all points outside the segment
         if sum(mask)>0:                                                  # Check that points actually fall within the angular segment
            binrad.append(np.mean(rad[mask]))                             # Average all points within that phase segment
            bintheta.append(2*pi*((2*i+1)/2.0)/cbin)
            averad.append(np.mean(rad))
            avetheta.append(2*pi*((2*i+1)/2.0)/cbin)

      inhom=np.sqrt(np.mean(((np.array(binrad)-np.mean(rad))/np.mean(rad))**2)) # Inhomogeneity is a measure of standard deviation from the average circle

      print 'Inhomogeneity =',inhom

      bintheta.append(bintheta[0])
      binrad.append(binrad[0])
      avetheta.append(avetheta[0])
      averad.append(averad[0])

      ax.plot(theta,rad,color='#999999',linestyle='', marker='x')         # Print datapoints
      ax.plot(avetheta,averad,'k-')                                       # Print average
      ax.plot(bintheta,binrad,'b-')                                       # Print inhomogeneity line
      ax.set_title('Circle Diagram ('+str(cftime)+'s folding)')
      ax.set_rmax(1.05*limany)
      plot_save(saveplots,show_block)
      pl.close()


   #-----'circanim' Option----------------------------------------------------------------------------------------------

   elif plotopt=='circanim':

      cistep=1                                                            # Track the number of steps taken
      dst=times[0]                                                        # Fetch largest and smallest fluxes to use to force same scale on all graphs
      det=times[-1]

      circsloc=raw_input('Folder to save images: ')
      try:
         castt=float(raw_input('Lowest Fold time  (s): '))
      except:
         castt=binning
      try:
         caend=float(raw_input('Highest Fold time (s): '))
      except:
         caend=(det-dst)/4.0
      print ''

      if os.path.exists(circsloc):                                        # Create the folder
         print 'Folder "'+circsloc+'" already exists...'
      else:
         print 'Creating folder "'+circsloc+'"...'
         os.makedirs(circsloc)

      here=os.getcwd()                                                    # Get current working directory (to move back to later)
      os.chdir(circsloc)                                                  # Change working directory to animation location

      caend=min((det-dst)/4.0,caend)

      foldinx=castt

      inhlist=[]
      maxinh=0                                                            # Set max inhomogeneity to zero initially
      maxinl=0

      while foldinx<caend:                                                # Set the maximum old index at one quarter of the observation length

         print "Creating",str(foldinx)+"s folded Circle Diagram"

         theta,rad=pan.circfold(times,flux,foldinx)

         if cistep==1:
            limany=max(rad)

         pl.figure()
         ax = pl.subplot(111, polar=True)

         # The binned curve #

         bintheta=[]
         binrad=[]
         avetheta=[]
         averad=[]

         for i in range(0,int(cbin)):                                     # For each angular segment...
            mask=((theta>=(2*pi*i/cbin)) & (theta<(2*pi*(i+1)/cbin)))     # Filter all points outside the segment
            if sum(mask)>0:                                               # Check that points actually fall within the angular segment
               binrad.append(np.mean(rad[mask]))                          # Average all points within that phase segment
               bintheta.append(2*pi*((2*i+1)/2.0)/cbin)
               averad.append(np.mean(rad))
               avetheta.append(2*pi*((2*i+1)/2.0)/cbin)

         inhom=np.sqrt(np.mean(((np.array(binrad)-np.mean(rad))/np.mean(rad))**2)) # Inhomogeneity is a measure of standard deviation from the average circle
         inhlist.append(inhom)

         if max(inhlist)>maxinh:
            maxinh=max(inhlist)
            maxinl=foldinx

         bintheta.append(bintheta[0])
         binrad.append(binrad[0])
         avetheta.append(avetheta[0])
         averad.append(averad[0])        


         ax.plot(theta,rad,color='#999999',linestyle='', marker='x')      # Print datapoints
         ax.plot(avetheta,averad,'k-')                                    # Print average
         ax.plot(bintheta,binrad,'b-')                                    # Print inhomogeneity line
         ax.set_title('Circle Diagram ('+str(foldinx)+'s folding)')
         ax.set_rmax(1.05*limany)
         ins=pl.axes([0.79,0.74,0.2,0.2])                                 # Save the inset inhomogeneity tracker
         pl.setp(ins,xticks=[],yticks=[0.05,0.1,0.15,0.2],xlim=[0,100],ylim=[0,max(inhlist)],title='Inhom.')
         ins.plot(inhlist)
         pl.savefig(str("%04d" % cistep)+'.png')                          # Save the figure with leading zeroes to preserve order when int convereted to string
         pl.close()

         cistep+=1
         foldinx+=((caend-castt)/100.0)

      print 'Cleaning up...'

      os.system ("convert -delay 10 -loop 0 *.png animation.gif")         # Use the bash command 'convert' to create the animated gif

      print ''
      print 'Maximum inhomogeneity =',maxinh,'at',str(maxinl)+'s'
      print "Animation saved to",circsloc+'/animation.gif!'
      os.chdir(here)        


   #-----'hidxy' Option------------------------------------------------------------------------------------------------

   elif plotopt[:3]=='hid':                                                # Plot x/y HID

      ht=plotopt[3:]                                                       # Collect the xy token from the user

      if nfiles>1:
         if not (ht in ['12','13','21','23','31','32']):                   # Check that the token is 2 long and contains two different characters of the set [1,2,3]

            print 'Invalid command!'
            print ''
            print 'Did you mean...'
            print ''
            print 'HID options:'
            print '* "hid21" for 2/1 colour'
            print '* "hid12" for 1/2 colour'
            if nfiles==3:
               print '* "hid32" for 3/2 colour'
               print '* "hid23" for 2/3 colour'
               print '* "hid31" for 3/1 colour'
               print '* "hid13" for 1/3 colour'

         elif ('3' in ht) and (nfiles<3):

            print 'Not enough infiles for advanced HID!'                  # If token contains a 3 but only 2 infiles are used, abort!

         else:

            h1=int(ht[0])                                                 # Extract numerator file number
            h2=int(ht[1])                                                 # Extract denominator file number
            ht=int(ht)
            pl.figure()
            doplot(col[ht],cole[ht],flux,fluxe)                           # Collect colours from col library and plot
            pl.ylabel(flux_axis)
            pl.xlabel('('+ch[h1]+'/'+ch[h2]+') colour')
            pl.title(fldtxt+'Hardness Intensity Diagram'+qflav)
            plot_save(saveplots,show_block)
            print 'File'+str(h1)+'/File'+str(h2)+' HID plotted!'

      else:
         print 'Not enough infiles for HID!'


   #-----'hardnessxy' Option------------------------------------------------------------------------------------------------

   elif plotopt[:8]=='hardness':                                           # Plot x/y hardness/time plot

      ht=plotopt[8:]                                                       # Collect the xy token from the user

      if nfiles>1:
         if not (ht in ['12','13','21','23','31','32']):                   # Check that the token is 2 long and contains two different characters of the set [1,2,3]

            print 'Invalid command!'
            print ''
            print 'Did you mean...'
            print ''
            print 'HID options:'
            print '* "hardness21" for 2/1 hardness over time'
            print '* "hardness12" for 1/2 hardness over time'
            if nfiles==3:
               print '* "hardness32" for 3/2 hardness over time'
               print '* "hardness23" for 2/3 hardness over time'
               print '* "hardness31" for 3/1 hardness over time'
               print '* "hardness13" for 1/3 hardness over time'

         elif ('3' in ht) and (nfiles<3):

            print 'Not enough infiles for advanced HID!'                  # If token contains a 3 but only 2 infiles are used, abort!

         else:

            h1=int(ht[0])                                                 # Extract numerator file number
            h2=int(ht[1])                                                 # Extract denominator file number
            ht=int(ht)
            pl.figure()
            doplot(times,timese,col[ht],cole[ht],ovr=True,per2=folded)    # Collect colours from col library and plot
            pl.ylabel('('+ch[h1]+'/'+ch[h2]+') colour')
            pl.xlabel(taxis)
            pl.title(fldtxt+'Hardness/Time plot'+qflav)
            plot_save(saveplots,show_block)
            print 'File'+str(h1)+'/File'+str(h2)+' Hardness/time diagram plotted!'

      else:
         print 'Not enough infiles for HID!'



   #-----'compbandsxy' Option------------------------------------------------------------------------------------------

   elif plotopt[:9]=='compbands':                                          # Plot two bands against each other

      ht=plotopt[9:]                                                       # Collect the xy token from the user
      if nfiles>1:
         if not (ht in ['12','13','21','23','31','32']):                   # Check that the token is 2 long and contains two different characters of the set [1,2,3]

            print 'Invalid command!'
            print ''
            print 'Did you mean...'
            print ''
            print 'HID options:'
            print '* "compbands21" to plot lightcurve 2 against 1'
            print '* "compbands12" to plot lightcurve 1 against 2'
            if nfiles==3:
               print '* "compbands32" to plot lightcurve 3 against 2'
               print '* "compbands23" to plot lightcurve 2 against 3'
               print '* "compbands31" to plot lightcurve 3 against 1'
               print '* "compbands13" to plot lightcurve 1 against 3'

         elif ('3' in ht) and (nfiles<3):

            print 'Not enough infiles for advanced lightcurve comparison!'# If token contains a 3 but only 2 infiles are used, abort!

         else:

            h1=int(ht[0])                                                 # Extract numerator file number
            h2=int(ht[1])                                                 # Extract denominator file number
            pl.figure()
            doplot(ys[h1],yes[h1],ys[h2],yes[h2])                         # Collect colours from col library and plot
            pl.ylabel(r'Band '+str(h1)+r' rate (cts s$^{-1}$ PCU$^{-1}$)')
            pl.xlabel(r'Band '+str(h2)+r' rate (cts s$^{-1}$ PCU$^{-1}$)')
            pl.title(fldtxt+'Lightcurve Comparison Diagram'+qflav)
            plot_save(saveplots,show_block)
            print 'File'+str(h1)+'/File'+str(h2)+' LC Comparison plotted!'

      else:
         print 'Not enough infiles for lightcurve comparison!'


   #-----'autocor' Option------------------------------------------------------------------------------------------------

   elif plotopt=='autocor':                                                # Plot autocorrelation function

      ccor=sig.correlate(flux-np.mean(flux),flux-np.mean(flux),mode='same')/(len(flux)*(np.std(flux)**2))
         
      nlen=len(times)

      cct=range(nlen)                                                # Set up the lag axis
      cct=((np.array(cct)-nlen/2.0)*binning).tolist()

      pl.figure()
      pl.plot(cct,ccor)
      pl.xlabel('Lag (s)')
      pl.ylabel('Cross-Correlation')
      plot_save(saveplots,show_block)
      print 'Autocorrelation diagram plotted!'


   #-----'crosscor' Option------------------------------------------------------------------------------------------------

   elif plotopt[:8]=='crosscor':                                           # Plot cross-correlation function

      ht=plotopt[8:]                                                       # Collect the xy token from the user

      if nfiles>1:
         if not (ht in ['12','13','21','23','31','32']):                   # Check that the token is 2 long and contains two different characters of the set [1,2,3]

            print 'Invalid command!'
            print ''
            print 'Did you mean...'
            print ''
            print 'HID options:'
            print '* "crosscor21" for 2 against 1 cross-correlation'
            if nfiles==3:
               print '* "crosscor32" for 3 against 2 cross-correlation'
               print '* "crosscor31" for 3 against 1 cross-correlation'

         elif ('3' in ht) and (nfiles<3):

            print 'Not enough infiles for advanced cross-correlation!'    # If token contains a 3 but only 2 infiles are used, abort!

         else:

            nlen=len(times)
            h1=int(ht[0])                                                 # Extract file 1 number
            h2=int(ht[1])                                                 # Extract file 2 number
            h1a=ys[h1]                                                    # Extract dataset 1
            h2a=ys[h2]                                                    # Extract dataset 2

            if nlen%2==0:                                                 # Cross-correlation only gives a point at zero lag if an even number of points are input
               nlen-=1
               h1a=h1a[:-1]
               h2a=h2a[:-1]

            ccor=sig.correlate(h1a-np.mean(h1a),h2a-np.mean(h2a),mode='same')/(nlen*np.std(h1a)*np.std(h2a))
            
            cct=range(nlen)                                               # Set up the lag axis
            cct=(((np.array(cct)-nlen/2.0)*binning)+(binning/2.0)).tolist()

            pl.figure()
            pl.plot(cct,ccor)
            pl.xlabel('Band '+str(h1)+' lag (s) wrt Band '+str(h2))
            pl.ylabel('Cross-Correlation')
            plot_save(saveplots,show_block)
            print 'Band '+str(h1)+' lag against band '+str(h2)+' (cross-correlation) diagram plotted!'

      else:
         print 'Not enough infiles for cross-correlation!'


   #-----'timeres crosscor' Option----------------------------------------------------------------------------------------

   elif plotopt[:16]=='timeres crosscor':                                 # Plot time-resolved cross-correlation function

      ht=plotopt[16:]                                                     # Collect the xy token from the user

      if nfiles>1:
         if not (ht in ['12','13','21','23','31','32']):                  # Check that the token is 2 long and contains two different characters of the set [1,2,3]

            print 'Invalid command!'
            print ''
            print 'Did you mean...'
            print ''
            print 'HID options:'
            print '* "timeres crosscor21" for 2 against 1 time-resolved cross-correlation'
            if nfiles==3:
               print '* "timeres crosscor32" for 3 against 2 time-resolved cross-correlation'
               print '* "timeres crosscor31" for 3 against 1 time-resolved cross-correlation'

         elif ('3' in ht) and (nfiles<3):

            print 'Not enough infiles for advanced cross-correlation!'    # If token contains a 3 but only 2 infiles are used, abort!

         else:

            try:
               nlen=int(float(raw_input('Length of segments (s): '))/binning)
               if nlen%2!=0: nlen-=1
               frang=float(raw_input('Max lag (s) : '))
            except:
               print 'Not a number!  Aborting!'
               continue

            maxtims=[]
            maxlocs=[]

            n=len(times)/nlen

            matrix=[]                                                     # Set up the 2d matrix
            mtimes=[]                                                     # Set up the time axis

            xn=range(nlen)                                                # Set up the lag axis
            xn=((np.array(xn)-nlen/2.0)*binning).tolist()

            h1=int(ht[0])                                                 # Extract file 1 number
            h2=int(ht[1])                                                 # Extract file 2 number

            #zps=int(frang/binning)-1                                      # Get index of zero in txn array (constructed later)

            for i in range(n):
               y1clip=np.array(ys[h1][nlen*i:nlen*(i+1)])                 # Clip the first lightcurve into chunks
               y2clip=np.array(ys[h2][nlen*i:nlen*(i+1)])                 # Clip the second lightcurve into chunks
               cclip=sig.correlate(y1clip-np.mean(y1clip),y2clip-np.mean(y2clip),mode='same')/(len(y2clip)*np.std(y1clip)*np.std(y2clip))
                                                                          # ^ do the cross-correlation

               cclip=np.array(cclip)[np.array(xn)<=frang]                 # Cut the matrix line to the user requested lag range
               txn=np.array(xn)[np.array(xn)<=frang]
               cclip=cclip[txn>-frang]
               txn=txn[txn>-frang]  

               #print txn[zps]
               #zeroint=ipo.CubicSpline([zps-2,zps-1,zps+1,zps+2],[cclip[zps-2],cclip[zps-1],cclip[zps+1],cclip[zps+2]])(zps)
               #print zeroint

               #pl.figure()
               #pl.plot(txn,cclip)

               #cclip[zps]=zeroint

               #pl.plot(txn,cclip)
               #pl.show(block=True)

               maxc=max(cclip)
               maxloc=txn[cclip.tolist().index(maxc)]

               matrix.append(cclip)                                       # Grow the matrix!
               mtimes.append(nlen*i*binning)
               maxtims.append(nlen*(i+0.5)*binning)
               maxlocs.append(maxloc)

            pl.figure()
            pl.pcolor(np.array(mtimes+[mtimes[-1]+binning*nlen])-binning/2.0,txn,np.array(matrix).T)
            pl.xlabel('Time (s)')
            pl.ylabel('Band '+str(h1)+' lag (s) wrt Band '+str(h2))
            pl.plot(maxtims,maxlocs,':k',label='Peak lag')
            pl.legend()
            plot_save(saveplots,show_block)
            print 'Band '+str(h1)+' lag against band '+str(h2)+' (time-resolved cross-correlation) diagram plotted!'

      else:
         print 'Not enough infiles for cross-correlation!'


   #-----'colxy' Option------------------------------------------------------------------------------------------------

   elif plotopt[:3]=='col':                                                # Plot x/y colour/t

      ht=plotopt[3:]                                                       # Collect the xy token from the user

      if nfiles>1:
         if not (ht in ['12','13','21','23','31','32']):                   # Check that the token is 2 long and contains two different characters of the set [1,2,3]

            print 'Invalid command!'
            print ''
            print 'Did you mean...'
            print ''
            print 'Col/t plot options:'
            print '*"col21" for 2/1 colour'
            print '*"col12" for 1/2 colour'
            if nfiles==3:
               print '*"col32" for 3/2 colour'
               print '*"col23" for 2/3 colour'
               print '*"col31" for 3/1 colour'
               print '*"col13" for 1/3 colour'

         elif ('3' in ht) and (nfiles<3):

            print 'Not enough infiles for advanced Col/t plot!'           # If token contains a 3 but only 2 infiles are used, abort!

         else:

            taxis='Phase' if folded else 'Time (s)'

            h1=int(ht[0])                                                 # Extract numerator file number
            h2=int(ht[1])                                                 # Extract denominator file number
            ht=int(ht)
            pl.figure()
            doplot(times,timese,col[ht],cole[ht],ovr=True)                # Collect colours from col library and plot
            pl.xlabel(taxis)
            pl.ylabel('('+ch[h1]+'/'+ch[h2]+') colour')
            pl.ylim(ymin=0)
            pl.title(fldtxt+'Colour over Time Diagram'+qflav)
            plot_save(saveplots,show_block)
            print 'File'+str(h1)+'/File'+str(h2)+' Colour over Time Diagram plotted!'

      else:
         print 'Not enough infiles for Col/t plot!'


   #-----'ccd' Option--------------------------------------------------------------------------------------------------

   elif plotopt=='ccd':                                                   # Plot Colour-Colour diagram

      if nfiles==3:
         pl.figure()
         doplot(col[31],cole[31],col[21],cole[21])
         pl.xlabel('('+ch[2]+'/'+ch[1]+') colour')
         pl.ylabel('('+ch[3]+'/'+ch[1]+') colour')
         pl.xlim(0,2)
         pl.ylim(0,2)
         pl.title(fldtxt+'Colour-Colour Diagram'+qflav)
         plot_save(saveplots,show_block)
         print 'CCD plotted!'
      else:
         print 'Not enough infiles for CCD!'


   #-----'all' Option--------------------------------------------------------------------------------------------------

   elif plotopt=='all':

      pl.figure()
      if   nfiles==3: colexp=2; rowexp=2; gexp=4                          # If 3 files given, 4 graphs will be plotted in a 2x2 grid
      elif nfiles==2: colexp=1; rowexp=2; gexp=2                          # If 2 files given, 2 grapgs will be plotted in a 2x1 grid
      else:           colexp=1; rowexp=1; gexp=1                          # If 1 file given, only one graph can be plotted
      
      print 'Plotting Lightcurve...'

      taxis='Phase' if folded else 'Time (s)'

      pl.subplot(rowexp,colexp,1)                                         # Create subplot in the first slot
      doplot(times,timese,flux,fluxe,ovr=True)                            # Always plot the lightcurve
      pl.xlabel(taxis)
      pl.ylabel(flux_axis)
      pl.ylim(ymin=0)
      pl.title(fldtxt+'Lightcurve'+qflav)

      if nfiles>1:                                                        # If 2+ files given, plot 2+ file data products

         print 'Plotting Soft Hardness-Intensity Diagram...'

         pl.subplot(rowexp,colexp,2)                                      # Create subplot in the second slot
         doplot(col[21],cole[21],flux,fluxe)                              # Plot Soft HID
         pl.xlim(0,2)
         pl.ylim(0,300)
         pl.ylabel(flux_axis)
         pl.xlabel('('+ch[2]+'/'+ch[1]+') colour')
         pl.title(fldtxt+'Soft HID'+qflav)

      if nfiles==3:                                                       # If 3 files given, plot 3 file data products

         print 'Plotting Hard Hardness-Intensity Diagram...'
         print 'Plotting Colour-Colour Diagram...'

         pl.subplot(rowexp,colexp,3)                                      # Create subplot in the third slot
         doplot(col[31],cole[31],flux,fluxe)                              # Plot Hard HID
         pl.xlim(0,2)
         pl.ylim(0,300)
         pl.ylabel(flux_axis)
         pl.xlabel('('+ch[3]+'/'+ch[1]+') colour')
         pl.title(fldtxt+'Hard HID'+qflav)

         pl.subplot(rowexp,colexp,4)                                      # Create subplot in the fourth slot
         doplot(col[31],cole[31],col[21],cole[21])                        # Plot CCD
         pl.xlim(0,2)
         pl.ylim(0,2)
         pl.ylabel('('+ch[2]+'/'+ch[1]+') colour')
         pl.xlabel('('+ch[3]+'/'+ch[1]+') colour')
         pl.title(fldtxt+'CCD'+qflav)

      print ''
      plot_save(saveplots,show_block)
      print 'All products plotted!'


   #-----'band' Option------------------------------------------------------------------------------------------------

   elif plotopt=='band':                                                  # Plot lightcurve of individual band

      if nfiles==1:
         user_b_band='1'                                                  # Select energy band to plot
      else:
         if nfiles==3:
            is_band_3=', 3'
         else:
            is_band_3=''
         user_b_band=raw_input('Select Energy Band [1, 2'+is_band_3+']: ')

      avail_b_band=['1','2']                                              # Define valid user inputs
      if nfiles==3:
         avail_b_band.append('3')                                         # Add '3' as a valid input if 3 bands present

      if user_b_band in avail_b_band:

         taxis='Phase' if folded else 'Time (s)'
         pl.figure()

         if user_b_band=='1':
            doplot(times,timese,y1[gmask],ye1[gmask],ovr=True,per2=folded) # Plot flux/time using doplot from pan_lib
            b_band_name='Band 1'
         elif user_b_band=='2':
            doplot(times,timese,y2[gmask],ye2[gmask],ovr=True,per2=folded)
            b_band_name='Band 2'
         else:
            b_band_name='Band 3'
            doplot(times,timese,y3[gmask],ye3[gmask],ovr=True,per2=folded)

         pl.xlabel(taxis)                                                 # Format plot
         pl.ylabel(flux_axis)
         pl.ylim(ymin=0)
         pl.title(fldtxt+b_band_name+'Lightcurve'+qflav)
         plot_save(saveplots,show_block)
         print b_band_name+' lightcurve plotted!'


   #-----'bands' Option------------------------------------------------------------------------------------------------

   elif plotopt=='bands':                                                 # Plot lightcurves of individual bands apart

      taxis='Phase' if folded else 'Time (s)'

      pl.figure()
      pl.subplot(nfiles,1,1) 
      doplot(times,timese,y1[gmask],ye1[gmask],ovr=True,per2=folded)      # Plot the lowest band
      pl.xlabel(taxis)
      pl.ylabel(flux_axis)
      pl.title(fldtxt+ch[1]+' Lightcurve'+qflav)

      if nfiles>1:
         pl.subplot(nfiles,1,2)
         doplot(times,timese,y2[gmask],ye2[gmask],ovr=True,per2=folded)   # Plot the second band
         pl.xlabel(taxis)
         pl.ylabel(flux_axis)
         pl.title(fldtxt+ch[2]+' Lightcurve'+flv2)

      if nfiles>2:
         pl.subplot(nfiles,1,3)
         doplot(times,timese,y3[gmask],ye3[gmask],ovr=True,per2=folded)   # Plot the third band
         pl.xlabel(taxis)
         pl.ylabel(flux_axis)
         pl.title(fldtxt+ch[3]+' Lightcurve'+flv3)

      plot_save(saveplots,show_block)
      print 'Banded lightcurves plotted!'


   #-----'xbands' Option------------------------------------------------------------------------------------------------

   # 'Same axes bands'

   elif plotopt=='xbands':                                                # Plot lightcurves of individual bands together

      taxis='Phase' if folded else 'Time (s)'

      pl.figure()
      leg=[ch[1]]                                                         # Create a legend array to populate with channel names
      doplot(times,timese,y1[gmask],ye1[gmask],ovr=True,ft='-b',per2=folded)    # Plot the lowest band
      if nfiles>1:
         doplot(times,timese,y2[gmask],ye2[gmask],ovr=True,ft='-g',per2=folded) # Plot the second band
         leg.append(ch[2])                                                # Append name of second channel to key
      if nfiles>2:
         doplot(times,timese,y3[gmask],ye3[gmask],ovr=True,ft='-r',per2=folded) # Plot the third band
         leg.append(ch[3])                                                # Append name of third channel to key
      pl.legend(leg)                                                      # Create key on plot
      pl.xlabel(taxis)
      pl.ylabel(flux_axis)
      pl.title(fldtxt+'Lightcurve'+qflav)
      plot_save(saveplots,show_block)
      print 'Banded lightcurves plotted!'


   #-----'burst get'--------------------------------------------------------------------------------------------------

   elif plotopt in ['burst get','burstget','burst_get','bursts get','burstsget','bursts_get']:

      if folded:
         print 'Cannot perform burst analsysis on folded data!'
         continue

      while True:
         try:
            iq_lo=float(raw_input('Low Threshold:  '))
            iq_hi=float(raw_input('High Threshold: '))
            assert iq_hi>iq_lo
            assert iq_hi<=100
            assert iq_lo>=0
            break
         except:
            print 'Invalid Entry!  Valid entry is of the form High>Low.'

      bursts={}
      bursts['endpoints']=pan.get_bursts(flux,q_lo=iq_lo,q_hi=iq_hi,just_peaks=False,alg=burst_alg)

      pl.figure()
      doplot(times,timese,flux,fluxe,ovr=True)                            # Plot flux/time using doplot from pan_lib

      t_lo=np.percentile(flux,iq_lo)                                      # Fetch thresholds used in get_bursts
      t_hi=np.percentile(flux,iq_hi)

      col_toggle=True

      for i in bursts['endpoints']:
         if col_toggle:
            col_toggle=False
            burst_colour='#c7c7c7'
         else:
            col_toggle=True
            burst_colour='#e7e7e7'
         pl.axvspan(times[i[0]],times[i[1]], facecolor=burst_colour, edgecolor='none')

      pl.plot([times[0],times[-1]],[t_lo,t_lo],'g')
      pl.plot([times[0],times[-1]],[t_hi,t_hi],'b')
      pl.legend(['Flux','Low Pass Threshold','High Pass Threshold'])
      pl.xlabel('Time (s)')
      pl.ylabel(flux_axis)
      pl.ylim(ymin=0)
      pl.title(fldtxt+'Lightcurve with Bursts Highlighted'+qflav)
      plot_save(saveplots,show_block)

      print ''
      print 'Bursts plotted!'
      print len(bursts['endpoints']),'bursts found!'
      print ''
      print 'Analysing Bursts...'

      bursts['peaks']=[]
      #bursts['trghs']=[]
      bursts['rises']=[]
      bursts['falls']=[]
      bursts['duras']=[]

      for endpoints in bursts['endpoints']:
         burst=flux[endpoints[0]:endpoints[1]]
         btime=times[endpoints[0]:endpoints[1]]
         peak,trough,pk_time,rise_time,fall_time=pan.eval_burst(btime,burst)
         bursts['peaks'].append(peak)
         #bursts['trghs'].append(troughs)
         bursts['rises'].append(rise_time)
         bursts['falls'].append(fall_time)
         bursts['duras'].append(btime[-1]-btime[0])

      print ''
      print 'Analysis Complete!'
      print 'Burst Products now available!'


   #-----'burst peaks'-------------------------------------------------------------------------------------------------

   elif plotopt in ['burst peaks','burstpeaks','burst_peaks','bursts peaks','burstspeaks','bursts_peaks']:

      burstplot('peaks','Peak Heights','cts/s/PCU')


   #-----'burst risetimes'---------------------------------------------------------------------------------------------

   elif plotopt in ['burst risetimes','burstrisetimes','burst_risetimes','bursts risetimes','burstsrisetimes','bursts_risetimes']:

      burstplot('rises','Rise Times','s')


   #-----'burst falltimes'---------------------------------------------------------------------------------------------

   elif plotopt in ['burst falltimes','burstfalltimes','burst_falltimes','bursts falltimes','burstsfalltimes','bursts_falltimes']:

      burstplot('falls','Fall Times','s')


   #-----'burst lengths'-----------------------------------------------------------------------------------------------

   elif plotopt in ['burst lengths','burstlengths','burst_lengths','bursts lengths','burstslengths','bursts_lengths']:

      burstplot('duras','Durations','s')


   #-----'bursts help'-------------------------------------------------------------------------------------------------

   elif plotopt in ['burst help','bursthelp','burst_help','bursts help','burstshelp','bursts_help']:
      print 'Help coming soon.'

   #-----'burst alg'---------------------------------------------------------------------------------------------------

   elif plotopt in ['burst alg','burstalg','burst_alg','bursts alg','burstsalg','bursts_alg']:
      print 'Available Burst-Finding Algorithms:'
      print ' * CUBIC SPLINE'
      print ' * LOAD'
      print ''
      inp_burst_alg=raw_input('Select Burst Algorithm: ').lower()
      if inp_burst_alg in ('cubic spline','load'):
         burst_alg=inp_burst_alg
         print 'Burst-Finding Algorithm set to "'+burst_alg+'"!'
      else:
         print 'Invalid Burst-Finding Algorithm!'


   #-----'lombscargle' Option------------------------------------------------------------------------------------------

   elif plotopt=='lombscargle':

      if folded:                                                          # If data is folded, abort
         print 'Cannot perform Lomb-Scargle on folded data!'
         continue

      if nfiles==1:
         user_scargl_bands='1'                                           # Select energy bands to LombScargle
      else:
         if nfiles==3:
            is_band_3=', 3'
         else:
            is_band_3=''
         user_scargl_bands=raw_input('Select Energy Band [1, 2'+is_band_3+', All]: ').lower()

      ls_st=max(4.0/(times[-1]-times[0]),0.005)
      ls_end=0.5/binning
      lsx=np.arange(ls_st,ls_end,(ls_end-ls_st)/2500.0)

      avail_scargl_bands=['1','2','all']                                  # Define valid user inputs
      if nfiles==3:
         avail_scargl_bands.append('3')                                   # Add '3' as a valid input if 3 bands present

      if user_scargl_bands in avail_scargl_bands:
         if user_scargl_bands=='1':
            lsy=pan.lomb_scargle(times,y1[gmask],ye1[gmask],lsx)          # Perform LombScargle of band 1 using lombscargle function defined in header
            s_band_name='band 1'
         elif user_scargl_bands=='2':
            lsy=pan.lomb_scargle(times,y2[gmask],ye2[gmask],lsx)          # Perform LombScargle of band 2 using lombscargle function defined in header
            s_band_name='band 2'
         elif user_scargl_bands=='3':
            lsy=pan.lomb_scargle(times,y3[gmask],ye3[gmask],lsx)          # Perform LombScargle of band 3 using lombscargle function defined in header
            s_band_name='band 3'
         else:
            lsy=pan.lomb_scargle(times,flux,fluxe,lsx)                    # Perform LombScargle of all bands using lombscargle function defined in header
            s_band_name='all bands'
            if user_scargl_bands!='all':
               print 'Invalid band!  Using all.'

         pl.figure()
         pl.plot(lsx,lsy,'k')                                          # Plot lombscargle
         pl.xlabel('Frequency (Hz)')
         pl.ylabel('Power')
         pl.xlim(0,max(lsx))
         pl.ylim(1,100000)
         pl.yscale('log')
         pl.title('Lomb-Scargle Periodogram of '+s_band_name+qflav)
         plot_save(saveplots,show_block)

         print ''
         print 'Lomb-Scargle Diagram of '+s_band_name+' plotted!'

      else:
         print 'Invalid energy band!'


   #-----'errors' Option-----------------------------------------------------------------------------------------------

   elif plotopt in ['errors','error']:                                    # Toggle Errors

      if es:
         es=False
         print 'Errors suppressed!'
      else:
         es=True
         print 'Errors displayed!'


   #-----'ckey' Option-------------------------------------------------------------------------------------------------

   # 'Colour Key'

   elif plotopt=='ckey':                                                  # Toggle Colour-key

      if cs:
         cs=False
         print 'Colour key suppressed!'
      else:
         cs=True
         print 'Colour key displayed!'


   #-----'lines' Option-------------------------------------------------------------------------------------------------

   elif plotopt=='lines':                                                 # Toggle Delineation

      if ls:
         ls=False
         print 'Plot Lines suppressed!'
      else:
         ls=True
         print 'Plot Lines displayed!'


   #-----'info' Option-------------------------------------------------------------------------------------------------

   elif plotopt=='info':

      dst=times[0]
      det=times[-1]

      print 'PlotDemon.py version',version
      print ''
      print nfiles,'files loaded:'
      print ''
      filn1,loca1=pan.xtrfilloc(file1)
      print 'File 1:'
      print ' Filename       = ',filn1
      print ' Location       = ',loca1
      print ' Mission        = ',mis1
      print ' Object         = ',obsd1[0]
      print ' Obs_ID         = ',obsd1[1]
      if mis1 in ['SUZAKU']:
         print ' Energy         = ',ch[1],'eV'
      else:
         print ' Channel        = ',ch[1]
      print ' Resolution     = ',str(bsz1)+'s'
      print ' BG Subtracted  = ',bsub1
      print ' No. of PCUs    = ',pcus1
      print ' Flavour        = ',flv1
      print ' FITSGenie Ver. = ',v1
      if nfiles>1:
         filn2,loca2=pan.xtrfilloc(file2)
         print ''
         print 'File 2:'
         print ' Filename       = ',filn2
         print ' Location       = ',loca2
         print ' Mission        = ',mis2
         print ' Object         = ',obsd2[0]
         print ' Obs_ID         = ',obsd2[1]
         if mis2 in ['SUZAKU']:
            print ' Energy         = ',ch[2],'eV'
         else:
            print ' Channel        = ',ch[2]
         print ' Resolution     = ',str(bsz2)+'s'
         print ' BG Subtracted  = ',bsub2
         print ' No. of PCUs    = ',pcus2
         print ' Flavour        = ',flv2
         print ' FITSGenie Ver. = ',v2
      if nfiles==3:
         filn3,loca3=pan.xtrfilloc(file3)
         print ''
         print 'File 3:'
         print ' Filename       = ',filn3
         print ' Location       = ',loca3
         print ' Mission        = ',mis3
         print ' Object         = ',obsd3[0]
         print ' Obs_ID         = ',obsd3[1]
         if mis3 in ['SUZAKU']:
            print ' Energy         = ',ch[3],'eV'
         else:
            print ' Channel        = ',ch[3]
         print ' BG Subtracted  = ',bsub3
         print ' Resolution     = ',str(bsz3)+'s'
         print ' No. of PCUs    = ',pcus3
         print ' Flavour        = ',flv3
         print ' FITSGenie Ver. = ',v3
      print ''
      print 'Other Info:'
      print ' Global Flavour = ',flavour
      print ' Obs. Starttime = ',str(tst1)+'s (0.0s)'                     # The start of the observation
      print ' Obs. Endtime   = ',str(oet+tst1)+'s ('+str(oet)+'s)'        # The end of the observation
      print ' Data Starttime = ',str(dst+tst1)+'s ('+str(dst)+'s)'        # The start of the data set (i.e. after GTI considerations and clipping)
      print ' Data Endtime   = ',str(det+tst1)+'s ('+str(det)+'s)'        # The end of the data set
      print ' Bin-size       = ',str(binning)+'s'
      print ' Background     = ',str(bg)+'cts/s/PCU'
      print ' Folded         = ',folded
      if folded:
         print ' Fold Period    = ',period
      print ' Errorbars      = ',es
      print ' Delineated     = ',ls
      print ' Colour-coded   = ',cs


   #-----'reflav' Option-----------------------------------------------------------------------------------------------

   elif plotopt=='reflav':

      print 'Please give a new flavour.'

      try:
         nflavour=raw_input('Flavour: ')
         assert nflavour!=''
         flavour=nflavour
         if flavour=='':
            qflav=''
         else:
            qflav=' "'+flavour+'"'
         print 'Flavour set to "'+flavour+'"'
      except:
         print 'Invalid flavour!  Flavour remains "'+flavour+'"'


   #-----'help' Option-------------------------------------------------------------------------------------------------

   elif plotopt in ['help','?']:                                          # Display instructions

      print 'Instructions:'
      print ''

      give_inst()                                                         # Re-call the instructions list, defined as the get_inst() function in initialisation


   #-----'quit' Option-------------------------------------------------------------------------------------------------

   elif plotopt not in ['quit','exit']:                                   # Invalid command if none of the if statements triggered and no 'q' given

      print 'Invalid command!'

   if plotopt not in ['quit','exit']:
      print ''
      print ' --------------------'


#-----Exiting Interactive Mode-----------------------------------------------------------------------------------------

print ''
print 'Goodbye!'                                           


#-----Footer-----------------------------------------------------------------------------------------------------------

pan.signoff()


