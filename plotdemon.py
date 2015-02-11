#!/usr/bin/python

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

#-----Importing Modules------------------------------------------------------------------------------------------------

import sys,os
import pylab as pl
import xtele_lib as xtl

from math import floor, log10, sqrt
from numpy import array, zeros
from numpy import append as npappend                                      # Importing numpy append as npappend to avoid confusion with in-built append function

#-----User-set Parameters----------------------------------------------------------------------------------------------

minbin=0.0625                                                             # The minimum bin size the code is allowed to attempt to use.  This can prevent long hang-ups


#-----Welcoming Header-------------------------------------------------------------------------------------------------

print ''
print '-------Running Plot Demon: J.M.Court, 2014------'
print ''


#-----Opening Files----------------------------------------------------------------------------------------------------

args=sys.argv                                                             # Fetching arguments; softest energy band first please
xtl.argcheck(args,2)

try:
   float(args[-1])                                                        # If the final argument can be converted to integer, assume user intends it as a binning 
   isbininp=True                                                          # "IS BINsize given as an INPut?"
except:
   isbininp=False

nfiles=len(args)-isbininp-1                                               # Fetch number of files the user is attempting to plot (total args minus one or two iff binsize given)
if nfiles>3: nfiles=3

file1=args[1]
xtl.flncheck(file1,'plotd')

x1r,y1r,ye1r,tst1,bsz1,gti,pcus1,bg1,flavour,ch1=xtl.plotdld(file1)       # Opening file 1
y1r=y1r/float(pcus1)                                                      # Normalising flux by dividing by the number of active PCUs and the binsize
ye1r=ye1r/float(pcus1)

if nfiles>1:
   file2=args[2]
   if xtl.flncheck(file2,'plotd'):
      x2r,y2r,ye2r,tst2,bsz2,null,pcus2,null,null,ch2=xtl.plotdld(file2)  # Opening file 2
      del null
      y2r=y2r/float(pcus2)                                                # Normalising flux by dividing by the number of active PCUs and the binsize
      ye2r=ye2r/float(pcus2)
   else:
      nfiles=1
      print 'Warning: File 2 of incorrect file type!'
      print 'Only loading File 1.'
else: x2r=y2r=ye2r=tst2=bsz2=None

if nfiles>2:
   file3=args[3]
   if xtl.flncheck(file2,'plotd'):
      x3r,y3r,ye3r,tst3,bsz3,null,pcus3,null,null,ch3=xtl.plotdld(file3)  # Opening file 3
      del null
      y3r=y3r/float(pcus3)                                                # Normalising flux by dividing by the number of active PCUs and the binsize
      ye3r=ye3r/float(pcus3)
   else:
      nfiles=2
      print 'Warning: File 3 of incorrect file type!'
      print 'Only loading Files 1 & 2.'
else: x3r=y3r=ye3r=tst3=bsz3=None

if nfiles>1:                                                              # Checking that start-times of files 1 & 2 match
   if tst1!=tst2:
      print 'Starting times for files 1 & 2 do not match!  Aborting!'
      xtl.signoff()
      exit()

if nfiles>2:                                                              # Checking that start-times of files 1 & 3 match (and thus 2 & 3 also match)
   if tst1!=tst3:
      print 'Starting times for files 1 & 3 do not match!  Aborting!'
      xtl.signoff()
      exit()


#-----Binning----------------------------------------------------------------------------------------------------------

if isbininp:
   binning=float(args[-1])                                                # Collect binsize input if given
else:
   goodbin=False
   while not goodbin:                                                     # Keep asking until a good response is given
      try:
         binning=float(raw_input("Enter bin size (s): "))                 # Ask for binsize in dialogue box
         goodbin=True
      except:
         print 'Invalid bin size input!'

binning=max(binning,bsz1,bsz2,bsz3,minbin)                                # Prevent overbinning by setting minimum binning to the maximum of the binnings of the files

print ''
print 'Bin size='+str(binning)+'s'  
print 'Binning File 1...'

x1,y1,ye1=xtl.binify(x1r,y1r,ye1r,binning)                                # Bin File 1 using 'binify' in xtel_lib
if nfiles>1:
   print 'Binning File 2...'
   x2,y2,ye2=xtl.binify(x2r,y2r,ye2r,binning)                             # Bin File 2 using 'binify' in xtel_lib
   if nfiles>2:
      print 'Binning File 3...'
      x3,y3,ye3=xtl.binify(x3r,y3r,ye3r,binning)                          # Bin File 3 using 'binify' in xtel_lib

print 'Binning complete!'
print ''

gmask=xtl.gtimask(x1,gti)                                                 # A mask to blank values that fall outside of the GTIs

#-----Fetch Colours----------------------------------------------------------------------------------------------------

def colorget():
   print 'Analysing Data...'
   times=x1[gmask]
   if nfiles==1:
      flux=y1[gmask]
      fluxe=ye1[gmask]
      col21=col21e=col32=col32e=col31=col31e=None
   elif nfiles==2:
      flux,fluxe,col21,col21e=xtl.pdcolex2(y1,y2,ye1,ye2,gmask)
      col32=col32e=col31=col31e=None
   elif nfiles==3:
      flux,fluxe,col21,col21e,col32,col32e,col31,col31e=xtl.pdcolex3(y1,y2,y3,ye1,ye2,ye3,gmask)
   else:
      print 'Error!  Too much data somehow.'
      exit()
   return times,flux,fluxe,col21,col21e,col32,col32e,col31,col31e

times,flux,fluxe,col21,col21e,col32,col32e,col31,col31e=colorget()

print 'Done!'
print ''



#-----Setting up plot environment--------------------------------------------------------------------------------------

plotopt=''
es=True
cs=False
ls=False

def doplot(x,xe,y,ye):                                                    # Defining short function to determine whether errorbars are needed on the fly

   if es:
      pl.errorbar(x,y,xerr=xe,yerr=ye)
   else:
      pl.plot(x,y)


#-----User Menu--------------------------------------------------------------------------------------------------------

def give_inst():                                                          # Define printing this list of instructions as a function
   print 'COMMANDS: Enter a command to manipulate data.'
   print ''
   print 'DATA:'
   print '* "rebin" to reset the data and load it with a different binning'
   print '* "fold" to fold data over a period of your choosing'
   print ''
   print '1+ DATASET PLOTS:'
   print '* "lc" to plot a simple graph of flux over time'
   if nfiles>1:
      print ''
      print '2+ DATASET PLOTS:'
      print '* "hid21" to plot a hardness-intensity diagram of file2/file1 colour against total flux'
      print '* "bands" to plot lightcurves of all bands on adjacent axes'
      print '* "xbands" to plot lightcurves of all bands on the same axes'
   if nfiles=3:
      print ''
      print '3 DATASET PLOTS:'
      print '* "hid32" to plot a hardness-intensity diagram of file3/file2 colour against total flux'
      print '* "hid31" to plot a hardness-intensity diagram of file3/file1 colour against total flux'
      print '* "ccd" to plot a colour-colour diagram (3/1 colour against 2/1 colour)
   print ''
   print 'TOGGLE OPTIONS:'
   print '* "errors" to toggle whether to display errors in plots'
   print ''
   print 'OTHER COMMANDS:'
   print '* "help" or "?" to display this list of Instructions again'
   print '* "quit" to Quit'

give_inst()                                                               # Print the list of instructions
print ''
print ' --------------------'


#-----Entering Interactive Mode----------------------------------------------------------------------------------------

while plotopt not in ['quit','exit']:                                     # If the previous command given was not quit, continue

   print ''
   plotopt=raw_input('Give command [? for help]: ')                       # Fetch command from user


   #-----'Rebin' option------------------------------------------------------------------------------------------------

   if plotopt=='rebin':                                                   # Rebin data

      goodbin=False
      while not goodbin:                                                  # Keep asking until a good response is given
         try:
            binning=float(raw_input("Enter bin size (s): "))              # Ask for binsize in dialogue box
            goodbin=True
         except:
            print 'Invalid bin size input!'

      print 'Binning File 1...'
      x1,y1,ye1=xtl.binify(x1r,y1r,ye1r,binning)                          # Bin File 1 using 'binify' in xtel_lib
      if nfiles>1:
         print 'Binning File 2...'
         x2,y2,ye2=xtl.binify(x2r,y2r,ye2r,binning)                       # Bin File 2 using 'binify' in xtel_lib
         if nfiles>2:
            print 'Binning File 3...'
            x3,y3,ye3=xtl.binify(x3r,y3r,ye3r,binning)                    # Bin File 3 using 'binify' in xtel_lib

      gmask=xtl.gtimask(x1,gti)                                           # Re-establish gmask

      print 'Binning complete!'
      print ''

      times,flux,fluxe,col21,col21e,col32,col32e,col31,col31e=colorget()  # Re-get colours
      print 'Done!'
      print ''


   #-----'lc' Option---------------------------------------------------------------------------------------------------

   elif plotopt=='fold':                                                  # Fold lightcurve
      print 'Not yet implemented!'

   #-----'lc' Option---------------------------------------------------------------------------------------------------

   elif plotopt=='lc':                                                    # Plot lightcurve

      pl.figure()
      doplot(times,zeros(len(times)),flux,fluxe)
      pl.xlabel('Time (s)')
      pl.ylabel('Flux (counts/s/PCU)')
      pl.title('Lightcurve "'+flavour+'"')
      pl.show(block=False)
      print ''
      print 'Lightcurve plotted!'


   #-----'hid21' Option------------------------------------------------------------------------------------------------

   elif plotopt=='hid21':                                                 # Plot 2/1 HID

      pl.figure()
      if nfiles>1:
         doplot(col21,col21e,flux,fluxe)
         pl.ylabel('Flux (counts/s/PCU)')
         pl.xlabel('('+ch2+'/'+ch1+') colour')
         pl.title('Hardness Intensity Diagram "'+flavour+'"')
         pl.show(block=False)
         print ''
         print 'File2/File1 HID plotted!'
      else:
         print 'Not enough infiles for HID!'


   #-----'hid32' Option------------------------------------------------------------------------------------------------

   elif plotopt=='hid32':                                                 # Plot 3/2 HID

      pl.figure()
      if nfiles=3:
         doplot(col32,col32e,flux,fluxe)
         pl.ylabel('Flux (counts/s/PCU)')
         pl.xlabel('('+ch3+'/'+ch2+') colour')
         pl.title('Hardness Intensity Diagram "'+flavour+'"')
         pl.show(block=False)
         print ''
         print 'File3/File2 HID plotted!'
      else:
         print 'Not enough infiles for hard HID!'


   #-----'hid21' Option------------------------------------------------------------------------------------------------

   elif plotopt=='hid31':                                                 # Plot 3/1 HID

      pl.figure()
      if nfiles=3:
         doplot(col31,col31e,flux,fluxe)
         pl.ylabel('Flux (counts/s/PCU)')
         pl.xlabel('('+ch3+'/'+ch1+') colour')
         pl.title('Hardness Intensity Diagram "'+flavour+'"')
         pl.show(block=False)
         print ''
         print 'File3/File1 HID plotted!'
      else:
         print 'Not enough infiles for hard HID!'


   #-----'ccd' Option--------------------------------------------------------------------------------------------------

   elif plotopt=='ccd':                                                   # Plot 3/1 HID

      pl.figure()
      if nfiles=3:
         doplot(col31,col31e,col21,col21e)
         pl.xlabel('('+ch2+'/'+ch1+') colour')
         pl.xlabel('('+ch3+'/'+ch1+') colour')
         pl.title('Colour-Colour Diagram "'+flavour+'"')
         pl.show(block=False)
         print ''
         print 'CCD plotted!'
      else:
         print 'Not enough infiles for CCD!'


   #-----'bands' Option------------------------------------------------------------------------------------------------

   elif plotopt=='bands':                                                 # Plot lightcurves of individual bands apart

      pl.figure()
      pl.subplot(nfiles,1,1) 
      doplot(times,zeros(len(times)),y1[gmask],ye1[gmask])
      pl.xlabel('Time (s)')
      pl.ylabel('Flux (counts/s/PCU)')
      pl.title(ch1+' Lightcurve "'+flavour+'"')
      if nfiles>1:
         pl.subplot(nfiles,1,2)
         doplot(times,zeros(len(times)),y2[gmask],ye2[gmask])
         pl.xlabel('Time (s)')
         pl.ylabel('Flux (counts/s/PCU)')
         pl.title(ch2+' Lightcurve "'+flavour+'"')
      if nfiles>2:
         pl.subplot(nfiles,1,3)
         doplot(times,zeros(len(times)),y3[gmask],ye3[gmask])
         pl.xlabel('Time (s)')
         pl.ylabel('Flux (counts/s/PCU)')
         pl.title(ch3+' Lightcurve "'+flavour+'"')
      pl.show(block=False)
      print ''
      print 'Banded lightcurves plotted!'


   #-----'xbands' Option------------------------------------------------------------------------------------------------

   elif plotopt=='xbands':                                                # Plot lightcurves of individual bands together

      pl.figure()
      leg=[ch1]
      doplot(times,zeros(len(times)),y1[gmask],ye1[gmask])
      if nfiles>1:
         doplot(times,zeros(len(times)),y2[gmask],ye2[gmask])
         leg.append(ch2)
      if nfiles>2:
         doplot(times,zeros(len(times)),y3[gmask],ye3[gmask])
         leg.append(ch3)
      pl.legend(leg)
      pl.xlabel('Time (s)')
      pl.ylabel('Flux (counts/s/PCU)')
      pl.title('Lightcurve "'+flavour+'"')
      pl.show(block=False)
      print ''
      print 'Banded lightcurves plotted!'

   #-----'Errors' Option-----------------------------------------------------------------------------------------------

   elif plotopt=='errors':                                                # Toggle Errors

      print ''

      if es:
         es=False
         print 'Errors suppressed!'
      else:
         es=True
         print 'Errors displayed!'


   #-----'Instructions' Option-----------------------------------------------------------------------------------------

   elif plotopt in ['help','?']:                                          # Display instructions

      print 'Instructions:'
      print ''

      give_inst()                                                         # Re-call the instructions list, defined as the get_inst() function in initialisation


   #-----'Quit' Option-------------------------------------------------------------------------------------------------

   elif plotopt not in ['quit','exit']:                                   # Invalid command if none of the if statements triggered and no 'q' given

      print 'Invalid command!'

   if plotopt not in ['quit','exit']:
      print ''
      print ' --------------------'


#-----Exiting Interactive Mode-----------------------------------------------------------------------------------------

print ''
print 'Goodbye!'                                           


#-----Footer-----------------------------------------------------------------------------------------------------------

xtl.signoff()


