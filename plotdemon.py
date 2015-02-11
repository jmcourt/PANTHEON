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

fac=1 ################################################################# I don't know why I need this factor but it doesn't work without it!!!

x1r,y1r,ye1r,tst1,bsz1,gti,pcus1,bg1,flavour,ch1=xtl.plotdld(file1)       # Opening file 1
y1r=fac*y1r/float(pcus1)                                                  # Normalising flux by dividing by the number of active PCUs and the binsize
ye1r=fac*ye1r/float(pcus1)

if nfiles>1:
   file2=args[2]
   if xtl.flncheck(file2,'plotd'):
      x2r,y2r,ye2r,tst2,bsz2,null,pcus2,null,null,ch2=xtl.plotdld(file2)  # Opening file 2
      del null
      y2r=fac*y2r/float(pcus2)                                            # Normalising flux by dividing by the number of active PCUs and the binsize
      ye2r=fac*ye2r/float(pcus2)
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
      y3r=fac*y3r/float(pcus3)                                            # Normalising flux by dividing by the number of active PCUs and the binsize
      ye3r=fac*ye3r/float(pcus3)
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


#-----Fetch Colours----------------------------------------------------------------------------------------------------

def colorget():
   print 'Analysing Data...'
   if nfiles==1:
      flux=y1
      fluxe=ye1
      col21=col21e=col32=col32e=col31=col31e=None
   elif nfiles==2:
      flux,fluxe,col21,col21e=xtl.pdcolex2(y1,y2,ye1,ye2)
      col32=col32e=col31=col31e=None
   elif nfiles==3:
      flux,fluxe,col21,col21e,col32,col32e,col31,col31e=xtl.pdcolex3(y1,y2,y3,ye1,ye2,ye3)
   else:
      print 'Error!  Too much data somehow.'
      exit()
   return flux,fluxe,col21,col21e,col32,col32e,col31,col31e

flux,fluxe,col21,col21e,col32,col32e,col31,col31e=colorget()

print 'Done!'
print ''


#-----User Menu--------------------------------------------------------------------------------------------------------

def give_inst():                                                          # Define printing this list of instructions as a function
   print 'COMMANDS: Enter a command to manipulate data.'
   print ''
   print 'DATA:'
   print '* "rebin" to reset the data and load it with a different binning'
   print ''
   print '1-SET PLOTS:'
   print '* "lc" to plot a simple graph of flux over time'
   if nfiles>1:
      print ''
      print '2-SET PLOTS:'
      print '* "hid21" to plot a hardness-intensity diagram of file2/file1 colour against total flux'
   print ''
   print 'TOGGLE OPTIONS:'
   print '* "errors" to toggle whether to display errors in plots'

give_inst()                                                               # Print the list of instructions
print ''
print ' --------------------'


#-----Setting up plot environment--------------------------------------------------------------------------------------

plotopt=''
es=True
cs=False
ls=False

def doplot(x,xe,y,ye,xax,yax,ttl):

   if es:
      pl.errorbar(x,y,xerr=xe,yerr=ye)
   else:
      pl.plot(x,y)
   pl.xlabel(xax)
   pl.ylabel(yax)
   pl.title(ttl)
   pl.show(block=False)


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

      print 'Binning complete!'
      print ''

      flux,fluxe,col21,col21e,col32,col32e,col31,col31e=colorget()        # Re-get colours
      print 'Done!'
      print ''


   #-----'lc' Option---------------------------------------------------------------------------------------------------

   elif plotopt=='lc':                                                    # Plot lightcurve

      doplot(x1,zeros(len(x1)),flux,fluxe,'Time (s)','Flux (counts)','t')
      print ''
      print 'Lightcurve plotted!'


   #-----'hid21' Option------------------------------------------------------------------------------------------------

   elif plotopt=='hid21':                                                 # Plot 2/1 HID

      pl.errorbar(x1,col21,yerr=col21e)
      pl.show(block=False)


   #-----'Errors' Option-----------------------------------------------------------------------------------------------

   elif plotopt=='errors':                                                # Toggle Errors

      if es:
         es=False
         print 'Errors suppressed!'
      else:
         es=True
         print 'Errors displayed!'


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


