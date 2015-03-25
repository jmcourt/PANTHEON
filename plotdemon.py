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

minbin=0.015625                                                           # The minimum bin size the code is allowed to attempt to use.  This can prevent long hang-ups
version=3.2                                                               # The version of PlotDemon


#-----Welcoming Header-------------------------------------------------------------------------------------------------

print ''
print '-------Running Plot Demon: J.M.Court, 2014------'
print ''


#-----Importing Modules------------------------------------------------------------------------------------------------

try:

   import sys,os
   import pylab as pl
   import pan_lib as pan

   from math import floor, log10, sqrt
   from numpy import array, ones, zeros
   from numpy import append as npappend                                   # Importing numpy append as npappend to avoid confusion with in-built append function

except:

   print 'Modules missing!  Aborting!'
   print ''
   print '------------------------------------------------'
   print ''
   exit()


#-----Opening Files----------------------------------------------------------------------------------------------------

args=sys.argv                                                             # Fetching arguments; softest energy band first please
pan.argcheck(args,2)

try:
   float(args[-1])                                                        # If the final argument can be converted to integer, assume user intends it as a binning 
   isbininp=True                                                          # "IS BINsize given as an INPut?"
except:
   isbininp=False

nfiles=len(args)-isbininp-1                                               # Fetch number of files the user is attempting to plot (total args minus one or two iff binsize given)
if nfiles>3: nfiles=3

file1=args[1]
pan.flncheck(file1,'plotd')

ch={}                                                                     # Save channel info in a library

print 'Opening',file1                                                     # Opening file 1
x1r,y1r,ye1r,tst1,bsz1,gti,pcus1,bg,flv1,ch[1],mis1,obsd1,v1=pan.plotdld(file1)
y1r=y1r/float(pcus1)                                                      # Normalising flux by dividing by the number of active PCUs and the binsize
ye1r=ye1r/float(pcus1)

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

flavour=flv1
if flavour=='':
   qflav=''
else:
   qflav=' "'+flavour+'"'

if nfiles>1:
   file2=args[2]
   if pan.flncheck(file2,'plotd'):
      print 'Opening',file2                                               # Opening file 2
      x2r,y2r,ye2r,tst2,bsz2,null,pcus2,null,flv2,ch[2],mis2,obsd2,v2=pan.plotdld(file2)
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
   if pan.flncheck(file2,'plotd'):
      print 'Opening',file3                                               # Opening file 3
      x3r,y3r,ye3r,tst3,bsz3,null,pcus3,null,flv3,ch[3],mis3,obsd3,v3=pan.plotdld(file3)
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
      pan.signoff()
      exit()

if nfiles>2:                                                              # Checking that start-times of files 1 & 3 match (and thus 2 & 3 also match)
   if tst1!=tst3:
      print 'Starting times for files 1 & 3 do not match!  Aborting!'
      pan.signoff()
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

x1,y1,ye1=pan.binify(x1r,y1r,ye1r,binning)                                # Bin File 1 using 'binify' in pan_lib
if nfiles>1:
   print 'Binning File 2...'
   x2,y2,ye2=pan.binify(x2r,y2r,ye2r,binning)                             # Bin File 2 using 'binify' in pan_lib
   if nfiles>2:
      print 'Binning File 3...'
      x3,y3,ye3=pan.binify(x3r,y3r,ye3r,binning)                          # Bin File 3 using 'binify' in pan_lib

print 'Binning complete!'
print ''

print 'Fetching GTI mask...'

gmask=pan.gtimask(x1,gti)                                                 # A mask to blank values that fall outside of the GTIs

print str(int(100*sum(gmask)/len(gmask)))+'% of data within GTI!'
print ''


#-----Fetch Colours----------------------------------------------------------------------------------------------------

def colorget(verbose=True):                                               # Define colorget to easily re-obtain colours if base data is modified
   if verbose:
      print 'Analysing Data...'
   times=x1[gmask]
   timese=None
   col={}
   cole={}
   if nfiles==1:                                                          # If only one file given, flux and flux_error are just the flux and error of this one file
      flux=y1[gmask]                                                      # Use gmask to clip out the areas outside of GTI
      fluxe=ye1[gmask]
   elif nfiles==2:
      flux,fluxe,col,cole=pan.pdcolex2(y1,y2,ye1,ye2,gmask)               # Get 2/1 and 1/2 colour information using PDColEx in pan_lib
   elif nfiles==3:
      flux,fluxe,col,cole=pan.pdcolex3(y1,y2,y3,ye1,ye2,ye3,gmask)        # Get ALL colour values with 3D PDColEx
   else:
      print 'Error!  Too much data somehow.'                              # This warning should never come up...
      pan.signoff()
      exit()
   return times,timese,flux,fluxe,col,cole

times,timese,flux,fluxe,col,cole=colorget()                               # Use colorget

print 'Done!'
print ''


#-----Setting up plot environment--------------------------------------------------------------------------------------

plotopt=''
es=True                                                                   # Options to keep track of what form the data is in.  'es': with error bars.
cs=False                                                                  # 'cs' with colour key
ls=False                                                                  # 'ls' with delineation
folded=False                                                              # 'folded' has been folded over some period

def doplot(x,xe,y,ye,ovr=False,ft='-k'):                                  # Defining short function to determine whether errorbars are needed on the fly
                                                                          # 'ovr' allows to override colour and line options, so lightcurves can be made differently

   if ovr: formst=ft                                                      # If override given, accept input format; if none given, just plot lines
   elif ls: formst='-ok'                                                  # If deLineate mode on, connect points with lines and mark points
   else: formst='ok'                                                      # If neither deLineate nor override on, just plot points.  No lines here, buddy

   if es:
      pl.errorbar(x,y,xerr=xe,yerr=ye,fmt=formst)                         # Plot errorbar plot if errors turned on
   else:
      pl.plot(x,y,formst)                                                 # Else plot regular graph
   if cs and not ovr:                                                     # If coloured mode on, colour first 5 data points unless override given
      if len(x)<5:                                                        # Abort if less than 5 data points present
         print 'Not enough data to colour!'
      else:
         pl.plot(x[0],y[0],'or')                                          # Plot a round marker over each of the first five points with colour ascending red->blue
         pl.plot(x[1],y[1],'oy')
         pl.plot(x[2],y[2],'og')
         pl.plot(x[3],y[3],'oc')
         pl.plot(x[4],y[4],'ob')


#-----User Menu--------------------------------------------------------------------------------------------------------

def give_inst():                                                          # Define printing this list of instructions as a function
   print 'COMMANDS: Enter a command to manipulate data.'
   print ''
   print 'DATA:'
   print '* "rebin" to reset the data and load it with a different binning'
   print '* "clip" to clip the data'
   print '* "fold" to fold data over a period of your choosing'
   print ''
   print '1+ DATASET PLOTS:'
   print '* "lc" to plot a simple graph of flux over time'
   print '* "animate" to create an animation of the lightcurve as the binning is increased'
   if nfiles>1:                                                           # Only display 2-data-set instructions if 2+ datasets given
      print ''
      print '2+ DATASET PLOTS:'
      print '* "hid21" to plot a hardness-intensity diagram of file2/file1 colour against total flux'
      print '* "hid12" to plot a hardness-intensity diagram of file1/file2 colour against total flux'
      print '* "col21" to plot file2/file1 colour against time.'
      print '* "col12" to plot file1/file2 colour against time.'
      print '* "bands" to plot lightcurves of all bands on adjacent axes.'
      print '* "xbands" to plot lightcurves of all bands on the same axes.'
      print '* "all" to plot all available data products.'
   if nfiles==3:                                                           # Only display 3-data-set instructions if 3 datasets given
      print ''
      print '3 DATASET PLOTS:'
      print '* "hid32" to plot a hardness-intensity diagram of file3/file2 colour against total flux.'
      print '* "hid23" to plot a hardness-intensity diagram of file2/file3 colour against total flux.'
      print '* "hid31" to plot a hardness-intensity diagram of file3/file1 colour against total flux.'
      print '* "hid13" to plot a hardness-intensity diagram of file1/file3 colour against total flux.'
      print '* "col32" to plot file3/file2 colour against time.'
      print '* "col23" to plot file2/file3 colour against time.'
      print '* "col31" to plot file3/file1 colour against time.'
      print '* "col13" to plot file1/file3 colour against time.'
      print '* "ccd" to plot a colour-colour diagram (3/1 colour against 2/1 colour).'
   print ''
   print 'TOGGLE OPTIONS:'
   print '* "errors" to toggle whether to display errors in plots.'
   print '* "lines" to toggle lines joining points in graphs.'
   print '* "ckey" to toggle colour key (red-blue) for the first five points in all plots.'
   print ''
   print 'OTHER COMMANDS:'
   print '* "info" to display a list of facts and figures about the current PlotDemon session.'
   print '* "reflav" to rewrite the flavour text used for graph titles.'
   print '* "help" or "?" to display this list of instructions again.'
   print '* "quit" to quit.'

give_inst()                                                               # Print the list of instructions
print ''
print ' --------------------'


#-----Entering Interactive Mode----------------------------------------------------------------------------------------

while plotopt not in ['quit','exit']:                                     # If the previous command given was not quit, continue

   print ''
   plotopt=raw_input('Give command [? for help]: ')                       # Fetch command from user
   print ''


   #-----'rebin' option------------------------------------------------------------------------------------------------

   if plotopt=='rebin':                                                   # Rebin data

      goodbin=False
      while not goodbin:                                                  # Keep asking until a good response is given
         try:
            binning=float(raw_input("Enter bin size (s): "))              # Ask for binsize in dialogue box
            assert binning>=minbin
            goodbin=True
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

      gmask=pan.gtimask(x1,gti)                                           # Re-establish gmask

      print 'Binning complete!'
      print ''

      times,timese,flux,fluxe,col,cole=colorget()                         # Re-get colours
      folded=False                                                        # Re-allow clipping
      print 'Done!'
      print ''


   #-----'fold' Option-------------------------------------------------------------------------------------------------

   elif plotopt=='fold':                                                  # Fold lightcurve

      goodfold=False                                                      # Keep asking user until they give a sensible period
      while goodfold==False:
         try:
            period=float(raw_input('Input period to fold over (s): '))    # Fetch period from user
            goodfold=True
         except:
            print "Invalid period!"                                       # Keep trying until they give a sensible input

      x1=x1*gmask;y1=y1*gmask;ye1=ye1*gmask                               # Zeroing all data points outside of GTI
      x1,y1,ye1=pan.smfold(x1,y1,ye1,period,binning,'ch. '+ch[1])         # Fold using smfold function from pan_lib

      if nfiles>1:
         x2=x2*gmask;y2=y2*gmask;ye2=ye2*gmask                            # Zeroing all data points outside of GTI
         x2,y2,ye2=pan.smfold(x2,y2,ye2,period,binning,'ch. '+ch[2])      # Fold data of file 2 if present

      if nfiles==3:
         x3=x3*gmask;y3=y3*gmask;ye3=ye3*gmask                            # Zeroing all data points outside of GTI
         x3,y3,ye3=pan.smfold(x3,y3,ye3,period,binning,'ch. '+ch[3])      # Fold data of file 3 if present

      gmask=ones(len(x1),dtype=bool)                                      # Re-establish gmask
      times,timese,flux,fluxe,col,cole=colorget()
      folded=True

      print 'Folding Complete!'
      print ''


   #-----'clip' Option-------------------------------------------------------------------------------------------------

   elif plotopt=='clip':                                                  # Clipping data

      if folded:
         print 'Cannot clip folded data!'

      else:

         print 'Clipping data'
         print ''

         print 'Time range is '+str(x1[0])+'s - '+str(x1[-1])+'s'

         print 'Please choose new range of data:'
         mint,maxt=pan.srinr(times,binning,'time')                        # Fetch new time domain endpoints using srinr function from pan_lib

         print 'Clipping...'

         x1=x1[mint:maxt]                                                 # Clip file 1
         y1=y1[mint:maxt]
         ye1=ye1[mint:maxt]

         if nfiles>1:

            x2=x2[mint:maxt]                                              # Clip file 2
            y2=y2[mint:maxt]
            ye2=ye2[mint:maxt]

         if nfiles==3:

            x3=x3[mint:maxt]                                              # Clip file 3
            y3=y3[mint:maxt]
            ye3=ye3[mint:maxt]

         gmask=pan.gtimask(x1,gti)                                        # Re-establish gmask
         times,timese,flux,fluxe,col,cole=colorget()

         print 'Data clipped!'


   #-----'lc' Option---------------------------------------------------------------------------------------------------

   elif plotopt=='lc':                                                    # Plot lightcurve

      pl.figure()
      doplot(times,timese,flux,fluxe,ovr=True)
      pl.xlabel('Time (s)')
      pl.ylabel('Flux (counts/s/PCU)')
      pl.title('Lightcurve'+qflav)
      pl.show(block=False)
      print 'Lightcurve plotted!'


   #-----'animate' Option----------------------------------------------------------------------------------------------

   elif plotopt=='animate':

      animsloc=raw_input('Folder to save images: ')
      print ''

      if os.path.exists(animsloc):                                        # Create the folder
         print 'Folder "'+animsloc+'" already exists...'
      else:
         print 'Creating folder "'+animsloc+'"...'
         os.makedirs(animsloc)

      here=os.getcwd()

      os.chdir(animsloc)

      animbin=0.0025                                                      # Start with an arbitrarily low binsize

      while animbin<max(bsz1,bsz2,bsz3,minbin):                           # Find lowest allowable binsize of the form 0.01*2^N
         animbin=animbin*2

      anstep=1                                                            # Track the number of steps taken

      dst=times[0]
      det=times[-1]

      while animbin<(det-dst)/4.0:                                        # Set the maximum binsize at one quarter of the observation length

         print "Creating",str(animbin)+"s binned lightcurve"

         x1,y1,ye1=pan.binify(x1r,y1r,ye1r,animbin)                       # Bin File 1 using 'binify' in pan_lib
         if nfiles>1:
            x2,y2,ye2=pan.binify(x2r,y2r,ye2r,animbin)                    # Bin File 2 using 'binify' in pan_lib
            if nfiles>2:
               x3,y3,ye3=pan.binify(x3r,y3r,ye3r,animbin)                 # Bin File 3 using 'binify' in pan_lib

         mina,maxa=pan.srinr(x1,binning,'time',minv=dst,maxv=det)         # Clip each individual lightcurve

         x1=x1[mina:maxa]                                                 # Clip file 1
         y1=y1[mina:maxa]
         ye1=ye1[mina:maxa]

         if nfiles>1:

            x2=x2[mina:maxa]                                              # Clip file 2
            y2=y2[mina:maxa]
            ye2=ye2[mina:maxa]

         if nfiles==3:

            x3=x3[mint:maxt]                                              # Clip file 3
            y3=y3[mint:maxt]
            ye3=ye3[mint:maxt]

         gmask=pan.gtimask(x1,gti)                                        # Re-establish gmask

         times,timese,flux,fluxe,col,cole=colorget(verbose=False)

         if anstep==1:
            if es:
              merr=max(fluxe)
            else:
              merr=0
            maxany=max(flux)+merr                                         # Calculate the range of all plots based on the range of the first plot
            minany=min(flux)-merr
            if minany<0: minany=0

         pl.figure()
         doplot(times,timese,flux,fluxe,ovr=True)
         pl.xlabel('Time (s)')
         pl.ylabel('Flux (counts/s/PCU)')
         pl.title('Lightcurve ('+str(animbin)+'s binning)')
         pl.xlim(dst,det)
         pl.ylim(minany,maxany)
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
      times,timese,flux,fluxe,col,cole=colorget(verbose=False)            # Re-get colours

      print ''
      print "Animation saved to",animsloc+'/animation.gif!'
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
            print '*"hid12" for 2/1 colour'
            print '*"hid21" for 1/2 colour'
            if nfiles==3:
               print '*"hid32" for 3/2 colour'
               print '*"hid23" for 2/3 colour'
               print '*"hid31" for 3/1 colour'
               print '*"hid13" for 1/3 colour'

         elif ('3' in ht) and (nfiles<3):

            print 'Not enough infiles for advanced HID!'                  # If token contains a 3 but only 2 infiles are used, abort!

         else:

            h1=int(ht[0])                                                 # Extract numerator file number
            h2=int(ht[1])                                                 # Extract denominator file number
            ht=int(ht)
            pl.figure()
            doplot(col[ht],cole[ht],flux,fluxe)                           # Collect colours from col library and plot
            pl.ylabel('Flux (counts/s/PCU)')
            pl.xlabel('('+ch[h1]+'/'+ch[h2]+') colour')
            pl.title('Hardness Intensity Diagram'+qflav)
            pl.show(block=False)
            print 'File'+str(h1)+'/File'+str(h2)+' HID plotted!'

      else:
         print 'Not enough infiles for HID!'


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
            print '*"col12" for 2/1 colour'
            print '*"col21" for 1/2 colour'
            if nfiles==3:
               print '*"col32" for 3/2 colour'
               print '*"col23" for 2/3 colour'
               print '*"col31" for 3/1 colour'
               print '*"col13" for 1/3 colour'

         elif ('3' in ht) and (nfiles<3):

            print 'Not enough infiles for advanced Col/t plot!'           # If token contains a 3 but only 2 infiles are used, abort!

         else:

            h1=int(ht[0])                                                 # Extract numerator file number
            h2=int(ht[1])                                                 # Extract denominator file number
            ht=int(ht)
            pl.figure()
            doplot(times,timese,col[ht],cole[ht],ovr=True)                # Collect colours from col library and plot
            pl.xlabel('Time (s)')
            pl.ylabel('('+ch[h1]+'/'+ch[h2]+') colour')
            pl.title('Colour over Time Diagram'+qflav)
            pl.show(block=False)
            print 'File'+str(h1)+'/File'+str(h2)+' Colour over Time Diagram plotted!'

      else:
         print 'Not enough infiles for Col/t plot!'


   #-----'ccd' Option--------------------------------------------------------------------------------------------------

   elif plotopt=='ccd':                                                   # Plot 3/1 HID

      if nfiles==3:
         pl.figure()
         doplot(col[31],cole[31],col[21],cole[21])
         pl.xlabel('('+ch[2]+'/'+ch[1]+') colour')
         pl.ylabel('('+ch[3]+'/'+ch[1]+') colour')
         pl.title('Colour-Colour Diagram'+qflav)
         pl.show(block=False)
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

      pl.subplot(rowexp,colexp,1)                                         # Create subplot in the first slot
      doplot(times,timese,flux,fluxe,ovr=True)                            # Always plot the lightcurve
      pl.xlabel('Time (s)')
      pl.ylabel('Flux (counts/s/PCU)')
      pl.title('Lightcurve'+qflav)

      if nfiles>1:                                                        # If 2+ files given, plot 2+ file data products

         print 'Plotting Soft Hardness-Intensity Diagram...'

         pl.subplot(rowexp,colexp,2)                                      # Create subplot in the second slot
         doplot(col[21],cole[21],flux,fluxe)                              # Plot Soft HID
         pl.ylabel('Flux (counts/s/PCU)')
         pl.xlabel('('+ch[2]+'/'+ch[1]+') colour')
         pl.title('Soft HID'+qflav)

      if nfiles==3:                                                       # If 3 files given, plot 3 file data products

         print 'Plotting Hard Hardness-Intensity Diagram...'
         print 'Plotting Colour-Colour Diagram...'

         pl.subplot(rowexp,colexp,3)                                      # Create subplot in the third slot
         doplot(col[31],cole[31],flux,fluxe)                              # Plot Hard HID
         pl.ylabel('Flux (counts/s/PCU)')
         pl.xlabel('('+ch[3]+'/'+ch[1]+') colour')
         pl.title('Hard HID'+qflav)

         pl.subplot(rowexp,colexp,4)                                      # Create subplot in the fourth slot
         doplot(col[31],cole[31],col[21],cole[21])                        # Plot CCD
         pl.ylabel('('+ch[2]+'/'+ch[1]+') colour')
         pl.xlabel('('+ch[3]+'/'+ch[1]+') colour')
         pl.title('CCD'+qflav)

      print ''
      pl.show(block=False)
      print 'All products plotted!'


   #-----'bands' Option------------------------------------------------------------------------------------------------

   elif plotopt=='bands':                                                 # Plot lightcurves of individual bands apart

      pl.figure()
      pl.subplot(nfiles,1,1) 
      doplot(times,timese,y1[gmask],ye1[gmask],ovr=True)                  # Plot the lowest band
      pl.xlabel('Time (s)')
      pl.ylabel('Flux (counts/s/PCU)')
      pl.title(ch[1]+' Lightcurve'+qflav)
      if nfiles>1:
         pl.subplot(nfiles,1,2)
         doplot(times,timese,y2[gmask],ye2[gmask],ovr=True)               # Plot the second band
         pl.xlabel('Time (s)')
         pl.ylabel('Flux (counts/s/PCU)')
         pl.title(ch[2]+' Lightcurve'+qflav)
      if nfiles>2:
         pl.subplot(nfiles,1,3)
         doplot(times,timese,y3[gmask],ye3[gmask],ovr=True)               # Plot the third band
         pl.xlabel('Time (s)')
         pl.ylabel('Flux (counts/s/PCU)')
         pl.title(ch[3]+' Lightcurve'+qflav)
      pl.show(block=False)
      print 'Banded lightcurves plotted!'


   #-----'xbands' Option------------------------------------------------------------------------------------------------

   # 'Same axes bands'

   elif plotopt=='xbands':                                                # Plot lightcurves of individual bands together

      pl.figure()
      leg=[ch[1]]                                                         # Create a legend array to populate with channel names
      doplot(times,timese,y1[gmask],ye1[gmask],ovr=True,ft='-b')          # Plot the lowest band
      if nfiles>1:
         doplot(times,timese,y2[gmask],ye2[gmask],ovr=True,ft='-g')       # Plot the second band
         leg.append(ch[2])                                                # Append name of second channel to key
      if nfiles>2:
         doplot(times,timese,y3[gmask],ye3[gmask],ovr=True,ft='-r')       # Plot the third band
         leg.append(ch[3])                                                # Append name of third channel to key
      pl.legend(leg)                                                      # Create key on plot
      pl.xlabel('Time (s)')
      pl.ylabel('Flux (counts/s/PCU)')
      pl.title('Lightcurve'+qflav)
      pl.show(block=False)
      print 'Banded lightcurves plotted!'


   #-----'errors' Option-----------------------------------------------------------------------------------------------

   elif plotopt=='errors':                                                # Toggle Errors

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
         print ' Resolution     = ',str(bsz3)+'s'
         print ' No. of PCUs    = ',pcus3
         print ' Flavour        = ',flv3
         print ' FITSGenie Ver. = ',v3
      print ''
      print 'Other Info:'
      print ' Main Flavour   = ',flavour
      print ' Obs. Starttime = ',str(tst1)+'s (0.0s)'                     # The start of the observation
      print ' Obs. Endtime   = ',str(oet+tst1)+'s ('+str(oet)+'s)'        # The end of the observation
      print ' Data Starttime = ',str(dst+tst1)+'s ('+str(dst)+'s)'        # The start of the data set (i.e. after GTI considerations and clipping)
      print ' Data Endtime   = ',str(det+tst1)+'s ('+str(det)+'s)'        # The end of the data set
      print ' Bin-size       = ',str(binning)+'s'
      print ' Background     = ',str(bg)+'cts/s/PCU'
      print ' Folded         = ',folded
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


