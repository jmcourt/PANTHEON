#!/usr/bin/python

# |----------------------------------------------------------------------|
# |------------------------------PLOT DEMON------------------------------|
# |----------------------------------------------------------------------|

# Call as ./plotdemon.py GRAPH_TYPE FILE1 [FILE2] [FILE3] BINNING [FLAVOUR] [PLOTSAVE]

# Takes 1-3 csv files of the form [ID TIME FLUX FLUX_ERROR] with a 3-row header and plots relevant astrometric plots
#
# Arguments:
#
#  GRAPH_TYPE
#   The type of graph(s) the user wants to output.  
#   Valid plot types are:
#
#   lc: Lightcurve
#    A plot of total flux of all files entered against time
#
#   shid: Soft Hardness Intensity Diagram
#    A plot of File2/File1 flux (Soft Colour) against total flux
#
#   hhid: Hard Hardness Intensity Diagram
#    A plot of File3/File1 flux (Hard Colour) against total flux
#
#   ccd: Colour-Colour Diagram
#    A plot of File3/File1 flux (Hard Colour) against File2/File1 flux (Soft Colour)
#
#   all: All
#    All of the above
#
#   all2: Lightcurve and SHID
#    Plots all graph defined for only 2 input spectra; namely LC and SHID
#
#   bands: Lightcurves per bin
#    Plots a separate lightcurve for each of 3 files passed into the script
#
#   bands2: Lightcurves per bin for 2 bins
#    Plots a separate lightcurve for each of 2 files passed into the script
#
#  FILE1
#   The absolute path to the first file to be used (generally the lowest energy band)
#
#  [FILE2]
#   Required if and only if 'lc' not used for GRAPH_TYPE.  The absolute path to the second file to be
#   used, (2nd lowest energy band)
#
#  [FILE3]
#   Required if and only if 'all', 'hhid', 'ccd' or 'bands' used for GRAPH_TYPE.  The absolute path to
#   the third file to be used, (highest energy band)
#
#  BINNING
#   The size, in seconds, of bins into which data will be sorted.
#
#  [FLAVOUR]
#   Strictly optional.  A string to be appended to graph titles and saved plot file name, aids with
#   identification and prevents accidental overwrite
# 
#  [PLOTSAVE]
#   Strictly optional, and can only be set if FLAVOUR is also specified.  A string to pass various options
#   into the the script as listed below.  Characters can be given in any order.  Plotsave defaults to '-'
#   if not set by user.
#
#   Meaningful character inputs are:
#    c - Colour    : colours the first five points of CCDs/HIDs red, yellow, green, cyan, blue respectively;
#                     useful for working out which direction a loop is executed in when used alongside l
#
#    e - Error      : create all graph(s) with error bars
#    f - Fold       : attempt to find periodicity in lightcurve data and hence fold it, then plot the folded
#                      data.  Plots two periods for the lightcurve.
#    l - deLineate  : joins plot points on any CCDs and HIDs to show how they progress in time
#    s - Save       : save graph(s), else display them
#    u - Unity      : enable 'unity mode' such that all graphs are created as subplots on a master plot
#    x - eXtract    : use just part of the data; prompt user to select what part to use
#
#

#-----Importing Modules------------------------------------------------------------------------------------------------

import sys,os
import pylab as pl

from math import floor, log10, sqrt
from numpy import array
from numpy import append as npappend                                      # Importing numpy append as npappend to avoid confusion with in-built append function
from xtele_lib import argcheck, binify, foldify, smfold, tnorm            # Custom module xtele_lib is included with the xtelextract package


#-----Welcoming Header-------------------------------------------------------------------------------------------------

print ''
print '-------Running Plot Demon: J.M.Court, 2014------'
print '' 


#-----Checking Validity of Arguments, Fetching Filenames---------------------------------------------------------------

infiles=1                                                                 # Smallest allowable number of input files
bargs=3                                                                   # Smallest allowable number of non-file arguments (Graph_Type, Binning and the function call)

args=sys.argv                                                             # Fetching arguments; softest energy band first please
argcheck(args,infiles+bargs)                                              # Check enough arguments have been given

grpl=args[1]
inf1=args[2]

if grpl not in ['lc','shid','hhid','ccd','all','all2','bands','bands2']:  # Checking graph type specified is valid, telling programme how many infiles to expect
   print 'Please enter one of the following as graph'\
      ' type: lc, shid, ccid, ccd, all, all2.'
   exit()
if grpl in ['hhid','ccd','all','bands','shid','all2','bands2']:           # Assume 2 input files is user requests shid or all2                                  
   infiles=2
   argcheck(args,infiles+bargs)
   inf2=args[3]

if grpl in ['hhid','ccd','all','bands']:                                  # Assume 3 input files is user requests hhid, ccd, all or bands3
   infiles=3
   argcheck(args,infiles+bargs)
   inf3=args[4]

binning=float(args[infiles+2])                                            # Fetching binning number


#-----Fetching Flavour and Plotsave------------------------------------------------------------------------------------

psdef='-'                                                                 # Default value for plotsave
plotsave=psdef

if len(args) in [infiles+bargs+1,infiles+bargs+2]:                        # Checking for the two optional inputs
   flavour=str(args[infiles+bargs])                                       # Allowing flavour to be set
   print "Flavour is '"+str(args[infiles+bargs])+"'"
   if len(args)==infiles+bargs+2:                                         # If only 1 extra argument is given, it is interpreted as Flavour
      plotsave=args[-1]
      print "Plotsave='"+plotsave+"':"
                                                  
   else:
      print 'Plotsave='+"'"+psdef+"' (default)"
elif len(args)>infiles+bargs+2:                                           # Preventing flavour overload!
   print 'Too many args!  No flavour used'
   print 'Plotsave='+"'"+psdef+"' (default)"
   flavour=''
else:
   print 'No flavour given'
   print 'Plotsave='+"'"+psdef+"' (default)"
   flavour=''
print ''


#-----Deciphering Plotsave, setting switches---------------------------------------------------------------------------

cs='c' in plotsave                                                        # Enables coloured mode to show direction of loop execution
es='e' in plotsave                                                        # Switch which causes error bars to be displayed on all plots
fs='f' in plotsave                                                        # Enables folding mode, in which code attempts to fold the data.  Requires user input
ls='l' in plotsave                                                        # Switch which causes adjacent-in-time data points on CCDs/HIDs to be joined with lines
ss='s' in plotsave                                                        # Switch which saves plots instead of displaying them
us='u' in plotsave                                                        # Enable unity mode if 'u' given as char of plotsave
xs='x' in plotsave

if cs: print '* First 5 data points will be colour-coded'
if es: print '* Error bars will be shown'
if fs: print '* Data will be folded'
if ls: print '* CCD/HID data points will be connected by lines'
if ss: print '* Plot will be saved instead of shown'
if us: print '* Graphs will be combined into one image'
if xs: print '* Data will be clipped to user specified range'

print ''

#-----Setting up empty Lists-------------------------------------------------------------------------------------------

x1=[];y1=[];ye1=[]                                                        # Setting up empty arrays to append file1 data to
x2=[];y2=[];ye2=[]                                                        # Setting up empty arrays to append file2 data to
x3=[];y3=[];ye3=[]                                                        # Setting up empty arrays to append file3 data to

time=[]                                                                   # Setting up empty array to append Time data to
flux=[];fluxe=[]                                                          # Setting up empty arrays to append Flux data to
scol=[];scole=[]                                                          # Setting up empty arrays to append Soft Colour data to
hcol=[];hcole=[]                                                          # Setting up empty arrays to append Hard Colour data to


#-----Opening Files, Extracting Raw Data-------------------------------------------------------------------------------

def getfile(name,x,y,ye):                                                 # Defining a script to open a 4 column csv file and paste its [2:4] columns into 3 arrays
   rownum=0
   csvf=open(name)                                                        # Open the file
   for row in csvf.readlines():                                           # Read the file
      if rownum>3:                                                        # Ignore the first 3 lines (ie a 3-line header)
         null,xi,yi,yei=row.split()                                       # Column1>null, Column2>xi, Column3>yi, Column4>yei
         x.append(float(xi))                                              # Append the values found above to the arrays passed into the function
         y.append(float(yi))
         ye.append(float(yei))
      rownum+=1

getfile(inf1,x1,y1,ye1)                                                   # Extracting file 1 

if infiles>1:  
   getfile(inf2,x2,y2,ye2)                                                # Extracting file 2, if given

if infiles==3:
   getfile(inf3,x3,y3,ye3)                                                # Extracting file 3, if given

totres=(x1[-1]-x1[0])/(len(x1)-1)                                         # Determine the resolution of the data in seconds
totres=round(totres, -int(floor(log10(totres))-4))                        # Round this resolution to 5 s.f. to remove error caused by subtraction of large numbers

if totres>binning:

   binning=totres                                                         # Prevent overbinning


#-----Clipping Data----------------------------------------------------------------------------------------------------

if xs:                                                                    # If eXtract mode is enabled...

   print 'Please select Data_ID range over which to plot:'
   pl.figure()
   pl.plot(range(0,len(y1)),y1,'-k')                                      # Plot simple lightcurve of file 1 for user's information
   pl.xlabel('Data_ID')
   pl.ylabel('Flux')
   pl.show(block=False)

   try:
      minif=int(raw_input('Min_ID: '))                                    # Allow the user to select the min of the range over which Fourier analysis will take place
      if minif not in range(0,len(y1)):
         print 'Out of range!  Setting to 0'
         minif=0                                                          # If minif given is out of range, set it to zero
   except: minif=0                                                        # If minif is not given as a number, or left blank, set it to zero

   try:
      maxif=int(raw_input('Max_ID: '))                                    # Allow the user to select the max of the range over which Fourier analysis will take place
      if maxif not in range(0,len(y1) or maxif<minif):
         print 'Out of range!  Setting to '+str(len(y1)-1)
         maxif=len(y1)-1                                                  # If maxif given is out of range, set it to the id of the last data point in y1
   except: maxif=len(y1)-1                                                # If maxif is not given as a number, or left blank, set it to the id of the last data point in y1

   x1=x1[minif:maxif]                                                     # Crop the arrays obtained from file 1
   y1=y1[minif:maxif]
   ye1=ye1[minif:maxif]
   if infiles>1:
      x2=x2[minif:maxif]                                                  # Crop the arrays obtained from file 2
      y2=y2[minif:maxif]
      ye2=ye2[minif:maxif]
      if infiles==3:
         x3=x3[minif:maxif]                                               # Crop the arrays obtained from file 3
         y3=y3[minif:maxif]
         ye3=ye3[minif:maxif]

   pl.close('all')
   print ''


#-----Folding Data-----------------------------------------------------------------------------------------------------

if fs:                                                                    # Only fold data if asked to with an 'f' in plotsave

   valp=False                                                             # Keep asking user until they give a sensible period
   while valp==False:
      try:
        period=float(raw_input('Input period to fold over: '))            # Fetch period from user
        valp=True
      except:
        print "Invalid period!"                                           # Keep trying until they give a sensible input

   x1,y1,ye1=smfold(x1,y1,ye1,period,totres,'1')                          # Fold using smfold function from xtele_lib
      
   if infiles>1:
      x2,y2,ye2=smfold(x2,y2,ye2,period,totres,'2')

      if infiles>2:
         x3,y3,ye3=smfold(x3,y3,ye3,period,totres,'3')

   print 'Folding Complete!'
   print ''


#-----Using Binning Subscript------------------------------------------------------------------------------------------

print 'Bin size='+str(binning)+'s'  
print 'Binning File 1...'
x1,y1,ye1=binify(x1,y1,ye1,binning)                                       # Bin File 1 using 'binify' in xtel_lib
if infiles>1:
   print 'Binning File 2...'
   x2,y2,ye2=binify(x2,y2,ye2,binning)                                    # Bin File 2 using 'binify' in xtel_lib
   if infiles>2:
      print 'Binning File 3...'
      x3,y3,ye3=binify(x3,y3,ye3,binning)                                 # Bin File 3 using 'binify' in xtel_lib

print 'Binning complete!'
print ''


#-----Turning Raw Data into Fluxes and Colours-------------------------------------------------------------------------

breaker=0                                                                 # A switch to allow breaking of nested loops

if infiles>1 and (not grpl in ['bands','bands2']):                        # Procedure to get flux and colour from the extracted data

   print 'Matching timestamps:'
   kmax=lmax=0;                                                           # Allowing all 3 iterative lists to start at element 0
   rowmatches=0                                                           # Tally of row matches
   for i in range(len(x1)):                                               # Iterate over all rows of file 1
      t1=x1[i];f1=y1[i];e1=ye1[i]
      for k in range(kmax,len(x2)):                                       # Iterate over all rows of file 2
         t2=x2[k];f2=y2[k];e2=ye2[k]
         if breaker==1 or t2>t1:                                          # Break the loop if the match was already found
            breaker=0                                                     # Reset breaker
            break
         if t2>t1:                                                        # If t3 overshoots t2, break
            kmax=k                                                        # Start at this value of t2 next iteration
            break
         if infiles>2:                                                    ## CCD,HHID,ALL: Fills Time,Flux,Scol,Hcol arrays only for timestamps present in all three files

            if t1==t2:                                                         # If timestamp in 2 matches stamp in 1, also try to match 3

               for l in range(lmax,len(x3)):
                  t3=x3[l];f3=y3[l];e3=ye3[l]
                  if t2==t3:                                                   # Search file 3 for timestamp matching that in file1[i]
                     kmax=k+1;lmax=l+1                                         # Tells the code where to start off searching files 2 and 3 on the next iteration
                     rowmatches+=1                                             # Tally of how many rows match up
                     if f1<=0 or f2<=0 or f3<=0:
                        break                                                  # Break if any flux is zero, as colours will give div0 errors
                     time.append(t1)                                           # Fill time array
                     flux.append(f1+f2+f3)                                     # Fill flux array
                     scol.append(f2/f1)                                        # Fill soft colour array
                     hcol.append(f3/f1)                                        # Fill hard colour array
                     fluxe.append(sqrt((e1**2)+(e2**2)+(e3**2)))               # Calculate flux error
                     scole.append((f2/f1)*sqrt( ((e1/f1)**2)+((e2/f2)**2) ))   # Calculate soft colour error
                     hcole.append((f3/f1)*sqrt( ((e1/f1)**2)+((e3/f3)**2) ))   # Calculate hard colour error
                     breaker=1                                                 # Switch to tell Python to also break the other loop when the match has been found
                     break
                  elif t3>t2:                                                  # If t3 overshoots t2, break
                     lmax=l                                                    # Start at this value of t3 next iteration
                     break

         else:                                                            ## SHID: Fills Time,Flux,Scol arrays only for timestamps present in both files

            if t1==t2:                                                         # Search file 2 for timestamp matching that in file1[i]
               kmax=k+1                                                        # Tells the code where to start of searching file 2 on the next iteration
               if f1<=0 or f2<=0:
                  break                                                        # Break if any flux is zero, as colours will give div0 errors
               rowmatches+=1                                                   # Tally of how many rows match up
               time.append(t1)                                                 # Fill time array
               flux.append(f1+f2)                                              # Fill flux array
   	       scol.append(f2/f1)                                              # Fill soft colour array
               fluxe.append(sqrt((e1**2)+(e2**2)))                             # Calculate flux error
               scole.append((f2/f1)*sqrt( ((e1/f1)**2)+((e2/f2)**2) ))         # Calculate soft colour error
               break
   print rowmatches,'matches found of a possible',max(len(x1),len(x2),len(x3))

else:                                                                     ## LC: If only one file entered, set Flux and Time as the flux and time from that one file.

   time=array(x1)                                                              # Use x values from file 1 as Time array
   flux=array(y1)                                                              # Use y values from file 1 as Flux array
   fluxe=array(ye1)                                                            # USe y error values from file 1 as Flux Error array

#-----Converting Lists-------------------------------------------------------------------------------------------------

flux=array(flux)                                                          # Redefining each list as a numpy array so they behave better with pylab
fluxe=array(fluxe)
scol=array(scol)
scole=array(scole)
hcol=array(hcol)
hcole=array(hcole)
time=array(time)

if grpl in ['bands','bands2']:                                            # Putting raw data lists into arrays for bands mode
   x1=tnorm(x1,binning)
   y1=array(y1)
   ye1=array(ye1)
   x2=tnorm(x2,binning)
   y2=array(y2)
   ye2=array(ye2)
   if grpl=='bands':
      x3=tnorm(x3,binning)
      y3=array(y3)
      ye3=array(ye3)

y2=array(y2)
y3=array(y3)

ye2=array(ye2)

if not fs:                                                                # This is already done when fold mode is on
   time=tnorm(time,binning)

print str(time[-1])+'s of Data matched!'
print ''


#-----Plot Counter-----------------------------------------------------------------------------------------------------

pos=1                                                                     # The position ID of the first graph to be plotted

if grpl=='all':                                                           # Calculating the number of plots to produce
   expplot=4                                                              # Expect 4 plots if user chose 'all'
elif grpl=='bands':
   expplot=3                                                              # Expect 3 plots if user chose 'bands'
elif grpl in ['all2','bands2']:
   expplot=2                                                              # Expect 2 plots if user chose 'all2' or 'bands2'
else:
   expplot=1                                                              # Expect 1 plots if user chose 'ls, shid, hhid, ccd'


#-----Plot Maker-------------------------------------------------------------------------------------------------------

def doplot(typg,xname,x,xe,yname,y,ye,pos,form='-k'):                          # Define the plotting script as a function for clarity

   if (not us) or pos==1:
      pl.figure()                                                         # Start a new figure for every plot if unity is off, and only for the first plot otherwise

   if us:                                                                 # If unity tag selected
      pl.subplot(rowp,colp,pos)                                           # Place current plot in the next available space on the plot grid

   if fs:
      isfold='Folded '                                                    # Add 'folded' to graph titles if data is folded
      if 'Lightcurve' in typg:                                            # If plotting a lightcurve of folded data, plot two periods instead of one
         pl.plot([x[-1],x[-1]],[min(y),max(y)*1.05],'k:')                 # Plot dotted line to show where period repeats
         x=npappend(x,x+x[-1]+binning)                                    # Double the length of the timespan being plotted over
         y=npappend(y,y)                                                  # Duplicate the y values and errors to populate the second half of the timespan
         ye=npappend(ye,ye)
         numperi=' (2 periods)'                                           # Make sure lightcurve is labelled as representing 2 periods
      else:
         x=npappend(x,x[0])                                               # For CCDs, HIDs, repeat first point at the end of the array to complete loops
         y=npappend(y,y[0])
         xe=npappend(xe,xe[0])
         ye=npappend(ye,ye[0])
         numperi=''
   else:
      isfold=''
      numperi=''

   print 'Plotting',typg,'graph'
   if es:                                                                 # If plotsave contains character 'e' show errors, else suppress
      pl.errorbar(x,y,xerr=xe,yerr=ye,fmt=form)
   else:
      pl.plot(x,y,form)
   if cs and 'Lightcurve' not in typg:                                    # If coloured mode on, colour first 6 data points
      if len(x)<5:                                                        # Abort if less than 6 data points present
         print 'Not enough data to colour!'
      else:
         pl.plot(x[0],y[0],'or')                                          # Plot a round marker over each of the first five points with colour ascending red->blue
         pl.plot(x[1],y[1],'oy')
         pl.plot(x[2],y[2],'og')
         pl.plot(x[3],y[3],'oc')
         pl.plot(x[4],y[4],'ob')
   pl.xlabel(xname)
   pl.ylabel(yname)
   pl.title(flavour+' '+isfold+typg+numperi)

   if us:                                                                 # If in unity mode:
      if pos==expplot:                                                    # If more plots expected, do nothing
         if ss:                                                           # If plotsave contains character 's' save plot, else show
            pl.savefig(flavour+'_all.png')
            print 'Plot saved as '+flavour+'all.png'
         else:
            pl.show(block=True)                                           # Block script to keep graph open
      
   else:                                                                  # If not in unity mode:

      if ss:                                                              # If plotsave contains character 's' save plot, else show
         pl.savefig(flavour+'_'+typg+'.png')
         print 'Plot saved as '+flavour+'_'+typg+'.png'
      else:
         if pos==expplot:
            pl.show(block=True)                                           # If number of plots created = number of plots selected, block script to keep graph open
         else:
            pl.show


#-----Using Plotmaker--------------------------------------------------------------------------------------------------

if us:                                                                    # Calculating size of subplot if UNITY tag is given
   if expplot == 4:
      colp=2                                                              # 2x2 Grid for 4 graphs
      rowp=2
   elif expplot == 3:
      colp=1                                                              # 1x3 Grid for 3 graphs
      rowp=3
   elif expplot == 2 and grpl == 'bands2':
      colp=1                                                              # 1x2 Grid for 2 graphs if plotting 'Bands2'
      rowp=2
   elif expplot == 2:
      colp=2                                                              # 2x1 Grid for 2 graphs otherwise
      rowp=1
   else:
      colp=1                                                              # 1x1 Grid for 1 graph
      rowp=1

pl.rc('lines', linestyle='None')                                          # Removing lines from all scatter plots

if ls:
   clpform='-ok'                                                          # Add lines to CCDs if asked to by user
else:
   clpform='ok'

tmlabel=r'Time (s)'                                                       # Setting labels for axes to represent Time, Flux, and Colours
fxlabel=r'Flux (cts/s)'
sclabel=r'Soft Colour'
hclabel=r'Hard Colour'

if grpl in ['lc','all','all2']:
   doplot('Lightcurve',tmlabel,time,0,fxlabel,flux,fluxe,pos)            # Plot the Lightcurve, overriding pl.rc and letting there be lines
   pos+=1

if grpl in ['shid','all','all2']:
   doplot('Soft HID',sclabel,scol,scole,fxlabel,flux,fluxe,pos,clpform)       # Plot the Soft HID
   pos+=1

if grpl in ['hhid','all']:
   doplot('Hard HID',hclabel,hcol,hcole,fxlabel,flux,fluxe,pos,clpform)       # Plot the Hard HID
   pos+=1

if grpl in ['ccd','all']:
   doplot('CCD',sclabel,scol,scole,hclabel,hcol,hcole,pos,clpform)            # Plot the CCD

if grpl in ['bands','bands2']:
   doplot('Lightcurve, Band 1','Time',x1,0,'Flux',y1,ye1,pos)            # Plot the Lightcurve of the first bin, overriding pl.rc and letting there be lines
   pos+=1
   doplot('Lightcurve, Band 2','Time',x2,0,'Flux',y2,ye2,pos)            # Plot the Lightcurve of the second bin, overriding pl.rc and letting there be lines
   pos+=1

if grpl == 'bands':
   doplot('Lightcurve, Band 3','Time',x3,0,'Flux',y3,ye3,pos)            # Plot the Lightcurve of the third bin, overriding pl.rc and letting there be lines


#-----Footer-----------------------------------------------------------------------------------------------------------

print ''
print '------------------------------------------------'
print ''


