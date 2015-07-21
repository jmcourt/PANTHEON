#! /usr/bin/env python

# |----------------------------------------------------------------------|
# |------------------------------SPEC ANGEL------------------------------|
# |----------------------------------------------------------------------|

# Call as ./specangel.py FILE1 [LBINNING]

# Takes 1 RXTE FITS Event file and produces an interactive spectrogram
#
# Arguments:
#
#  FILE1
#   The absolute path to the file to be used.
#
#  [LBINNING]
#   Optional- the logarithmic binning factor 'x'; frequency data will be binned into bins which have
#   their left-hand edges defined by the formula 10**(ix) for integer i.
#
#

#-----User-set Parameters----------------------------------------------------------------------------------------------

lplresdefault=0.005                                                       # The best resolution, in log10 space, in which the data will be analysed
version=4.2                                                               # The version of SpecAngel


#-----Welcoming Header-------------------------------------------------------------------------------------------------

print ''
print '-------Running Spec Angel: J.M.Court, 2015------'
print ''


#-----Importing Modules------------------------------------------------------------------------------------------------

try:

   import sys,os
   import pylab as pl
   import pan_lib as pan
   from astropy.io import fits
   from numpy import array, arange, mean, meshgrid, transpose, zeros
   from numpy import exp, histogram, linspace, log10, nonzero, sqrt
   from numpy import append as npappend                                   # Importing numpy append as npappend to avoid confusion with in-built append function
   from numpy import min as npmin
   from numpy import max as npmax
   from numpy import sum as npsum
   from scipy import delete
   from scipy.fftpack import fft
   from scipy.stats import spearmanr
   from scipy.optimize import brentq

except:

   print 'Modules missing!  Aborting!'
   print ''
   print '------------------------------------------------'
   print ''
   exit()


#-----Checking Validity of Arguments, Fetching File--------------------------------------------------------------------

args=sys.argv
pan.argcheck(args,2)                                                      # Must give at least 2 args (Filename and the function call)

filename=args[1]                                                          # Fetch file name from arguments
pan.flncheck(filename,'speca')


#-----Extracting data from file-----------------------------------------------------------------------------------------

print 'Opening '+str(filename)                                            # Use SpecaLd from pan_lib to load data from file
loadmatrix,good,rates,prates,trates,phcts,bg,bsz,foures,bgest,flv,cs,mis,obsd,wtype,slide,binfac,v=pan.specald(filename)
flavour=flv
if flavour=='':
   qflav=''
else:
   qflav=' "'+flavour+'"'


#-----Initially normalising data----------------------------------------------------------------------------------------

print ''
print 'Normalizing and binning...'

if len(args)>2:                                                           # Check for lplres input, else request one
   try:
      lplres=float(args[2])
      assert lplres>0
   except:
      lplres=lplresdefault
else:
   try:
      lplres=float(raw_input('Logarithmic binning factor: '))
   except:
      lplres=lplresdefault

if len(args)>3:                                                           # Check for knorm input, else request one
   try:
      knorm=str(args[3])
   except:
      knorm='nupnu'
else:
   try:
      knorm=str(raw_input('Input normalisation [leahy, rms, nupnu]: '))
   except:
      knorm='nupnu'

numstep=len(loadmatrix)

nleahy=float(sum(good))
lspec=npsum(loadmatrix,axis=0)/nleahy                                     # Create the average Leahy spectrum
const=pan.lhconst(lspec)                                                  # Calculate the normalisation of noise

def constmi(k):                                                           # Define nuP(nu) noise average as a function of  Leahy constant.
   nlspec=pan.lh2rms(lspec,mean(rates),mean(bg),k)
   nlrang=arange(len(nlspec))
   nlrang=nlrang[int(4*len(nlspec)/5.0):]
   nlspec=nlspec[int(4*len(nlspec)/5.0):]
   return (mean(nlrang*nlspec))

const=(brentq(constmi,const-0.1,const+0.1))                               # Minimise the above function to improve the constant

datres=int(foures/bsz)
tfl=linspace(0.0, (1.0/2.0)*datres/float(foures), (datres/2)+1)           # Create linearly spaced frequency domain up to the Nyquist frequency 1/2 (N/T)
tfl=tfl[:-1]
nulldat=zeros((datres/2))                                                 # Create null data with the same number of points as tfl
tf,null,null=pan.lbinify(tfl[1:],nulldat[1:],nulldat[1:],lplres)          # Fetch new array of bins to be output after lbinning
del null

print ''

def lbin(lplres,norm='nupnu',prt=False):                                  # Defining a log-binning function that just depends on bin resolution and normalisation

   if True:                                                               # Instructions with a list of possible normalisations
      if norm not in ['rms','nupnu','leahy']:
         print 'Unknown normalisation selected!'
         print ''
         print 'Available normalisations are:'
         print '* "leahy" for Leahy-normalised power'
         print '* "rms" for RMS-normalised power'
         print '* "nupnu" for RMS-normalised power multiplied by frequency'
         print ''
         print 'Using "leahy" normalisation:'
         norm='leahy'
      else:
         print 'Using "'+norm+'" normalisation:'

   errgr=[]                                                               # Set up matrix of errors
   fourgr=[]
   for i in range(numstep):
      tsfdata=loadmatrix[i]                                               # Load a row of data
      errs=pan.lh2rms(tsfdata,rates[i],bg[i],0)                           # Errors of a Leahy spectrum = the Leahy spectrum

      if norm in ['rms','nupnu']:
         if norm=='nupnu':
            sconst=const                                                  # For nuP(nu) normalisation, the Leahy constant will need subtracting
         else:
            sconst=0                                                      # In RMS norm, it can stay
         tsfdata=pan.lh2rms(tsfdata,rates[i],bg[i],sconst)                # Convert to RMS-normalised data using the LH2RMS function from pan_lib
         if norm=='nupnu':
            tsfdata=tsfdata*tfl                                           # Multiply by frequency if nupnu normalisation requested
            errs=errs*tfl
      tf,fours,errs=pan.lbinify(tfl[1:],tsfdata[1:],errs[1:],lplres)      # Logarithmically bin the data using lbinify from pan_lib
      fourgr.append(fours)                                                # Populate the data matrix
      errgr.append(abs(errs))                                             # Populate the error matrix

      prog=i+1
      if prt and ((prog % 5)==0 or prog==numstep):
         print str(prog)+'/'+str(numstep)+' series re-binned...'          # Display progress every 5 series

   fourgr=transpose(fourgr)                                               # Flip the matrices (makes them easier to plot the correct way round in spectrogram)
   errgr=transpose(errgr)
   return fourgr,errgr,norm

fourgr,errgr,knorm=lbin(lplres,prt=True,norm=knorm)  

print ''
print 'Preparing spectrogram...'

deftitle='Spectrogram'+qflav                                              # Define default title for spectrogram

if knorm=='leahy':                                                        # Define default key label for spectrogram
   defzlabl='Leahy-Normalised Power (Hz^-1)'
elif knorm=='rms':
   defzlabl='RMS Normalised Power (Hz^-1)'
else:
   defzlabl='Frequency x RMS Normalised Power'

fourgrm=fourgr                                                            # Storing a copy of the matrix in memory so it can be reset
errgrm=errgr
td=arange(0,(numstep+1)*slide,slide)                                      # Creating the time domain as an array
tdg, tfg = meshgrid(td, tf)                                               # Making a grid from the time and frequency domains

specopt=''                                                                # Force spectrogram manipulation mode to trigger
speclog=False                                                             # Indicate that the spectrogram is not initially logarithmic
stitle=deftitle                                                           # Give an initial title
rtlabl=defzlabl                                                           # Give an initial key label

tmdbin=foures                                                             # Initial time binning
frqbin=(tf[-1]-tf[0])/(len(tf)-1)                                         # Initial freq binning
tdgd=tdg                                                                  # Saving default grid [Time Domain Grid- Default]
tfgd=tfg
tdlm=td                                                                   # Saving 1D arrays to re-form grids [time-domain linear, modifiable]
tflm=tf
ogood=good                                                                # Save copy of the 'good' list

fudge=npmin(abs(fourgr[nonzero(fourgr)]))                                 # Obtain smallest nonzero value in array to add on when using logarithm to prevent log(0)

print 'Done!'
print ''


#-----Setting up Spectrogram Environment-------------------------------------------------------------------------------

es=True                                                                   # Start with errors on by default

def spectrogram(td,tfc,fourgr,zlabel=defzlabl,title=deftitle):            # Defining the creation of the spectrogram plot 's' as a function for clarity later
   pl.close('Spectrogram')                                                # Close any previous spectrograms that may be open
   fg=pl.figure('Spectrogram')
   ax=fg.add_subplot(1,1,1)
   if speclog:
      pl.pcolor(td,tfc,(fourgr))
   else:
      pl.pcolor(td,tfc,fourgr,vmin=sgfloor,vmax=sgceil)                   # Plot spectrogram
   cbar=pl.colorbar()                                                     # Create colourbar key
   cbar.set_label(zlabel)
   pl.xlabel('Time(s)')
   pl.ylabel('Frequency(Hz)')
   pl.title(title)
   ax.set_yscale('log')                                                   # Make the freq-axis logarithmic too
   pl.ylim(tfc[0,0],tfc[-1,-1])                                           # Resize axes to fit the spectrogram
   pl.xlim(td[0,0],td[-1,-1])
   pl.show(block=False)

sxlab='Frequency (Hz)'
sylab=defzlabl
szlab=defzlabl                                                            # Give an initial key label, storing second copy

def give_inst():                                                          # Define printing this list of instructions as a function
   print 'COMMANDS: Enter a command to manipulate data.'
   print ''
   print 'DATA:'
   print '* "rebin" to reset the data and load it with a different normalisation and binning.'
   print '* "clip" to clip the range of data.'
   print '* "reset" to reset data.'
   print ''
   print 'SPECTROGRAM:'
   print '* "sg plot" to plot the spectrogram currently being worked on.'
   print '* "sg floor" to set a minimum value for the spectrogram'+"'"+'s z-axis colour key.'
   print '* "sg ceil" to set a maximum value for the spectrogram'+"'"+'s z-axis colour key.'
   print '* "sg auto" to automatically set colour floor and ceiling.'
   print '* "sg log" to toggle logarithmic spectrogram plotting.'
   print ''
   print 'POWER SPECTRA:'
   print '* "aspec" to plot the average spectrum and return the frequency of its highest peak.'
   print '* "gspec" to get an individual spectrum at any time and plot it.'
   print '* "peaks" to plot a graph of the frequency of the strongest oscillation against time.'
   print '* "rates" to get a simple lightcurve of the data.'
   print '* "fqflux" to plot "peaks" against "rates".'
   print ''
   print 'TOGGLE OPTIONS:'
   print '* "errors" to toggle errorbars on power spectra plots.'
   print ''
   print 'OTHER COMMANDS:'
   print '* "info" to display a list of facts and figures about the current SpecAngel session.'
   print '* "reflav" to rewrite the flavour text used for graph titles.'
   print '* "leahy dump" to create an ASCII file of the Leahy-normalised average power density spectrum.'
   print '* "help" or "?" to display this list of instructions again.'
   print '* "quit" to Quit'

give_inst()                                                               # Print the list of instructions
print ''
print ' --------------------'

sgfloor=max(npmin(fourgrm),0)
sgceil=npmax(fourgrm)


#-----Entering Interactive Mode----------------------------------------------------------------------------------------

while specopt not in ['quit','exit']:                                     # If the previous command given was not quit, continue

   print ''
   specopt=raw_input('Give command [? for help]: ')                       # Fetch command from user
   print ''


   #-----'rebin' Option------------------------------------------------------------------------------------------------

   if specopt=='rebin':                                                   # Rebinning data    

      print 'Data currently binned in '+str(tmdbin)+'s and 10^'+str(lplres)+'Hz bins.'

      speclog=False                                                       # Indicate that the spectrogram is not logarithmic
      stitle=deftitle                                                     # Restore initial title
      tmdbin=foures                                                       # Initial time binning
      frqbin=(tf[-1]-tf[1])/(len(tf)-2)                                   # Initial freq binning
      fourgrm=fourgr                                                      # Reload original, unmodified data
      errgrm=errgr
      tdgd=tdg                                                            # Reloading default grid
      tfgd=tfg
      tdlm=td                                                             # Loading 1D arrays to re-form grids [time-domain linear, modifiable]
      tflm=tf
      good=ogood

      knorm=raw_input('Select normalisation [leahy, rms, nupnu]: ')

      try:
         tbinmult=int(raw_input('Input time-domain binning factor: '))
         if tbinmult<1: tbinmult=1                                        # Prevent bin-sizes smaller than current bin-size
         if tbinmult>len(fourgrm)/2.0:tbinmult=1                          # Forces there to be at least 2 bins

      except:
         tbinmult=1                                                       # Cancel binning if a non-number is entered
         print 'Invalid time binning!'

      tmdbin=foures*tbinmult                                              # Recover current binning

      try:
         newfbin=float(raw_input('Input new freq-domain exponent: '))
         if newfbin>0: lplres=newfbin                                     # Prevent bin-sizes smaller than zero
             
      except:
         newfbin=0                                                        # Cancel binning if a non-number is entered
         print 'Invalid frequency binning!'

      print ''
      print 'Re-binning...'

      tflm,null,null=pan.lbinify(tfl[1:],nulldat,nulldat,lplres)          # Fetch new array of bins to be output after lbinning
      del null
      fourgrm,errgrm,knorm=lbin(lplres,prt=True,norm=knorm)               # Re log-bin data

      if tbinmult!=1:                                                     # Cancel binning if new bin is not greater than old bin

         fourgrm,errgrm,tdlm,good=pan.mxrebin(fourgrm,errgrm,tdlm,good,tbinmult)

      tdgd,tfgd=meshgrid(tdlm,tflm)                                       # Recreate grid from rescaled axes

      print ''
      print 'Data rebinned by '+str(tmdbin)+'s, [10^'+str(lplres)+'n]Hz.'
      print str(int(sum(good)))+'/'+str(len(good))+' power spectra are good'

      if knorm=='leahy':
         sylab='Leahy-Normalised Power (Hz^-1)'
      elif knorm=='rms':
         sylab='RMS Normalised Power (Hz^-1)'
      else:
         sylab='Frequency x RMS Normalised Power'

      defzlabl=sylab                                                      # Restore root z label for spectrogram              
      rtlabl=defzlabl                                                     # Restore initial key label
      szlab=defzlabl                                                      # Restore initial key label, storing second copy

      sgfloor=max(npmin(fourgrm),0)                                       # Reset spectrogram colour floor & ceil
      sgceil=npmax(fourgrm)


   #-----'sg plot' Option----------------------------------------------------------------------------------------------

   # 'Spectrogram'

   elif specopt=='sg plot':                                               # Plotting data

      proce=True                                                          # Assume plot will be made
      npl=len(fourgrm[0,:])*len(fourgrm[:,0])                             # Check how large this plot will be

      if npl>1000000:                                                     # Check for very large plots, ask user whether to proceed
          try:
             print "LARGE DATA WARNING! Plot will contain "+"{:,}".format(npl)+" elements."
             proc=raw_input("Proceed? [y/n]: ")
             if proc!='y': proce=False                                    # Cancel plot is user doesn't explicity say 'y'
          except:
             proce=False

      if proce==True:                                                     # If all is ok...
         print 'Plotting...'
         spectrogram(tdgd,tfgd,fourgrm,szlab,stitle)                      # Plot spectrogram 


   #-----'sg floor' Option---------------------------------------------------------------------------------------------

   # 'Colour Floor'

   elif specopt=='sg floor':                                              # Setting floor of spectrogram colour scale:

      sgfloor=raw_input('Input spectrogram colour floor: ')               # Ask user to input floor
      try:
         sgfloor=float(sgfloor)                                           # Check floor is a number
         print 'Colour floor set!'
      except:
         sgfloor=max(npmin(fourgrm),0)
         print 'Invalid floor!'


   #-----'sg ceil' Option----------------------------------------------------------------------------------------------

   # 'Colour Ceiling'

   elif specopt=='sg ceil':                                               # Setting ceiling of spectrogram colour scale:

      sgceil=raw_input('Input spectrogram colour ceiling: ')              # Ask user to input floor
      try:
         sgceil=float(sgceil)                                             # Check ceil is a number
         print 'Colour ceiling set!'
      except:
         sgceil=npmax(fourgrm)
         print 'Invalid ceiling!'


   #-----'sg auto' Option----------------------------------------------------------------------------------------------

   # 'Auto Recolour'

   elif specopt=='sg auto':

      sgfloor=max(npmin(fourgrm),0)                                       # Reset spectrogram colour floor & ceil
      sgceil=npmax(fourgrm)

      print 'Colour floor and ceiling automatically set!'


   #-----'sg log' Option-----------------------------------------------------------------------------------------------

   # 'Logarithm'

   elif specopt=='sg log':                                                # Taking or undoing log of data

      if speclog:                                                         # If the spectrogram was already logged, undo this with exp

         print 'Exponentiating spectrogram...'

         stitle=deftitle                                                  # Reset title to default
         rtlabl=rtlabl[4:]                                                # Remove 'log ' from the start of both saved z-labels
         szlab=szlab[4:]
         fourgrm=10**(fourgrm)-fudge                                      # Exponentiate every element of every power spectrum
         speclog=False                                                    # Indicate that spectrogram is no longer logged
         print 'Done!'

      else:

         print 'Taking logarithm of spectrogram...'

         stitle='Log '+deftitle                                           # Add 'log ' to start of titles and labels
         rtlabl='Log '+rtlabl
         szlab='Log '+szlab
         fourgrm=log10(abs(fourgrm)+fudge)                                # Take the log of every element of every power spectrum
         speclog=True                                                     # Indicate that spectrum is logged
         print 'Done!'


   #-----'sg' Catch-All Help Message-----------------------------------------------------------------------------------

   elif specopt=='sg':

      print 'SPECTROGRAM COMMANDS:'
      print '* "sg plot" to plot the spectrogram currently being worked on.'
      print '* "sg floor" to set a minimum value for the spectrogram'+"'"+'s z-axis colour key.'
      print '* "sg ceil" to set a maximum value for the spectrogram'+"'"+'s z-axis colour key.'
      print '* "sg auto" to automatically set colour floor and ceiling.'
      print '* "sg log" to toggle logarithmic spectrogram plotting.'


   #-----'clip' Option-------------------------------------------------------------------------------------------------

   elif specopt=='clip':                                                  # Clipping data

      print 'Clipping data'
      print ''

      print 'Time range is '+str(tdlm[0])+'s - '+str(tdlm[-2]+foures)+'s'
      print 'Freq range is '+str(tflm[0])+'Hz- '+str(tflm[-1]+foures)+'Hz'
      print ''

      print 'Please choose new range of data:'
      mint,maxt,null=pan.srinr(tdlm,tmdbin,'time')                        # Fetch new time domain endpoints using srinr function from pan_lib
      minf,maxf,null=pan.srinr(tflm,frqbin,'freq')                        # Fetch new freq domain endpoints using srinr function from pan_lib

      print 'Clipping...'

      tdlm=tdlm[mint:maxt]                                                # Clip the time-domain array
      tflm=tflm[minf:maxf]                                                # Clip the freq-domain array
      fourgrm=fourgrm[minf:maxf,mint:maxt]                                # Clip the spectrogram data
      errgrm=errgrm[minf:maxf,mint:maxt]
      tdgd,tfgd=meshgrid(tdlm,tflm)                                       # Recreate grid from rescaled axes

      print 'Data clipped!'


   #-----'reset' Option------------------------------------------------------------------------------------------------

   elif specopt=='reset':                                                 # Resetting data

      print 'Resetting spectrogram...'
               
      speclog=False                                                       # Indicate that the spectrogram is not logarithmic
      stitle=deftitle                                                     # Restore initial title
      rtlabl=defzlabl                                                     # Restore initial key label
      szlab=defzlabl                                                      # Restore initial key label, storing second copy
      tmdbin=foures                                                       # Initial time binning
      frqbin=(tf[-1]-tf[1])/(len(tf)-2)                                   # Initial freq binning
      fourgrm=fourgr                                                      # Reload original, unmodified data
      errgrm=errgr
      tdgd=tdg                                                            # Reloading default grid
      tfgd=tfg
      tdlm=td                                                             # Loading 1D arrays to re-form grids [time-domain linear, modifiable]
      tflm=tf
      good=ogood                                                          # Reload original 'good' list

      print 'Spectrogram reset!'


   #-----'aspec' Option------------------------------------------------------------------------------------------------

   # 'Average Spec'

   elif specopt=='aspec':                                                 # Find the time-averaged spectrum

      if slide!=foures:
         print "Warning!  Data taken with sliding window: data points not independent!"
         print ''

      print "Fetching time-averaged power spectrum..."

      spec=npsum(fourgrm, axis=1)/sum(good)                               # Sum all spectra in the matrix and divide by the number of good columns
      err=sqrt(npsum( array(errgrm)**2, axis=1))/sum(good)
      ttl='Average power density spectrum'+qflav
      pan.slplot(tflm,spec,err,sxlab,sylab,ttl,'spc',typ='log',errors=es) # SLPlot from the pan_lib plots data on standard and log-log axes
      scerr=spec-(err**0.5)

      print 'Maximum power found at '+str(tflm[scerr.argmax()])+'Hz!'     # Suggest a peak location
      print '  (Period of '+str(1.0/tflm[scerr.argmax()])+'s)'
      

   #-----'gspec' Option------------------------------------------------------------------------------------------------

   # 'Get Spec'

   elif specopt=='gspec':                                                 # Find the power spectrum at a specific point in time

      if slide!=foures:
         print "Warning!  Data taken with sliding window: data points not independent!"
         print ''

      print 'Getting a spectrum at a time (since start of observation) of your choice.'
      print ''

      while True:
         try:
            specid=float(raw_input('Enter time: '))                       # Fetch raw time suggestion from user
            break
         except:
            print 'Invalid time!'


      specid=int((specid-tdlm[0])/tmdbin)                                 # Work out which time bin this would correlate to

      if 0<=specid<len(fourgrm[0,:])-1:                                   # Check that this time bin actually exists in the matrix
         if good[specid]:
            print "Fetching power spectrum at "+str(specid*tmdbin)+"s: "

            gsp=fourgrm[:,specid]                                         # Extract the lightcurve from this bin
            ger=errgrm[:,specid]
            ttl='Power density spectrum "'+flavour+'" at +'+str(specid*tmdbin)+'s'
            pan.slplot(tflm,gsp,ger,sxlab,sylab,ttl,'spc',typ='log',errors=es)
            scerr=spec-(err**0.5)
            
            print 'Maximum power found at '+str(tflm[scerr.argmax()])+'Hz!' # Suggest a peak location
             
         else:
            print 'Time not in GTI!'
      else:
         print 'Time not in range!'


   #-----'peaks' Option------------------------------------------------------------------------------------------------

   elif specopt=='peaks':

      peaks=[]

      for i in range(len(fourgrm[1])):
         row=fourgrm[:,i]
         peaks.append(tflm[list(row).index(max(row))])

      pl.close('pk')
      pl.figure('pk')
      pl.semilogy(array(tdlm[:-1])[good],array(peaks)[good],'-ok')
      pl.xlabel('Time (s)')
      pl.ylabel('Frequency (Hz)')
      pl.title('Peak Frequency/Time Plot'+qflav)
      pl.show(block=False)


   #-----'rates' Option------------------------------------------------------------------------------------------------

   elif specopt=='rates':

      print 'Average rate of',str(mean(rates[ogood]))+'c/s.'
      print str(phcts),'total counts.'

      print ''

      datasel=raw_input('Select Rates to plot [ave/peak/trough]: ')

      titles['ave']='Flux (photons/s/PCU)'
      titles['peak']='Peak '+str(binfac*bsz)+'s Flux (photons/s/PCU)'
      titles['trough']='Trough '+str(binfac*bsz)+'s Flux (photons/s/PCU)'
      if datasel not in ('ave','peak','trough'):
         print 'Invalid selection!  Using Average flux.'
         datasel='ave'

      brates['ave']=rates
      brates['peak']=prates
      brates['trough']=trates

      pl.close('lc')
      pl.figure('lc')
      pl.plot(td[:-1][ogood],brates[datasel][ogood],'-ok')
      pl.xlabel('Time (s)')
      pl.ylabel('Flux (photons/s/PCU)')
      pl.title('Lightcurve'+qflav)
      pl.show(block=False)


   #-----'fqflux' Option-----------------------------------------------------------------------------------------------

   elif specopt=='fqflux':

      peaks=[]

      for i in range(len(fourgrm[1])):
         row=fourgrm[:,i]
         peaks.append(tflm[list(row).index(max(row))])

      brates={}
      titles={}

      titles['ave']='Flux (photons/s/PCU)'
      titles['peak']='Peak '+str(binfac*bsz)+'s Flux (photons/s/PCU)'
      titles['trough']='Trough '+str(binfac*bsz)+'s Flux (photons/s/PCU)'

      brates['ave']=pan.vcrebin(rates,len(rates)/len(peaks))
      brates['peak']=pan.vcrebin(prates,len(prates)/len(peaks))
      brates['trough']=pan.vcrebin(trates,len(trates)/len(peaks))

      datasel=raw_input('Select Rates to plot against [ave/peak/trough]: ')
      if datasel not in ('ave','peak','trough'):
         print 'Invalid selection!  Using Average flux.'
         datasel='ave'

      print 'Spearman Rank Coefficient: ',spearmanr(array(peaks)[good],brates[datasel][good])[1]

      pl.close('pr')
      pl.figure('pr')
      pl.semilogx(array(peaks)[good],brates[datasel][good],'ok')
      pl.ylabel(titles[datasel])
      pl.xlabel('Frequency (Hz)')
      pl.title('Flux/Peak Frequency Plot'+qflav)
      pl.show(block=False)


   #-----'errors' Option-----------------------------------------------------------------------------------------------

   elif specopt in ['error', 'errors']:                                   # Toggle Errors

      if es:
         es=False
         print 'Errors suppressed!'
      else:
         es=True
         print 'Errors displayed!'


   #-----'info' Option-------------------------------------------------------------------------------------------------

   elif specopt=='info':

      print 'SpecAngel.py version',version
      print ''
      print '1 file loaded:'
      print ''
      filn,loca=pan.xtrfilloc(filename)
      print 'File 1:'
      print ' Filename       = ',filn
      print ' Location       = ',loca
      print ' Mission        = ',mis
      print ' Object         = ',obsd[0]
      print ' Obs_ID         = ',obsd[1]
      if mis in ['SUZAKU']:
         print ' Energy         = ',cs,'eV'
      else:
         print ' Channel        = ',cs
      print ' Resolution     = ',str(bsz)+'s'
      print ' Flavour        = ',flv
      print ' FITSGenie Ver. = ',v
      print ''
      print 'Windowing:'
      print ' Shape          = ',wtype
      print ' Sliding        = ',slide!=foures
      print ' Length         = ',str(foures)+'s'
      if slide!=foures:
         print ' Separation     = ',str(slide)+'s'
      print ''
      print 'Normalisation:'
      print ' Normalisation  = ',knorm
      print ' Leahy constant = ',const
      print ''
      print 'Other Info:'
      print ' Main Flavour   = ',flavour
      print ' Obs length     = ',str(foures*numstep)+'s'
      print ' Time. Bin-size = ',str(tmdbin)+'s'
      print ' Freq. Bin-size = ',lplres
      print ' Num. Time Bins = ',str(int(len(good)))
      print ' Good Time Bins = ',str(int(sum(good)))
      print ' Avg. Rates     = ',str(mean(rates[ogood]))
      print ' Total photons  = ',phcts
      print ' Background     = ',str(bgest)+'cts/s/PCU'
      print ' Errorbars      = ',es


   #-----'reflav' Option-----------------------------------------------------------------------------------------------

   elif specopt=='reflav':

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


   #-----'leahy dump' Option-------------------------------------------------------------------------------------------

   elif specopt=='leahy dump':

      aflname=raw_input('Filename:')
      try:

         assert len(aflname)>0
         fle=open(aflname,'w')
         for i in range(len(tfl)):
            a=[str(tfl[i]),' ',str(lspec[i]),' ',str(lspec[i]/sqrt(nleahy)),'\n']
            fle.writelines(a)
         fle.close()

         print 'ASCII Leahy-normalised spectrum saved to',aflname+'!'

      except:

         print 'Invalid filename!'


   #-----'help' Option-------------------------------------------------------------------------------------------------

   elif specopt in ['help','?']:                                          # Display instructions

      print 'Instructions:'
      print ''

      give_inst()                                                         # Re-call the instructions list, defined as the get_inst() function in initialisation


   #-----'quit' Option-------------------------------------------------------------------------------------------------

   elif specopt not in ['quit','exit']:                                   # Invalid command if none of the if statements triggered and no 'q' given

      print 'Invalid command!'

   if specopt not in ['quit','exit']:
      print ''
      print ' --------------------'


#-----Exiting Interactive Mode-----------------------------------------------------------------------------------------

print ''
print 'Goodbye!'                                           


#-----Footer-----------------------------------------------------------------------------------------------------------

pan.signoff()


