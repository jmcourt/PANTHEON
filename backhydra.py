#! /usr/bin/env python

# |----------------------------------------------------------------------|
# |-----------------------------BACK HYDRA-------------------------------|
# |----------------------------------------------------------------------|

# Call as ./bckghydra.py DATA_FILE BACK_FILE SAVE_FILE

# Takes a .plotd file and a background file created with PCABACKEST and returns
#
# Arguments:
#
#  DATA_FILE
#   The absolute path to the file to be used as data.
#
#  BACK_FILE
#   The file to be used as background; does not need to be the same binning as File 1.
#   suggest using pcabackest from FTOOLS to produce this file.
#   FTOOLS can be found at http://heasarc.gsfc.nasa.gov/ftools/
#
#  SAVE_FILE
#   The location to save the resultant background-subtracted file
#

#-----Welcoming Header-------------------------------------------------------------------------------------------------

print ''
print '-------Running BackHydra: J.M.Court, 2015-------'
print ''


#-----Importing Modules------------------------------------------------------------------------------------------------

try:

   import sys
   import pan_lib as pan
   from astropy.io import fits
   import numpy as np
   import pylab as pl

except ImportError:

   print 'Modules missing!  Aborting!'
   print ''
   print '------------------------------------------------'
   print ''
   exit()


#-----Checking Validity of Filenames-----------------------------------------------------------------------------------

args=sys.argv
pan.argcheck(args,4)                                                      # Must give at least 3 args (Both filenames function call)

data_filename=args[1]                                                     # Fetch datafile name from arguments
back_filename=args[2]                                                     # Fetch background file name from arguments
save_filename=args[3]

print 'Loading Data...'
datafile_packed=pan.plotdld(data_filename)                                # Load datafile
print 'Loading Background...'
backfile_packed=fits.open(back_filename)


#-----Unpack Data------------------------------------------------------------------------------------------------------

# Collect all data from the data file

b_sub='True'
data_x        = datafile_packed[0]
data_f        = datafile_packed[1]
data_fe       = datafile_packed[2]
data_t_start  = datafile_packed[3]
data_bin_size = datafile_packed[4]
data_gti      = datafile_packed[5]
data_maxpcus  = datafile_packed[6]
data_flavour  = datafile_packed[10]
data_channels = datafile_packed[11]
data_mission  = datafile_packed[12]
data_obs_data = datafile_packed[13]
data_fitsg_v  = datafile_packed[14]

data_obsid=data_obs_data[1]

# Collect only relevant data from the background file

backfile_data=backfile_packed[1].data
back_mission=backfile_packed[1].header['TELESCOP']

if back_mission == 'XTE' :
   try:
      import xtepan_lib as inst                                           # Import XTE extraction functions
   except:
      print 'XTE PANTHEON Library not found!  Aborting!'
      pan.signoff()
      exit()

elif back_mission == 'SUZAKU':
   try:
      import szkpan_lib as inst                                           # Import SUZAKU extraction functions
   except:
      print 'Suzaku PANTHEON Library not found!  Aborting!'
      pan.signoff()
      exit()

low_chan,high_chan=data_channels.split('-')

back_x,back_f,back_fe = inst.getbg(backfile_data,int(low_chan),int(high_chan))
back_x=back_x[back_f>0]
back_fe=back_fe[back_f>0]
back_f=back_f[back_f>0]
back_t_start  = back_x[0]
back_bin_size = inst.getbin(backfile_packed,None)


#-----Check Background and Data files are compatible-------------------------------------------------------------------

same_mission  = ( data_mission == back_mission )                          # Check the mission names match

if not same_mission:                                                      # Abort if missions differ
   print 'Files are from different missions!'
   print 'Aborting!'
   pan.signoff()
   exit()

 
#-----Shift Arrays-----------------------------------------------------------------------------------------------------

back_x=pan.tnorm(back_x,back_bin_size)                                    # Force background x array to start at 0

shifted_data_x=data_x+data_t_start-back_t_start                           # Create shifted data axis to align with a background starting at 0s
shifted_back_x=back_x+back_t_start-data_t_start                           # Create shifted background axis to align with data starting at 0s

if back_x[0]>shifted_data_x[-1] or shifted_data_x[0]>back_x[-1]:          # Abort if the timescales don't overlap
   print 'WARNING! Files times do not overlap!'
   print ''
   print 'Estimating constant background.'
   print ''
   b_sub='Estimate'
   dump_file=open('backhydra_log.txt','w')
   dump_file.write('File and background times did not overlap!')
   dump_file.close()


#-----Define Background Subtraction------------------------------------------------------------------------------------

def backgr(i):                                                            # Function that returns the appropriate background counts at each point in the datafile

   timestamp=shifted_data_x[i]                                            # Collect the timestamp of the ith data element
   st_i=int(timestamp/back_bin_size)                                      # Collect the start and endpoints of the bg bin in which the timestamp falls
   ed_i=st_i+1

   if st_i<0:
      return back_f[0],back_fe[0]                                         # Return startpoint background if sampling before bg range
   elif ed_i>=len(back_x):
      return back_f[-1],back_fe[-1]                                       # Return endpoint background if sampling after bg range
   else:

      posit_in_bin=(timestamp % back_bin_size)/back_bin_size              # Work out where in the bin the timestamp falls
      f_est  = back_f[st_i]+posit_in_bin*(back_f[ed_i]-back_f[st_i])      # Linearly interpolate between two background points to return background estimate
      fe_est = back_fe[st_i]+posit_in_bin*(back_fe[ed_i]-back_fe[st_i])   # Collect error too

      return f_est,fe_est


#-----Perform Background Subtraction-----------------------------------------------------------------------------------

print 'Subtracting Background...'

for i in pan.eqrange(data_x):
   back,back_e=backgr(i)
   data_f[i]-=back
   data_fe[i]=(data_fe[i]**2+back_e**2)**0.5


#-----Re-save Data-----------------------------------------------------------------------------------------------------

print 'Saving...'

new_bg_est=np.mean(back_f)/data_maxpcus
bg_data=(shifted_back_x[(shifted_back_x>=data_x[0]) | (shifted_back_x<=data_x[-1])],back_f[(shifted_back_x>=data_x[0]) | (shifted_back_x<=data_x[-1])])

pan.plotdsv(save_filename,data_x,data_f,data_fe,data_t_start,data_bin_size,
            data_gti,data_maxpcus,new_bg_est,b_sub,bg_data,data_flavour,
            data_channels,data_mission,data_obs_data,data_fitsg_v)

print ''
print 'Background Subtracted file saved as "'+save_filename+'.plotd"!'


