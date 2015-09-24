#! /usr/bin/env python

# |----------------------------------------------------------------------|
# |-----------------------------Back HYDRA-------------------------------|
# |----------------------------------------------------------------------|

# Call as ./bckghydra.py DATA_FILE BACK_FILE SAVE_FILE

# Takes 2 .plotd infiles, one of data and one of the background of the same observation,
# and returns the first file after background subtraction.
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
   import numpy   as np
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
backfile_packed=pan.plotdld(back_filename)                                # Load background file


#-----Unpack Data------------------------------------------------------------------------------------------------------

# Collect all data from the data file

data_x        = datafile_packed[0]
data_f        = datafile_packed[1]
data_fe       = datafile_packed[2]
data_t_start  = datafile_packed[3]
data_bin_size = datafile_packed[4]
data_gti      = datafile_packed[5]
data_maxpcus  = datafile_packed[6]
data_bg_est   = datafile_packed[7]
data_flavour  = datafile_packed[8]
data_channels = datafile_packed[9]
data_mission  = datafile_packed[10]
data_obs_data = datafile_packed[11]
data_fitsg_v  = datafile_packed[12]

data_obsid=data_obs_data[1]

# Collect only relevant data from the background file

back_x        = backfile_packed[0]
back_f        = backfile_packed[1]
back_fe       = backfile_packed[2]
back_t_start  = backfile_packed[3]
back_bin_size = backfile_packed[4]
back_maxpcus  = backfile_packed[6]
back_channels = backfile_packed[9]
back_mission  = backfile_packed[10]
back_obs_data = backfile_packed[11]

back_obsid=back_obs_data[1]


#-----Check Background and Data files are compatible-------------------------------------------------------------------

same_channels = (data_channels == back_channels)                          # Check the channel ranges match
same_mission  = ( data_mission == back_mission )                          # Check the mission names match
same_obsid    = (   data_obsid == back_obsid   )                          # Check the OBSIDs match

if not same_mission:                                                      # Abort if missions differ
   print 'Files are from different missions!'
   print 'Aborting!'
   pan.signoff()
   exit()

if not same_obsid:                                                        # Abort if observations differ
   print 'Files are from different observations!'
   print 'Aborting!'
   pan.signoff()
   exit()

if not same_channels:                                                     # Abort if channel ranges differ
   print 'Files contain data from different energy ranges!'
   print 'Aborting!'
   pan.signoff()
   exit()
 
#-----Shift Arrays-----------------------------------------------------------------------------------------------------
  
shifted_data_x=data_x+(back_t_start-(data_t_start+back_x[0]))             # Created a 'data_x array' with its startpoint aligned with background
back_x=pan.tnorm(back_x,back_bin_size)                                    # Force background x array to start at 0
back_f=back_f*(back_bin_size/data_bin_size)                               # Rescale background counts to reflect the imminent rebinning


#-----Define Background Subtraction------------------------------------------------------------------------------------

def backgr(i):                                                            # Function that returns the appropriate background counts at each point in the datafile
   if shifted_data_x[i]<back_x[0]:
      return back_f[0],back_fe[0]                                         # Return startpoint background if sampling before bg range
   elif shifted_data_x[i]>=back_x[-1]:
      return back_f[-1],back_fe[0]                                        # Return endpoint background if sampling after bg range
   else:
      timestamp=shifted_data_x[i]                                         # Collect the timestamp of the ith data element
      st_i=int(timestamp/back_bin_size)                                   # Collect the start and endpoints of the bg bin in which the timestamp falls
      ed_i=st_i+1
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

pan.plotdsv(save_filename,data_x,data_f,data_fe,data_t_start,data_bin_size,
            data_gti,data_maxpcus,data_bg_est,data_flavour,data_channels,
            data_mission,data_obs_data,data_fitsg_v)

print ''
print 'Background Subtracted file saved as "'+save_filename+'"!'


