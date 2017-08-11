#! /usr/bin/env python

# |----------------------------------------------------------------------|
# |-----------------------------SWFPAN_LIB-------------------------------|
# |----------------------------------------------------------------------|

# A selection of useful Swift-specific functions which are placed here to reduce clutter in the other
# files of PANTHEON.  Supports event data extracted with xselect.
#
# Contents:
#
#  CHRANGE   - clips a set of event data given to it to screen out photons outside of some energy range
#              given by the user.
#
#  DISCNEV   - Discards non-photon events from a table of data
#
#  GETBIN    - gets the binning time of a FITS data table
#
#  GETDAT    - gets the data of a FITS data table
#
#  GETGTI    - gets the GTIs associated with a FITS data table
#
#  GETINI    - gets the start time of the observation associated with a FITS data table
#
#  GETPCU    - gets the number of PCUs that contributed events to a FITS data table
#
#  MAXEN     - returns the highest energy or channel valid for the instrument.
#

# Modes supported:

# Event data extracted with xselect


#-----Importing Modules------------------------------------------------------------------------------------------------

import pan_lib as pan
from numpy import array, zeros
from numpy import sum as npsum


#-----ChRange----------------------------------------------------------------------------------------------------------

def chrange(data,low,high,header):

   '''Channel Ranger

   Description:

    Takes xselected event data from Swift and removes all photons outside of an energy
    channel range given by the user.

   Inputs:

    data     - FITS DATA: the event data; raw photon arrival times and data words.
    low      -       INT: the lowest channel the user wants to consider.
    high     -       INT: the highest channel the user wants to consider.
    datamode -    STRING: the datamode in which the data was taken.

   Outputs:

    ch_data  - FITS DATA: the same as data, but with all photons outside of the given channel range
                          removed.

   -J.M.Court, 2015'''

   datamode=header['DATAMODE']
   energies=data.field(1)

   if low<=0 and high>=4096:
      return data                                                      # Don't bother searching through if the user wants full range


   low=evmchan(low)                                                 # Convert the channels given into range IDs
   high=evmchan(high)

   energies=array(energies)

   mask1=(energies>=low)

   mask2=(energies<=high)

   ch_data=data[mask1&mask2]

   return ch_data


#-----DiscNEv----------------------------------------------------------------------------------------------------------

def discnev(datas,datamode):

   '''Discard Non-Events

   Decription:
    Given a FITS data table, discards all non-photon events from the table.

   -J.M.Court, 2015'''

   return datas


#-----GetBin-----------------------------------------------------------------------------------------------------------

def getbin(event,datamode):

   '''Get Bin

   Description:

    Returns the time binning of a given set of FITS data.

   Inputs:

    event - FITS OBJECT: The FITS file that has been opened.

   Outputs:

    bsz   -       FLOAT: The bin size of the observation (s).

   -J.M.Court, 2015'''

   bsz=event[1].header['TIMEDEL']
   return bsz


#-----GetDat-----------------------------------------------------------------------------------------------------------

def getdat(event):

   '''Get data

   Description:

    Returns the data table of a given set of FITS data.

   Inputs:

    event - FITS OBJECT: The FITS file that has been opened.

   Outputs:

    data  -  FITS TABLE: The data table from the FITS file.

   -J.M.Court, 2015'''

   data=event[1].data                                                     # Extract GTI indices
   return data


#-----GetGTI-----------------------------------------------------------------------------------------------------------

def getgti(event):

   '''Get GTI

   Description:

    Returns the GTI of a given set of FITS data.

   Inputs:

    event - FITS OBJECT: The FITS file that has been opened.

   Outputs:

    gti   -     2D LIST: A list, each element of which contains the start and end point of a GTI.

   -J.M.Court, 2015'''

   gti=event[2].data                                                      # Extract GTI indices
   return gti


#-----GetIni-----------------------------------------------------------------------------------------------------------

def getini(event):

   '''Get Ini

   Description:

    Returns the initial (lowest) time of a given set of FITS data.

   Inputs:

    event - FITS OBJECT: The FITS file that has been opened.

   Outputs:

    ini   -       FLOAT: The starting time of the observation (s).

   -J.M.Court, 2015'''

   ini=event[1].header['TSTART']
   return ini


#-----Get Obs----------------------------------------------------------------------------------------------------------

def getobs(event,datamode,filepath):

   '''Get Obs

   Description:

    Fetches a tuple consisting of the object and obs_id of the observation.'''

   obsdata=(event[1].header['OBJECT'],event[1].header['OBS_ID'])

   return obsdata


#-----Get PCU----------------------------------------------------------------------------------------------------------

def getpcu(words,header):

   '''Get PCUs

   Description:

    Not valid for Swift.  Returns 1.

   Inputs:

    words    - 2D ARRAY: The array of event words, each of which is an array of True/False statements.
    datamode -   STRING: The DATAMODE in which the data to be analysed is stored.

   Outputs:

    pcus     -      INT: The number of PCUs active when the data was taken.  If a non-recognised event
                         word is given, 1 is returned instead along with an error message.

   -J.M.Court, 2015'''

   return 1


#-----Get Tim----------------------------------------------------------------------------------------------------------

def gettim(data,tstart,res,datamode):

   '''Get Times

   Description: 
   
    Returns observation start time

   -J.M.Court, 2015'''

   return data.field(0) 


#-----Get Wrd----------------------------------------------------------------------------------------------------------

def getwrd(data,datamode):

   '''Get Words

   Description: 

    Not valid for Swift.  Returns None.

   -J.M.Court, 2015'''

   return None


#-----Get Wrd Row------------------------------------------------------------------------------------------------------

def getwrdrow(words,mask,datamode):

   '''Get Word Rows

    Not valid for Swift.  Returns None.

   -J.M.Court, 2015'''

   return None


#-----MaxEn------------------------------------------------------------------------------------------------------------

def maxen(datamode):

   '''Max Energy

   Description:
    Returns the highest energy or channel valid for the instrument.

   -J.M.Court, 2015'''

   return 4096


