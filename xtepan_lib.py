#!/usr/bin/python

# |----------------------------------------------------------------------|
# |-----------------------------XTEPAN_LIB-------------------------------|
# |----------------------------------------------------------------------|

# A selection of useful XTE-specific functions which are placed here to reduce clutter in the other
# files of PANTHEON.  Supports GoodXenon_2s and E_125us_64M_0_1s DATAMODEs.
#
# Contents:
#
#  CHRANGE   - clips a set of event data given to it to screen out photons outside of some energy range
#              given by the user.
#
#  DISCNEV   - Discards non-photon events from a table of data
#
#  EVMCHAN   - converts an RXTE channel ID into an E_125us_64M_0_1s DATAMODE channel range ID.
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

#-----Importing Modules------------------------------------------------------------------------------------------------

import pan_lib as pan
from numpy import array


#-----ChRange----------------------------------------------------------------------------------------------------------

def chrange(data,low,high,datamode):

   '''Channel Ranger

   Description:

    Takes raw event or GoodXenon data from PCA on RXTE and removes all photons outside of an energy
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

   words=data.field(1)

   if low<=0 and high>=255:
      return data                                                         # Don't bother searching through if the user wants full range

   if datamode=='E_125us_64M_0_1s':
      low=evmchan(low)                                                    # Convert the channels given into range IDs
      high=evmchan(high)
      r=5,11                                                              # Identify where in the E_125 data word the channel is hidden
   elif datamode=='GoodXenon_2s':
      r=17,25                                                             # GoodXenon data contains the channels as written, no need to convert
   else:
      print datamode,'not yet supported, using full range!'
      return data

   words=array(pan.boolval((words[:,r[0]:r[1]]).tolist()))

   mask1=(words>=low)
   mask2=(words<=high)

   ch_data=data[mask1&mask2]

   return ch_data


#-----DiscNEv----------------------------------------------------------------------------------------------------------

def discnev(datas):

   '''Discard Non-Events

   Decription:
    Given a FITS data table, discards all non-photon events from the table.

   -J.M.Court, 2014'''

   mask=datas['Event'][:,0]==True                                         # Creating a mask to obscure any data not labelled as photons
   datas=datas[mask]                                                      # Applying the mask
   return datas


#-----EvMChan----------------------------------------------------------------------------------------------------------

def evmchan(chan):

   '''Event Mode Channel-Get

   Description:

    Converts a channel number into the ID of the range which contains that channel in E_125us_64M_0_1s
    data from PCA on RXTE.  This is the least interesting function I have ever written.

   Inputs:

    chan   - INT: the real channel number for PCA data from RXTE.

   Outputs:

    n_chan - INT: the ID of the range containing the relevant channel.

   -J.M.Court, 2015'''

   chan=int(chan)

   if chan<5:     n_chan=0                                                # This is just a list of ifs.  It checks if the value falls into each and, if not, carries on.
   elif chan<7:   n_chan=1
   elif chan<8:   n_chan=2
   elif chan<9:   n_chan=3
   elif chan<10:  n_chan=4
   elif chan<11:  n_chan=5
   elif chan<12:  n_chan=6
   elif chan<13:  n_chan=7
   elif chan<14:  n_chan=8
   elif chan<15:  n_chan=9
   elif chan<16:  n_chan=10
   elif chan<18:  n_chan=11
   elif chan<20:  n_chan=12
   elif chan<22:  n_chan=13
   elif chan<24:  n_chan=14
   elif chan<26:  n_chan=15
   elif chan<28:  n_chan=16
   elif chan<30:  n_chan=17
   elif chan<32:  n_chan=18
   elif chan<34:  n_chan=19
   elif chan<36:  n_chan=20
   elif chan<38:  n_chan=21
   elif chan<40:  n_chan=22
   elif chan<42:  n_chan=23
   elif chan<44:  n_chan=24
   elif chan<47:  n_chan=25
   elif chan<50:  n_chan=26
   elif chan<53:  n_chan=27
   elif chan<56:  n_chan=28
   elif chan<59:  n_chan=29
   elif chan<62:  n_chan=30
   elif chan<65:  n_chan=31
   elif chan<68:  n_chan=32
   elif chan<72:  n_chan=33
   elif chan<76:  n_chan=34
   elif chan<80:  n_chan=35
   elif chan<84:  n_chan=36
   elif chan<88:  n_chan=37
   elif chan<92:  n_chan=38
   elif chan<97:  n_chan=39
   elif chan<102: n_chan=40
   elif chan<107: n_chan=41
   elif chan<112: n_chan=42
   elif chan<117: n_chan=43
   elif chan<122: n_chan=44
   elif chan<127: n_chan=45
   elif chan<132: n_chan=46
   elif chan<138: n_chan=47
   elif chan<144: n_chan=48
   elif chan<150: n_chan=49
   elif chan<156: n_chan=50
   elif chan<162: n_chan=51
   elif chan<168: n_chan=52
   elif chan<175: n_chan=53
   elif chan<182: n_chan=54
   elif chan<189: n_chan=55
   elif chan<196: n_chan=56
   elif chan<203: n_chan=57
   elif chan<210: n_chan=58
   elif chan<218: n_chan=59
   elif chan<226: n_chan=60
   elif chan<234: n_chan=61
   elif chan<244: n_chan=62
   else:          n_chan=63

   return n_chan


#-----GetBin-----------------------------------------------------------------------------------------------------------

def getbin(event):

   '''Get Bin

   Description:

    Returns the time binning of a given set of FITS data.

   Inputs:

    event - FITS OBJECT: The FITS file that has been opened.

   Outputs:

    bsz   -       FLOAT: The bin size of the observation (s).

   -J.M.Court, 2014'''

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

   -J.M.Court, 2014'''

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

   -J.M.Court, 2014'''

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

   -J.M.Court, 2014'''

   ini=event[1].header['TSTART']
   return ini


#-----Get PCU----------------------------------------------------------------------------------------------------------

def getpcu(words,datamode):

   '''Get PCUs

   Description:

    If given a list of RXTE event words, and the data-mode in which their associated data is stored,
    returns the number of PCUs that contributed photons to the dataset.  Currently works for two
    datamodes: GoodXenon_2s and E_125us_64M_0_1s with the intention to eventually generalise it to
    any event data.

   Inputs:

    words    - 2D ARRAY: The array of event words, each of which is an array of True/False statements.
    datamode -   STRING: The DATAMODE in which the data to be analysed is stored.

   Outputs:

    pcus     -      INT: The number of PCUs active when the data was taken.  If a non-recognised event
                         word is given, 1 is returned instead along with an error message.

   -J.M.Court, 2015'''

   if datamode=='E_125us_64M_0_1s':
      r=1,4
   elif datamode=='GoodXenon_2s':
      r=7,10
   else:
      print datamode,'not yet supported!'
      return 1

   pcus=0

   words=(words[:,r[0]:r[1]]).tolist()

   if [False,False,False] in words: pcus+=1
   if [False,False,True] in words: pcus+=1
   if [False,True,False] in words: pcus+=1
   if [False,True,True] in words: pcus+=1
   if [True,False,False] in words: pcus+=1

   if pcus==0:
      print 'Error!  0 PCUs present in data!'
      return 1
   else:
      return pcus


#-----Get Tim----------------------------------------------------------------------------------------------------------

def gettim(data):

   '''Get Times

   Description: Returns the DATA Word or table column containing data on arrival times.

   -J.M.Court, 2015'''

   return data.field(0) 


#-----Get Wrd----------------------------------------------------------------------------------------------------------

def getwrd(data):

   '''Get Words

   Description: Returns the DATA Word or table column containing data on PCUs.

   -J.M.Court, 2015'''

   return data.field(1)


#-----MaxEn------------------------------------------------------------------------------------------------------------

def maxen():

   '''Max Energy

   Description:
    Returns the highest energy or channel valid for the instrument.

   -J.M.Court, 2015'''

   return 255

