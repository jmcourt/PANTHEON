#! /usr/bin/env python

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
#  BIHCHAN   - converts an RXTE channel ID into a channel range ID.
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

# 'E'
# 'B''
# 'SB'
# 'GoodXenon' 


#-----Importing Modules------------------------------------------------------------------------------------------------

import pan_lib as pan
from numpy import array, zeros, arange
from numpy import sum as npsum


#-----ChRange----------------------------------------------------------------------------------------------------------

def chrange(data,low,high,header):

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

   datamode=header['DATAMODE']

   if datamode[0]=='B':

      chconv=header['TDDES2'].split('&')[2].split('[')[1].split(']')[0].split(',')
      low=bihchan(datamode,low,chconv)
      high=bihchan(datamode,high,chconv)+1

      ndat=zeros(len(data[0]))

      for i in range(low,high):
         ndat+=array(data[i])

      return ndat

   elif datamode[:2]=='SB':

      print 'No energy information in this datamode!'                     # No energy information in SB datamodes, so don't filter
      return data

   else:

      words=data.field(1)


      if low<=0 and high>=255:
         return data                                                      # Don't bother searching through if the user wants full range

      if datamode[:10]=='GoodXenon_':
         r=(17,25)                                                        # GoodXenon data contains the channels as written, no need to convert

      elif datamode[:2]=='E_':

         enigma=header['TEVTB2'].split('^')[0][1:-1].split('}')[:-1]      # Fetch (and decode) the enigmatic TEVTB2 word
         tk,r=pan.tokenloc(enigma,'C')

         if tk==None:
            print 'No energy information in this datamode!'               # If no channels token is present, can't filter on energy
            return data

         low=bihchan(datamode,low,tk)
         high=bihchan(datamode,high,tk)+1

      else:
         print datamode,'not yet supported, using full range!'
         return data

      words=array(pan.boolval((words[:,r[0]:r[1]]).tolist()))
      mask1=(words>=low)

      mask2=(words<=high)

      ch_data=data[mask1&mask2]

      return ch_data


#-----DiscNEv----------------------------------------------------------------------------------------------------------

def discnev(datas,datamode):

   '''Discard Non-Events

   Decription:
    Given a FITS data table, discards all non-photon events from the table.

   -J.M.Court, 2014'''

   if datamode[:2] in ['B_','SB']:
      return discnevb(datamode,datas.field(1))

   mask=datas['Event'][:,0]==True                                         # Creating a mask to obscure any data not labelled as photons
   datas=datas[mask]                                                      # Applying the mask
   return datas


#-----DiscNEvB---------------------------------------------------------------------------------------------------------

def discnevb(datamode,datas):

   '''Discard Non-Events: Binned Data Version

   Decription:
    Given a FITS binned data table, discards all non-photon events from the table.

   -J.M.Court, 2014'''

   if datamode[:2]=='SB':

      datas=datas.reshape([len(datas)*len(datas[0])])
      return datas

   chan={}

   nbchan=int(datamode.split('_')[2][:-1])

   for j in range(nbchan):
      chan[j]=[]

   for i in range(len(datas)):
      for j in range(nbchan):
         chan[j]+=datas[i][j].tolist()

   outchan=[]

   for j in range(nbchan):
      outchan.append(chan[j])

   return outchan


#-----BiHChan----------------------------------------------------------------------------------------------------------

def bihchan(datamode,chan,chanconv):

   '''Channel-Get

   Description:

    Converts a channel number into the ID of the range which contains that channel in B_2ms_4B_0_35_H or B_8ms_16A_0_35_H
    data from PCA on RXTE.

   Inputs:

    chan   - INT: the real channel number for PCA data from RXTE.

   Outputs:

    n_chan - INT: the ID of the range containing the relevant channel.

   -J.M.Court, 2015'''


   n_nchan=len(chanconv)
   chan=int(chan)
   t_chan=[0]*n_nchan

   for i in range(len(chanconv)):

      t_chan[i]=map(int,chanconv[i].split('~'))

   if chan>t_chan[-1][-1]:                                                # Simple sanity check to prevent messy accidents

      print 'This data type does not store photons above Channel '+str(t_chan[-1][-1])+'!'
      pan.signoff()
      exit()

   elif chan<t_chan[0][0]:

      print 'This data type does not store photons below Channel '+str(t_chan[0][0])+'!'
      pan.signoff()
      exit()

   for i in range(len(t_chan)):
      if chan<=t_chan[i][-1]:                                             # This is just a list of ifs.  It checks if the value falls into each and, if not, carries on.
         return i


#-----GetBG------------------------------------------------------------------------------------------------------------

@pan.mjit()
def getbg(data,low_channel,high_channel):

   '''Get BG

   Description:

    Returns the counts from a standard XTE Background file created with pcabackest

   Inputs:

    event        - FITS OBJECT: The data array of FITS file that has been opened.
    low_channel  -         INT: The highest channel data is to be collected from.
    high_channel -         INT: The lowest channel data is to be collected from.

   Outputs:

    times  -              LIST: The time array associated with the fluxes
    fluxes -              LIST: The fluxes at each time

   -J.M.Court,2015'''

   times=data.field(0)

   PCU0=data.field(1)                                                     # Collect data from PCU0
   PCU1=data.field(2)                                                     # Collect data from PCU1
   PCU2=data.field(3)                                                     # Collect data from PCU2
   PCU3=data.field(4)                                                     # Collect data from PCU3
   PCU4=data.field(5)                                                     # Collect data from PCU4
   PCU=array([PCU0,PCU1,PCU2,PCU3,PCU4])
   PCU[:,:,:low_channel]=0                                                # Remove data from outside of channel arrays
   PCU[:,:,(high_channel+1):]=0
   PCU=npsum(PCU,axis=2)
   fluxes=npsum(PCU,axis=0)
   flux_ers=fluxes**0.5
   fluxes=fluxes/16.0                                                     # Change from counts to counts/s
   flux_ers=flux_ers/16.0                                                 # Change from counts to counts/s

   return times,fluxes,flux_ers

#-----GetBin-----------------------------------------------------------------------------------------------------------

def getbin(event,datamode):

   '''Get Bin

   Description:

    Returns the time binning of a given set of FITS data.

   Inputs:

    event - FITS OBJECT: The FITS file that has been opened.

   Outputs:

    bsz   -       FLOAT: The bin size of the observation (s).

   -J.M.Court, 2014'''

   if datamode[:2] in ['B_','SB']:
      bsz=event[1].header['1CDLT2']
   else:
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


#-----Get Obs----------------------------------------------------------------------------------------------------------

def getobs(event,datamode,filepath):

   '''Get Obs

   Description:

    Fetches a tuple consisting of the object and obs_id of the observation.'''

   if datamode[:10]=='GoodXenon_':
      try:
         obsid=(filepath.split('/')[-4])                                  # GoodXenon for XTE doesnt store obs_id for some reason
      except:
         obsid=''
   else:
      try:
         obsid=event[1].header['OBS_ID']
      except:
         obsid=''

   obsdata=(event[1].header['OBJECT'],obsid)

   return obsdata


#-----Get PCU----------------------------------------------------------------------------------------------------------

def getpcu(words,header,t_pcus=None,pculist=False):

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

   datamode=header['DATAMODE']

   if datamode[0]=='E':
      enigma=header['TEVTB2'].split('^')[0][1:-1].split('}')[:-1]      # Again fetch (and decode) the enigmatic TEVTB2 word
      tk,r=pan.tokenloc(enigma,'D')
   elif datamode=='GoodXenon_2s':
      r=7,10
   else:
      r=None
   if r==None and t_pcus==None:
      goodpcus=False

      while not goodpcus:
         try:
            pcus=int(raw_input('Number of Active PCUS: '))             # Ask the user how many there are
            assert pcus<6                                              # Check they give a number between 1 and 5 (inclusive)
            assert pcus>0
            goodpcus=True
         except:
            'Invalid number of PCUs!'

      return pcus                                                      # Use the number they give as the number of PCUs

   elif t_pcus!=None:

      return t_pcus

   words=(words[:,r[0]:r[1]]).tolist()

   pcus=[0,0,0,0,0]

   if [False,False,False] in words: pcus[0]=1
   if [False,False,True] in words: pcus[1]=1
   if [False,True,False] in words: pcus[2]=1
   if [False,True,True] in words: pcus[3]=1
   if [True,False,False] in words: pcus[4]=1

   if pcus==0:
      print 'Error!  0 PCUs present in data!'
      return 1
   elif not pculist:
      return sum(pcus)
   else:
      return (sum(pcus), pcus)


#-----Get Tim----------------------------------------------------------------------------------------------------------

def gettim(data,event,tstart,res,datamode):

   '''Get Times

   Description: Returns the DATA Word or table column containing data on arrival times.

   -J.M.Court, 2015'''

   if datamode[0]=='B':
      times=[]
      for i in range(len(data)):
         times+=[(i*res)+tstart]
      return array(times),data

   elif datamode[:2]=='SB':
      data=data.tolist()
      timeseeds=event.field(0)-event.field(0)[0]+tstart
      times=[]
      indx=0
      prevtimeseed=timeseeds[0]
      for i in timeseeds:                                                
         nx=prevtimeseed+1
         while i-nx>0.9:                                                  # When a 1s chunk is missing, zero pad the data here
            times+=(arange(0,1.0,res)+nx).tolist()
            data=data[:indx]+([(data[indx]+data[indx-1])/2.0]*int(1.0/res))+data[indx:]
            nx+=1
            indx+=int(1.0/res)
         for j in arange(0,1.0,res):
            times+=[i+j]
         indx+=int(1.0/res)
         prevtimeseed=i
      return array(times),array(data)
      

   else:
      return data.field(0),data 


#-----Get Wrd----------------------------------------------------------------------------------------------------------

def getwrd(data,datamode):

   '''Get Words

   Description: Returns the DATA Word or table column containing data on PCUs.

   -J.M.Court, 2015'''

   if datamode[:2] in ['B_','SB']:
      return None
   else:
      return data.field(1)


#-----Get Wrd Row------------------------------------------------------------------------------------------------------

def getwrdrow(words,mask,datamode):

   '''Get Word Rows

   Description: Returns Data Words filtered by a mask, for data that has datawords.

   -J.M.Court, 2015'''

   if datamode[:2] in ['B_','SB']:
      return None
   else:
      return words[mask]


#-----MaxEn------------------------------------------------------------------------------------------------------------

def maxen(datamode):

   '''Max Energy

   Description:
    Returns the highest energy or channel valid for the instrument.

   -J.M.Court, 2015'''

   if datamode[:2] in ['B_','SB']:
      mcha=datamode.split('_')[-2]
      return int(mcha)
   else:
      return 255


