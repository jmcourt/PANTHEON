#!/bin/bash

# |----------------------------------------------------------------------|
# |------------------------------XTE EV-EX-------------------------------|
# |----------------------------------------------------------------------|

# Call as ./xteevex.sh

# Requirements:
#  None!

# Takes raw RXTE PCA data in the current working directory, searches it for any labelled as .evt Event
# data, and copies it into a subdirectory named 'event0' as EVENT_ev#.fits
#
# WARNING!  If run in a directory containing a subdirectory named 'event0', the contents of that 
# subdirectory will be removed.
#


#-----Welcoming Header-------------------------------------------------------------------------------------------------

echo ''
echo '-------Running XTE Ev-Ex: J.M.Court, 2014--------'
echo ''


#-----Setting Switches-------------------------------------------------------------------------------------------------

ev0=0                                                                     # Counter of Event data files


#-----Creating event0 Directory----------------------------------------------------------------------------------------

echo 'Sorting FITS files, retrieving Event Data products'
echo ''
echo 'Creating directories:'

if [ -d event0 ]
then
   echo "Directory 'event0' found: removing contents."
   rm event0/*                                                            # Empty event0 directory if it is found
else
   mkdir event0                                                           # Make event0 directory if it is not found
fi

echo ''


#-----Seek and Copy .evt Event Files-----------------------------------------------------------------------------------

for f in *.evt;
do
   echo $f 'is Event data'
   cp $f ./event0/EVENT_ev$ev0                                         # If file extension is evt, move this file to /event0
   ev0=$((ev0+1))                                                      # Increasing Event counter

done


#-----Return Warnings and Clean Up if Event not Found------------------------------------------------------------------

if [ $ev0 == 0 ]
then
  rm -r event0                                                            # Remove event0 if no event data found
  echo 'No matching Event Data files found!'
  echo 'Removing event0 directory'
fi


#-----Footer-----------------------------------------------------------------------------------------------------------

echo ''
echo '------------------------------------------------'
echo ''


