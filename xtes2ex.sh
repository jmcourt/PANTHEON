#!/bin/bash

# |----------------------------------------------------------------------|
# |------------------------------XTE S2-EX-------------------------------|
# |----------------------------------------------------------------------|

# Call as ./xtegxex.sh

# Requirements:
#  fkeyprint and make_se from FTools: http://heasarc.gsfc.nasa.gov/ftools/

# Takes raw RXTE PCA data in the current working directory, searches it for any labelled as GoodXenon
# data, and copies it into a subdirectory named 'gx0'.  Event data is then made from these files
# and saved as EVENT_gx#.fits using make_se
#
# WARNING!  If run in a directory containing a subdirectory named 's20', the contents of that 
# subdirectory will be removed.
#


#-----Welcoming Header-------------------------------------------------------------------------------------------------

echo ''
echo '-------Running XTE S2-Ex: J.M.Court, 2015--------'
echo ''


#-----Setting Switches-------------------------------------------------------------------------------------------------

s2=0                                                                      # Counter of GoodXenon1 events


#-----Creating s20 Directory-------------------------------------------------------------------------------------------

echo 'Sorting FITS files, retrieving Standard2 Data products'
echo ''
echo 'Creating directories:'

if [ -d s20 ]
then
   echo "Directory 's20' found: removing contents."
   rm s20/*                                                               # Empty gx0 directory if it is found
else
   mkdir s20                                                              # Make gx0 directory if it is not found
fi

echo ''


#-----Seek and Copy Standard Files-------------------------------------------------------------------------------------

for f in *;
do
   dm=$(fkeyprint $f+1 DATAMODE 2> /dev/null )                            # Collecting DATAMODE from FITS file, suppressing error from opening non FITS files
   if [[ $dm =~ .*Standard2.* ]]
   then
      echo $f 'is Standard 2 data'                                        # Increasing Standard2 counter
      cp $f ./s20/$f
      rm *.bg >/dev/null                                                  # Kill it if its a bg file
      if [[ -a $f ]];
      then
         mv ./s20/$f ./s20/EVENT_s2$s2.fits                               # If DATAMODE is Standard2, move this file to /s2
         s2=$((s2+1)) 
      fi
   fi
done


#-----Return Warnings and Clean Up if GoodXenon not Found--------------------------------------------------------------

if [ $s2 == 0 ]
then
  rm -r s20                                                               # Remove s2 if no Standard2 data found
  echo 'No Standard2 files found!'
  echo 'Cannot make good xenon files!'
  echo 'Removing s2 directory'
else                                                                      # Remove any stray background files
  rm s20/*.bg
fi


#-----Footer-----------------------------------------------------------------------------------------------------------

echo ''
echo '------------------------------------------------'
echo ''


