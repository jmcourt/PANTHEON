#!/bin/bash

# |----------------------------------------------------------------------|
# |------------------------------FITS GX-EX------------------------------|
# |----------------------------------------------------------------------|

# Call as ./fitsgxex.sh

# Requirements:
#  fkeyprint and make_se from FTools: http://heasarc.gsfc.nasa.gov/ftools/

# Takes raw RXTE PCA data in the current working directory, searches it for any labelled as GoodXenon
# data, and copies it into a subdirectory named 'gx0'.  Event data is then made from these files
# and saved as EVENT_gx#.fits using make_se
#
# WARNING!  If run in a directory containing a subdirectory named 'gx0', the contents of that 
# subdirectory WILL BE REMOVED.
#


#-----Welcoming Header-------------------------------------------------------------------------------------------------

echo ''
echo '-------Running FITS Gx-Ex: J.M.Court, 2015-------'
echo ''


#-----Setting Switches-------------------------------------------------------------------------------------------------

gx1=0                                                                     # Counter of GoodXenon1 events
gx2=0                                                                     # Counter of GoodXenon2 events


#-----Creating gx0 Directory-------------------------------------------------------------------------------------------

echo 'Sorting FITS files, retrieving Good Xenon Data products'
echo ''
echo 'Creating directories:'

if [ -d gx0 ]
then
   echo "Directory 'gx0' found: removing contents."
   rm gx0/*                                                               # Empty gx0 directory if it is found
else
   mkdir gx0                                                              # Make gx0 directory if it is not found
fi

echo ''


#-----Seek and Copy Standard Files-------------------------------------------------------------------------------------

for f in *;
do
   dm=$(fkeyprint $f+1 DATAMODE 2> /dev/null )                            # Collecting DATAMODE from FITS file, suppressing error from opening non FITS files
   if [[ $dm =~ .*GoodXenon1.* ]]
   then                                                                   # Switching ce to 'on' (i.e. Good Xenon file found)
      echo $f 'is GoodXenon1 data'
      gx1=$((gx1+1))                                                      # Increasing GoodXenon1 counter
      cp $f ./gx0/GX1_$gx1.fits                                           # If DATAMODE is GoodXenon1, move this file to /gx0
   elif [[ $dm =~ .*GoodXenon2.* ]]
   then
      echo $f 'is GoodXenon2 data'
      gx2=$((gx2+1))                                                      # Increasing GoodXenon1 counter
      cp $f ./gx0/GX2_$gx2.fits                                           # If DATAMODE is GoodXenon2, move this file to /gx0
   fi
done


#-----Return Warnings and Clean Up if GoodXenon not Found--------------------------------------------------------------

if [ $gx1 == 0 ] || [ $gx2 == 0 ]
then
  rm -r gx0                                                               # Remove gx0 if no Good Xenon data found
  echo 'No matching GoodXenon Data files found!'
  echo 'Cannot make good xenon files!'
  echo 'Removing gx0 directory'


#-----Make GoodXenon File from GoodXenon files-------------------------------------------------------------------------

else
  
  cd gx0

  if [ $gx2 -gt $gx1 ]                                                    # Picking out minimum value of gx1, gx2
  then
     echo 'Mismatch warning!  More Xenon2 than Xenon1!'
     end=$gx1
  elif [ $gx1 -gt $gx2 ]
  then
     echo 'Mismatch warning!  More Xenon1 than Xenon2!'
     end=$gx2
  else
     end=$gx1
  fi
     
  for i in $(seq 1 $end)
  do
     echo 'GX1_'$i'.fits
        GX2_'$i'.fits' >> temp.txt                                        # Writing input file list required by make_se
  done
  echo 'temp.txt
     EVENT' > temp2.txt                                                   # Input for make_se
  make_se <temp2.txt                                                      # Make gx file as EVENT_gx#.fits
  rm temp*.txt

fi


#-----Footer-----------------------------------------------------------------------------------------------------------

echo ''
echo '------------------------------------------------'
echo ''


