#! /usr/bin/env python

# |----------------------------------------------------------------------|
# |------------------------------DATA FAIRY------------------------------|
# |----------------------------------------------------------------------|

# Call as ./datafairy.py

# A tool to create fake data readable by PlotDemon
#
#

#-----User-set Parameters----------------------------------------------------------------------------------------------

rang=1,4000                                                               # The range of data, in seconds, to be saved
reso=0.1                                                                  # The resolution, in seconds, of that data

## Change the functions below to change the form of the data ##

def band1(t):
   return random.random()

def band2(t):
   return counts1[t]*(sin(times[t]*2*pi/10)+1)

def band3(t):
   return counts1[i]*(cos(times[i]*2*pi/10)+1)


#-----Importing Modules------------------------------------------------------------------------------------------------

from numpy import *
import pan_lib as pan
import random
import pylab as pl


#-----Welcoming Header-------------------------------------------------------------------------------------------------

print ''
print '-------Running Plot Demon: J.M.Court, 2014------'
print ''


#-----Data Creation----------------------------------------------------------------------------------------------------

times=arange(rang[0],rang[1],reso)
counts1=[]
counts2=[]
counts3=[]
for i in range(len(times)): counts1.append(band1(i))                       #############################
for i in range(len(times)): counts2.append(band2(i))  ###--CHANGE FOR FAKE DATA--##
for i in range(len(times)): counts3.append(counts1[i]*(cos(times[i]*2*pi/10)+1))  #############################
errors1=sqrt(counts1)
errors2=sqrt(counts2)
errors3=sqrt(counts2)


#-----Data Saving------------------------------------------------------------------------------------------------------

pan.plotdsv('t__low',times,counts1,errors1,0,reso,[[rang[0],rang[1]]],1,0,'low_ fake','0-1','Fake',('Prithivi','00-00-00-00'),0)
pan.plotdsv('t__med',times,counts2,errors2,0,reso,[[rang[0],rang[1]]],1,0,'med_ fake','2-3','Fake',('Prithivi','00-00-00-00'),0)
pan.plotdsv('t_high',times,counts3,errors3,0,reso,[[rang[0],rang[1]]],1,0,'high fake','4-5','Fake',('Prithivi','00-00-00-00'),0)

print 'Done!'


#-----Footer-----------------------------------------------------------------------------------------------------------

pan.signoff()


