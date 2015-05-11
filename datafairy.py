#! /usr/bin/env python

# |----------------------------------------------------------------------|
# |------------------------------DATA FAIRY------------------------------|
# |----------------------------------------------------------------------|

# Call as ./datafairy.py

# A tool to create fake data readable by PlotDemon
#
#

#-----User-set Parameters----------------------------------------------------------------------------------------------

rang=0,100000                                                             # The range of data, in seconds, to be saved
reso=29.3023                                                              # The resolution, in seconds, of that data

## Change the functions below to change the form of the data ##

def band1(t):
   return sin(2*pi*t/16.12)+2+0.2*random.random()+0.3*sin(2*pi*t/20.04)+1/(abs(cos(2*pi*t/40))+0.2)

def band2(t):
   return 0

def band3(t):
   return 0


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
for i in times: counts1.append(band1(i))
for i in times: counts2.append(band2(i))
for i in times: counts3.append(band3(i))
errors1=sqrt(counts1)
errors2=sqrt(counts2)
errors3=sqrt(counts2)

pl.plot(times,counts1)
pl.show(block=True)
import pan_lib as pan

pan.circfold(times,counts1,40)

#-----Data Saving------------------------------------------------------------------------------------------------------

pan.plotdsv('t_0-14'  ,times,counts1,errors1,0,reso,[[rang[0],rang[1]]],1,0,'low_ fake','0-1','Fake',('None','00-00-00-00'),0)
pan.plotdsv('t_15-35' ,times,counts2,errors2,0,reso,[[rang[0],rang[1]]],1,0,'med_ fake','2-3','Fake',('None','00-00-00-00'),0)
pan.plotdsv('t_36-255',times,counts3,errors3,0,reso,[[rang[0],rang[1]]],1,0,'high fake','4-5','Fake',('None','00-00-00-00'),0)

print 'Done!'


#-----Footer-----------------------------------------------------------------------------------------------------------

pan.signoff()


