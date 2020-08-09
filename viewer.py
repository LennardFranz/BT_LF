#! /usr/bin/env python3

from __future__ import division
from __future__ import print_function
ls

import os, sys, inspect

cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile(inspect.currentframe()))[0],"src")))
if cmd_subfolder not in sys.path:
  sys.path.insert(0, cmd_subfolder)

from TLeCroy import *

# --------------------------------
# constants to be set by the user
# --------------------------------
channelNr = 4

traceNr = 18


 
#inName = "--testpulse--"
#inName = "--testpulse--autosave--"
inName = "--testpulse--autosave--segments--"


#inDir = "/Users/Lennard/Bachelor_Thesis_LF/BA_Lennard_Franz/cs137-longzylinder/"

inDir = "cs137-longzylinder/"


#inDir = "/Users/Lennard/desktop/osci/lennard/segment-flat-CeBr-CS137/"
#inDir = "/Users/Lennard/desktop/osci/lennard/segment-CeBr-cs137/"
#inDir = "/Users/Lennard/desktop/osci/lennard/cebr-big-4ch-cs137/"
#inDir = "/volumes/lennad/Lennard Franz/cs137-longzylinder-shaper/"
#inDir = "/Users/Lennard/desktop/osci/lennard/cs137-longzylinder/"
#inDir = "/Users/Lennard/desktop/osci/lennard/cs137-filter/"
#inDir = "/Users/Lennard/desktop/autosave-segments/"
#inDir = "/Users/Lennard/Desktop/osci/lennard/cebr-cs137/"

#inDir = "/home/lfranz/work/osci/lennard/segment-longcable-CeBr-CS137/"
#inDir = "/home/lfranz/work/osci/lennard/segment-flat-plexi-background/"
#inDir = "/home/lfranz/work/osci/lennard/segment-CeBr-ra226/"
#inDir = "/home/lfranz/work/osci/lennard/segment/"
#inDir = "/home/lfranz/work/autosave-segments/"
#inDir = "/ZIH.fast/projects/d-lab_data/MAPMT/2019_gsi/segment/"
#inDir = "/ZIH.fast/projects/d-lab_data/MAPMT/2019_gsi/segments-new/"
#inDir = "/ZIH.fast/projects/d-lab_data/MAPMT/2019_gsi/segments-new/Autosave-new/"
# inDir = "/ZIH.fast/projects/d-lab_data/MAPMT/2019_gsi/Autosave/"
# --------------------------------

fileName = inDir + "C" + str(channelNr).zfill(1) + inName + str(traceNr).zfill(5) + ".trc"
TLC = TLeCroy(fileName, debug=True)
TLC.PrintPrivate()
#TLC.PlotTrace(raw=False)

header = TLC.GetHeader()
x, y = TLC.GetTrace()
y = y - header['VERTICAL_OFFSET']
#y = -y
nano = 1e-9

import numpy as np
import scipy.integrate as integrate
 
#print(header['VERTICAL_OFFSET'])
#print('   vertical gain:', header['VERTICAL_GAIN'])
#print('----------------')
seqLength = header['WAVE_ARRAY_COUNT'] // header['SUBARRAY_COUNT']
seqFirst = header['FIRST_VALID_PNT']
print('  LENGTH', seqLength)
print(len(x)/1000,len(y)/1000)


import matplotlib.pyplot as plt
#plt.plot(x / nano, y)


idxHigh = np.zeros(header['SUBARRAY_COUNT'],dtype=np.int)
idxLow = np.zeros(header['SUBARRAY_COUNT'],dtype=np.int)
mean = np.zeros(header['SUBARRAY_COUNT'])
std = np.zeros(header['SUBARRAY_COUNT'])
area =  np.zeros(header['SUBARRAY_COUNT'])
area_all = np.zeros(header['SUBARRAY_COUNT'])
A = np.zeros(header['SUBARRAY_COUNT'])
area_bool = np.zeros(header['SUBARRAY_COUNT'])
area_mean = np.zeros(header['SUBARRAY_COUNT'])
area_all = np.zeros(header['SUBARRAY_COUNT'])

for seqNr in range(1, header['SUBARRAY_COUNT'] + 1):
  #bondarys of sequences
  idxHigh[seqNr-1] = seqFirst +( seqLength * seqNr - 1)
  idxLow[seqNr-1] = idxHigh[seqNr-1] - seqLength + 1
  #finding basline
  #y[...] has to be adjusted depending vertical position of the puls                     
  from scipy.stats import norm
  mean[seqNr-1], std[seqNr-1] = norm.fit(y[ idxLow[seqNr-1] : idxLow[seqNr-1] + int(0.35*seqLength)])

  #print(seqNr, idxLow, idxHigh)
mean_all = sum(mean)/header['SUBARRAY_COUNT']
std_all = sum(std)/header['SUBARRAY_COUNT'] 
#print(mean_all,std_all)
#plt.hist(mean, bins=100)


#plt.hist(y[idxLow[seqNr-1] : idxLow[seqNr-1]+int(0.35*seqLength)])


for seqNr in range(1, header['SUBARRAY_COUNT'] + 1):
 # if mean[seqNr -1] > 0.0285 and mean[seqNr - 1] < 0.0315 and np.min( y[idxLow[seqNr-1] : idxHigh[seqNr-1]]-mean[seqNr-1]) > -0.0035:
    
    #area of the pulse
  area[seqNr-1] = integrate.trapz((y[idxLow[seqNr - 1]:idxHigh[seqNr - 1]]- mean[seqNr-1]) > 2*std[seqNr-1] ,
                  x[idxLow[seqNr - 1]:idxHigh[seqNr - 1]]/nano)  
 
    # area_all[seqNr-1] = integrate.trapz((y[idxLow[seqNr - 1]:idxHigh[seqNr - 1]] - mean_all) > 2*std_all ,
                     # x[idxLow[seqNr - 1]:idxHigh[seqNr - 1]]/nano)
    
  #boolArr = np.where((y[idxLow[seqNr - 1]:idxHigh[seqNr - 1]]- mean[seqNr-1]) > (3*std[seqNr-1]))
  boolArr = (y[idxLow[seqNr - 1]:idxHigh[seqNr - 1]]- mean[seqNr-1]) > (3*std[seqNr-1])
  #print((boolArr))
  #area_bool[seqNr-1] = np.sum(y[boolArr[0]])*(header['HORIZ_INTERVAL']/nano) - (mean[seqNr-1])*(len(boolArr[0])*header['HORIZ_INTERVAL']/nano)  
  #area_bool[seqNr-1] = np.sum(y[boolArr[0]])
  area_bool[seqNr-1] = np.sum((y[idxLow[seqNr - 1]:idxHigh[seqNr - 1]]*(header['HORIZ_INTERVAL']/nano)),where = boolArr) - mean[seqNr-1]*np.sum(boolArr)*header['HORIZ_INTERVAL']/nano


  #print(np.trapz(y[boolArr[0]],x[boolArr[0]])/nano,area_bool[seqNr-1])
  #print(np.sum(y[idxLow[seqNr-1]:idxHigh[seqNr-1]]>(mean[seqNr-1]+3*std[seqNr-1])),len(boolArr[0])) 
  area_bool_positive = area_bool > 0.0 
  
  
                    
  #print(seqNr, idxLow, idxHigh)
  #area = integrate.trapz(y[idxLow:idxHigh],x[idxLow:idxHigh]/nano)
  #print(area)
  if seqNr > 18  and seqNr < 20 :
    #plt.plot(x[idxLow[seqNr - 1]:idxHigh[seqNr - 1]]  / nano, y[idxLow[seqNr - 1]:idxHigh[seqNr - 1]]-mean[seqNr-1], linewidth=0.5, color='k')
    plt.plot(x[idxLow[0]:idxHigh[0]]  / nano, (y[idxLow[seqNr - 1]:idxHigh[seqNr - 1]]-mean[seqNr-1])*1000, linewidth=1, color='k')
  
  #area_mean[seqNr-1]= mean[seqNr-1]*len(boolArr[0])*header['HORIZ_INTERVAL']/nano
  area_all[seqNr-1] = np.sum((y[idxLow[seqNr - 1]:idxHigh[seqNr - 1]]*(header['HORIZ_INTERVAL']/nano)),where = boolArr)
  #print(area_all[seqNr-1])


#plt.hist(area_mean)
#plt.hist(area_all)
#print(np.sum(y[idxLow[seqNr-1]:idxHigh[seqNr-1]]>(mean[seqNr-1]+3*std[seqNr-1])))
#print(area_bool)


plt.xlabel('time (ns)')
plt.ylabel('voltage (mV)')
plt.grid()
#plt.plot(x[idxLow[428]:idxHigh[428]]/nano , y[idxLow[428]:idxHigh[428]]-mean[428])

#x=np.linspace(0,1000.1e-9,num=10001,endpoint=False)
#print(x)
#plt.plot(x , y[idxLow[3]:idxHigh[3]]-mean[3])
#area_1 = integrate.trapz((y[idxLow[1]:idxHigh[1]] - mean_all) - 1*std_all,  x[idxLow[1]:idxHigh[1]]/nano)
    
            
#print(area_1)

#print("mean all:",mean_all,"std all:", std_all)

#print(mean,std) 

#print("-----------------------")
#print("area corr" , area)
#print("std" , std)
#print("mean" , mean)
#print("area uncorr" , #area_1)
#print("-----------------------")

plt.xlim(-50,150)



#ax1.hist(area,bins=50)
mean_E,std_E = norm.fit(area) 
#print("---------------------")
#print("mean:", mean_E, "std:", std_E)
#print("---------------------")

#plt.hist(area_bool ,bins=50)
#mean_E_all, std_E_all = norm.fit(area_bool)
#print("---------------------")
#print("mean:", mean_E_all, "std:", std_E_all)
#print("---------------------")


#import tikzplotlib
#tikzplotlib.save("pulse_mit_shaper.tex")



plt.show()

