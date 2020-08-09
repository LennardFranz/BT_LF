#! /usr/bin/env python3

from __future__ import division
from __future__ import print_function

import os, sys, inspect

cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile(inspect.currentframe()))[0],"src")))
if cmd_subfolder not in sys.path:
  sys.path.insert(0, cmd_subfolder)

from TLeCroy import *   # lecroy.py 
import matplotlib.pyplot as plt
from scipy.stats import norm
import numpy as np
import scipy.integrate as integrate
import matplotlib.mlab as mlab
from scipy.signal import find_peaks

# --------------------------------
# constants to be set by the user
# --------------------------------
channelNr = 4

 
StartTrace = 0
StopTrace = 25



#inName = "--testpulse--"
#inName = "--testpulse--autosave--"
inName = "--testpulse--autosave--segments--"

inDir = "/Users/Lennard/Bachelor_Thesis_LF/BA_Lennard_Franz/cs137-longzylinder/"

#inDir = "/volumes/lennad/segment/segment-CeBr-cs137/"
#inDir = "/volumes/lennad/Lennard Franz/cs137-longzylinder/"
#inDir = "/volumes/lennad/Lennard Franz/th232-cebr-4x4/"
#inDir = "/Users/Lennard/desktop/osci/lennard/segment-flat-CeBr-CS137/"
#inDir = "/Users/Lennard/desktop/osci/lennard/cs137-longzylinder-shaper/"
#inDir = "/Users/Lennard/destop/osci/lennard/cebr-big-4ch-cs137/"
#inDir = "/Users/Lennard/desktop/osci/lennard/segment-CeBr-am241/"
#inDir = "/Users/Lennard/desktop/osci/lennard/ra226-filter/"
#inDir = "/Users/Lennard/desktop/autosave-segments/"
#inDir = "/Users/Lennard/desktop/osci/lennard/gag-am241/"


#inDir = "/home/lfranz/work/osci/lennard/segment-longcable-CeBr-CS137/"
#inDir = "/home/lfranz/work/osci/lennard/segment-flat-CeBr-CS137/"
#inDir = "/home/lfranz/work/osci/lennard/segment-flat-plexi-background/"
#inDir = "/home/lfranz/work/osci/lennard/segment-CeBr-am241/"
#inDir = "/home/lfranz/work/osci/lennard/segment-CeBr-cs137/"
#inDir = "/home/lfranz/work/osci/lennard/segment-CeBr-ra226/"

#inDir = "/home/lfranz/work/osci/lennard/segment/"
#inDir = "/home/lfranz/work/autosave-segments/"
#inDir = "/ZIH.fast/projects/d-lab_data/MAPMT/2019_gsi/segment/"
#inDir = "/ZIH.fast/projects/d-lab_data/MAPMT/2019_gsi/segments-new/"
#inDir = "/ZIH.fast/projects/d-lab_data/MAPMT/2019_gsi/segments-new/Autosave-new/"
# inDir = "/ZIH.fast/projects/d-lab_data/MAPMT/2019_gsi/Autosave/"
# --------------------------------

import numpy as np
import scipy.integrate as integrate

#empty array for area of sequnces for every trace
A = np.empty([0])
print(inDir)
for traceNr in range(StartTrace, StopTrace + 1, 1):
  print(" INFO: working on trace: ", "C"+ str(channelNr).zfill(1) + inName + str(traceNr).zfill(5) + ".trc")
 
  fileName = inDir + "C" + str(channelNr).zfill(1) + inName + str(traceNr).zfill(5) + ".trc"
  TLC = TLeCroy(fileName, debug=True)
  #TLC.PrintPrivate()
  # TLC.PlotTrace(raw=False)

  header = TLC.GetHeader()
  x, y = TLC.GetTrace()
  y = y - header['VERTICAL_OFFSET']
  y = -y 
  nano = 1e-9

  #print('----------------')
  #print('   vertical gain:', header['VERTICAL_GAIN'])
  #print('----------------')
  seqLength = header['WAVE_ARRAY_COUNT'] // header['SUBARRAY_COUNT']
  seqFirst = header['FIRST_VALID_PNT']
  #print('  LENGTH', seqLength)

  #plt.plot(x / nano, y)
  
  #arrays  
  idxHigh = np.zeros(header['SUBARRAY_COUNT'],dtype=np.int)
  idxLow = np.zeros(header['SUBARRAY_COUNT'],dtype=np.int)
  mean = np.zeros(header['SUBARRAY_COUNT'])
  std = np.zeros(header['SUBARRAY_COUNT'])
  area =  np.zeros(header['SUBARRAY_COUNT'])
  area_all = np.zeros(header['SUBARRAY_COUNT'])
  area_all_sum = np.zeros(header['SUBARRAY_COUNT'])
  area_bool =np.zeros(header['SUBARRAY_COUNT'])
  
# determine background
  for seqNr in range(1, header['SUBARRAY_COUNT'] + 1):
    #bondarys of sequences
    idxHigh[seqNr-1] = seqFirst +( seqLength * seqNr - 1)
    idxLow[seqNr-1] = idxHigh[seqNr-1] - seqLength + 1
    #finding basline
    #y[...] has to be adjusted depending vertical position of the puls                     
    mean[seqNr-1], std[seqNr-1] = norm.fit(y[ idxLow[seqNr-1] : idxLow[seqNr-1] + int(0.35*seqLength)])

  mean_all = sum(mean)/header['SUBARRAY_COUNT']
  std_all = sum(std)/header['SUBARRAY_COUNT'] 
  #Histogram of values for baseline 
  #plt.hist(y[ idxLow[seqNr-1] : idxLow[seqNr-1] + int(0.35*seqLength)])
  

#determination of area of pulses
  for seqNr in range(1, header['SUBARRAY_COUNT'] + 1):
   # if mean[seqNr -1] > 0.0285 and mean[seqNr - 1] < 0.0315 and np.min( y[idxLow[seqNr-1] : idxHigh[seqNr-1]]-mean[seqNr-1]) > -0.0035:

#area of the pulse
    #area[seqNr-1] = np.trapz((y[idxLow[seqNr - 1]:idxHigh[seqNr - 1]]- mean[seqNr-1]) > 3*std[seqNr-1] ,
                 #x[idxLow[seqNr - 1]:idxHigh[seqNr - 1]]/nano)  
 
    #area_all[seqNr-1] = np.sum(y[idxLow[seqNr - 1] : idxHigh[seqNr - 1]])/(nano/header['HORIZ_INTERVAL']) - ((mean_all)*seqLength*header['HORIZ_INTERVAL']/nano)
   
    #print(area_all_sum)
    #print(area_all)
    #condition to differentiate values of pulse from background 
    #boolArr = np.where((y[idxLow[seqNr - 1]:idxHigh[seqNr - 1]]- mean[seqNr - 1]) > (3*std[seqNr -1]))
    boolArr = (y[idxLow[seqNr - 1]:idxHigh[seqNr - 1]]- mean[seqNr-1]) > (3*std[seqNr-1])

    #sum of all values which meet condition and scaling 
    #substracting of area belonging to background
    
    #print(y[boolArr[0]] - mean[seqNr-1])
    #area_bool[seqNr-1] = np.sum(y[boolArr[0]])*(header['HORIZ_INTERVAL']/nano) - (mean[seqNr-1])*(len(boolArr[0])*header['HORIZ_INTERVAL']/nano)
    area_bool[seqNr-1] = np.sum((y[idxLow[seqNr - 1]:idxHigh[seqNr - 1]]*(header['HORIZ_INTERVAL']/nano)),where = boolArr) - mean[seqNr-1]*np.sum(boolArr)*header['HORIZ_INTERVAL']/nano

    #area_bool[seqNr-1] = np.sum(y[boolArr[0]]) 
     #- (mean[seqNr-1])*(len(boolArr[0]))


    #print(len(boolArr[0])*header['HORIZ_INTERVAL']/nano)    
    #plt.hist(len(boolArr[0])*header['HORIZ_INTERVAL']/nano)

    #plt.plot(x[idxLow[seqNr - 1]:idxHigh[seqNr - 1]] / nano, y[idxLow[seqNr - 1]:idxHigh[seqNr - 1]]- mean_all)
 # plt.plot(x[idxLow[2]:idxHigh[2]] / nano, y[idxLow[2]:idxHigh[2]]- mean_all)
  


  #print("-----------------------")
  #print("mean all:",mean_all,"sigma all:", std_all)
  #print("array of mean of segments:",mean,"array of sigma of segments",std) 
  #print("-----------------------")
  #appending area of segments to A for every trace
  #length of A is number of sequnces * number of traces 
  #print(y[boolArr[0]] - mean[seqNr-1])
  A = np.append(A,area_bool)




#peaks = find_peaks_cwt(A)
#print(peaks)
#histogram of all sequences of all traces       
n,bins,patches=plt.hist(A,bins=1000,histtype='step', color = 'k')
#print(n,bins)
plt.grid(True)

plt.xlim(0.,10.0)

peaks = find_peaks(n,height=200, width=2, distance=5)
print( bins[peaks[0]])


plt.title("histogram of pulse energy")
plt.xlabel("energy [nVs]")
plt.ylabel("counts N")
#determination of mean and sigma of data in histogram
MEAN, STD = norm.fit(A)
#print("--------------------------")
#print("Entries:" , len(A))
#print("--------------------------")
#print("Mean of pulse energy :", MEAN ,"Sigma of pulse energy:", STD)
#print("--------------------------")

plt.show()
