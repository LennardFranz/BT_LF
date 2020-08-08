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
from scipy.signal import find_peaks, peak_widths, peak_prominences 

# --------------------------------
# constants to be set by the user
# --------------------------------
channelNr = 1

 
StartTrace = 0
StopTrace =20



#inName = "--testpulse--"
#inName = "--testpulse--autosave--"
inName = "--testpulse--autosave--segments--"

#inDir = "/volumes/lennad/segment/segment-CeBr-cs137/"
#inDir = "/volumes/lennad/Lennard Franz/th-zylinder/"
inDir = "/volumes/lennad/Lennard Franz/cs137-longzylinder-shaper/"
#inDir = "/Users/Lennard/desktop/osci/lennard/segment-flat-CeBr-CS137/"
#inDir = "/Users/Lennard/desktop/osci/lennard/cs137-longzylinder-shaper/"
#inDir = "/Users/Lennard/destop/osci/lennard/cebr-big-4ch-cs137/"
#inDir = "/Users/Lennard/desktop/osci/lennard/segment-CeBr-ra226/"
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
  tot = np.zeros(header['SUBARRAY_COUNT'])
  A = np.zeros(header['SUBARRAY_COUNT'])

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
    
    A[seqNr-1]= np.max(y[idxLow[seqNr - 1]:idxHigh[seqNr - 1]])
    
    area_bool[seqNr-1] = np.sum((y[idxLow[seqNr - 1]:idxHigh[seqNr - 1]]*(header['HORIZ_INTERVAL']/nano)),where = boolArr) - mean[seqNr-1]*np.sum(boolArr)*header['HORIZ_INTERVAL']/nano
    tot[seqNr-1] = (np.sum(boolArr)-1)*header['HORIZ_INTERVAL']/nano

    #print(area_bool[seqNr - 1],tot[seqNr - 1] )
#plt.plot(tot,area_bool, 'bo', ms=1) 


#print(A)

#E_real = np.array([569, 662,1173,1332,1063])
#E_fit = np.array([6.485, 7.2341,10.7222,11.7675,10.2165])


E_real = np.array ([662,1173,1332,569,1063])
E_fit = np.array ([0.7535521180502402,1.1657715187599953,1.273437535904615,0.701803960984309,1.1086378999695787])


from scipy import stats 
slope, intercept, r_value, p_value, std_err = stats.linregress(E_fit, E_real)

print(slope)
print(intercept)

E = intercept+slope*area_bool

print(E)


plt.plot(tot,area_bool,'bo',ms=1.5,color='k', label='Daten' )



#plt.hist2d(tot,area_bool, bins= 50)
#plt.yscale('log')
#plt.xscale('log')

from scipy.optimize import curve_fit

def fit(x,a,b):
    return a*np.exp(b*x)

popt, pcov = curve_fit(fit, tot ,area_bool, bounds= ([0,0],[1000,1]),)

perr = np.sqrt(np.diag(pcov))


print(popt)

print(perr)


x_fit = np.linspace(50,106,1000)


plt.plot(x_fit, fit(x_fit, *popt), label='exp. Fit', color='r')


plt.text(52,0.75, '$A(T)=a \cdot e^{(b \cdot x)}$')
plt.text(52,0.70, '$a= 0.053 \ {nVs}$')
plt.text(52,0.65, '$b= 0.027 \ ns^{-1} $')


plt.ylim(0.15,1)

plt.legend()
plt.grid()
plt.xlabel('ToT (ns)')
plt.ylabel('Pulsfläche A (nVs) ')


#import tikzplotlib
#tikzplotlib.save("area_tot.tex")



plt.show()
