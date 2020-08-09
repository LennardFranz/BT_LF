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
StopTrace = 20


#inName = "--testpulse--"
#inName = "--testpulse--autosave--"
inName = "--testpulse--autosave--segments--"


inDir = "/Users/Lennard/Bachelor_Thesis_LF/BA_Lennard_Franz/cs137-longzylinder/"

#inDir = "/Users/Lennard/Bachelor_Thesis_LF/BA_Lennard_Franz/ra226-longzylinder/"

#inDir = "/Users/Lennard/desktop/osci/lennard/cebr-ra226/"
inDir = "/volumes/lennad/osci-final/co60-zylinder/"
#inDir = "/volumes/lennad/segment/segment-CeBr-am241/"
#inDir = "/Users/Lennard/desktop/osci/lennard/cs137-longzylinder/"
#inDir = "/Users/Lennard/desktop/autosave-segments/"
#inDir = "/Users/Lennard/desktop/osci/lennard/cebr-big-4ch-cs137/"
#inDir = "/Users/Lennard/desktop/osci/lennard/gag-cs137/"
#inDir = "/Users/Lennard/desktop/osci/lennard/segment-CeBr-cs137/"

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
T = np.empty([0])
AA = np.empty([0])
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
  tot1 = np.zeros(header['SUBARRAY_COUNT'])
  Ampl = np.zeros(header['SUBARRAY_COUNT'])
  
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
    boolArr = np.where((y[idxLow[seqNr - 1]:idxHigh[seqNr - 1]]- mean[seqNr - 1]) > (5*std[seqNr -1]))
     
    boolArr1 = (y[idxLow[seqNr - 1]:idxHigh[seqNr - 1]] - mean[seqNr-1]) > (3*std[seqNr-1])
    #print(np.sum(boolArr1))
    #print(len(boolArr))
    
    #print(boolArr)
    #sum of all values which meet condition and scaling
    #substracting of area belonging to background
    area_bool[seqNr-1] = np.sum(y[boolArr[0]])/(nano/header['HORIZ_INTERVAL']) - (mean[seqNr-1])*(len(boolArr[0])*header['HORIZ_INTERVAL']/nano)
    
    tot[seqNr-1] = (len(boolArr[0])-1)*header['HORIZ_INTERVAL']/nano
    #tot1[seqNr-1] = (np.sum(boolArr1)-1)*header['HORIZ_INTERVAL']/nano
    #print(tot[seqNr-1],tot1[seqNr-1])
    #print(tot[seqNr-1]-tot1[seqNr-1])
    #Ampl[seqNr-1]= np.max(y[idxLow[seqNr - 1]:idxHigh[seqNr - 1]])
    #plt.plot(x[idxLow[seqNr - 1]:idxHigh[seqNr - 1]] / nano, y[idxLow[seqNr - 1]:idxHigh[seqNr - 1]]- mean_all)
    # plt.plot(x[idxLow[2]:idxHigh[2]] / nano, y[idxLow[2]:idxHigh[2]]- mean_all)



  #print("--------------print("mean all:",mean_all,"sigma all:", std_all)
  #print("array of mean of segments:",mean,"array of sigma of segments",std)
  #print("-----------------------")
  #appending area of segments to A for every trace
  #length of A is number of sequnces * number of traces
  
  #AA = np.append(AA,Ampl)
  
  A = np.append(A,area_bool)

  T = np.append(T,tot)


#print(T)

#mit shaper
a = 20.88665687
b = 0.04458985



E = a*np.exp(b*T)

print(E)




#peaks = find_peaks_cwt(A)
#print(peaks)
#histogram of all sequences of all traces

print(std_all, 6*std_all)


n,bins,patches=plt.hist(T,bins=200,  histtype= 'step', color='k')


#n,bins,patches=plt.hist(A,bins=500)
#print(n,bins)
plt.grid(True)

#plt.xscale('log')
#plt.yscale('log')

#peaks = find_peaks(n,width=7)
#print( bins[peaks[0]])

#plt.title("histogram of pulse energy")
plt.title('Histogramm der ToT')

#plt.xlabel("energy [nVs]")
plt.xlabel('ToT (ns)')

#plt.xlim(50, 55)


bin_width= (bins[1]-bins[0])

#plt.xlim(0.,15.0)

peaks,properties = find_peaks(n,width=10,height=250, rel_height=0.8)
print("bins of peaks:", peaks)
print("X Values of peaks:", bins[peaks] )
print(properties)

#print(bins[peaks])
#print(bins[np.int(peaks[6] - 1.5* properties["widths"][6]) :np.int( peaks[6]+ 1.5* (properties["widths"][6])) ])
#print(peaks[2] - 1.5* properties["widths"][2])
#print(np.int(peaks[2] - 1.5* properties["widths"][2]))

#print(n,bins)

#peakNr = 6
#print (bins[np.int(peaks[peakNr] - 1.5 * properties["widths"][peakNr]) :np.int( peaks[peakNr]+ 1.5 * (properties["widths"][peakNr])) ])
#print (n[np.int(peaks[peakNr] - 1.5* properties["widths"][peakNr]) :np.int( peaks[peakNr]+ 1.5* (properties["widths"][peakNr]))  ] )


from scipy.optimize import curve_fit
def gauss (x,a,sigma,mu):
  #return (1/np.sqrt(2*np.pi*sigma**2))* np.exp(- (x-mu)**2 / 2*sigma**2)
  return a*np.exp(-(x-mu)**2/(2*sigma**2))


bins = np.delete(bins,[len(bins)-1])


for peakNr in range(0, len(peaks) ):
  
  fit_width = 0.1
  x_values = bins[np.int(peaks[peakNr] - fit_width *  properties["widths"][peakNr]) :np.int( peaks[peakNr]+ (fit_width+0.7)* (properties["widths"][peakNr])) ]
  y_values = n[np.int(peaks[peakNr] - fit_width * properties["widths"][peakNr]) :np.int( peaks[peakNr]+ (fit_width+0.7) * (properties["widths"][peakNr]))  ]

  
  #print(x_values)
  #print(y_values)

  a_bound = [0, 1.1*np.max(y_values)]
  sigma_bound = [0, 100]
  mu_bound = [np.min(x_values), np.max(x_values)]

  #print(a_bound)
  #print(sigma_bound)
  #print(mu_bound)

  #popt, pcov = curve_fit(gauss, bins[np.int(peaks[peakNr] - 1.5 * properties["widths"][peakNr]) :np.int( peaks[peakNr]+ 1.5 * (properties["widths"][peakNr])) ] , n[np.int(peaks[peakNr] - 1.5* properties["widths"][peakNr]) :np.int( peaks[peakNr]+ 1.5* (properties["widths"][peakNr]))  ] )
  if peakNr >= 0:
    popt, pcov = curve_fit(gauss, x_values , y_values, bounds= ([a_bound[0],sigma_bound[0],mu_bound[0]], [a_bound[1],sigma_bound[1],mu_bound[1]] ) )
  
    #plt.plot(bins, gauss(bins, *popt), color='g')
    #plt.plot((popt[2], popt[2]), (0, 1350), 'r-')

    #statistic error of curve fit
    perr = np.sqrt(np.diag(pcov))

    #FWHM
    FWHM = 2* np.sqrt(2*np.log(2)) * popt[1]
    delta_FWHM = 2* np.sqrt(2*np.log(2)) * perr[1]


    print("Peak Nr:" , peakNr)
    print("Energy:" , popt[2], "statistic error of energy:" , perr[2], "FWHM:", FWHM, "statistic error of FWHM", delta_FWHM) 






plt.xlim(0,140)
plt.ylim(0,400)

plt.ylabel('Anzahl N')
#determination of mean and sigma of data in histogram
#MEAN, STD = norm.fit(A)
print("--------------------------")
print("Entries:" , len(A))
print("--------------------------")
#print("Mean of pulse energy :", MEAN ,"Sigma of pulse energy:", STD)
#print("--------------------------")



#import tikzplotlib
#tikzplotlib.save("mytikz.tex")



plt.show()
                                                                                                                                                                                                                                                                                                                                                       
