import numpy as np
import matplotlib.pyplot as plt

#inDir = "GSI/bi207-1660mV.txt"
#inDir = "GSI/cs137-1660mV-neu.txt"
inDir = "GSI/co60-1660mV.txt"

#inDir = "test-frequenzgenerator-pos-pulse.txt"
#inDir = "GSI/cs137-noshaper-1400mV.txt"


data = np.loadtxt(inDir)   
x = data[:,0]
y = data[:,3]
    
  
#a = 32.63 
#b = 0.169

#x = a*np.exp(b*x)
  
print(x)
  

N_bins = 500
bin_width = int(10000/N_bins)



idxLow = np.linspace(0,9999 - (bin_width-1), N_bins, dtype = np.int)
idxHigh = np.linspace((bin_width - 1),9999, N_bins, dtype = np.int)
idxMid = np.linspace((bin_width/2 - 1),9999 - (bin_width/2) , N_bins, dtype = np.int)



bins = np.zeros(N_bins)
n = np.zeros(N_bins) 

#print(idxLow)
#print(idxMid,len(idxMid))
#print(idxHigh)

for i in range(0,N_bins,1):
    #print(i)
    bins[i]= x[idxMid[i]]
    n[i] = np.sum(y[idxLow[i]:idxHigh[i]])
    

#print(n)
#print(bins)
    
bins_new = np.zeros(10000)
n_new = np.zeros(10000)
    


for i in range(0,N_bins,1):
    bins_new[idxLow[i] : idxHigh[i]] = x[idxMid[i]]
    n_new[idxLow[i]:idxHigh[i]] = np.sum(y[idxLow[i]:idxHigh[i]])

#print(bins_new)
#print(n_new)

a = 32.63 
b = 0.169

calib = a*np.exp(b*bins)

calib_bins = calib[calib<1600]

N = n[calib<1600]

#plt.plot(x,y)

n,bins,patches=plt.hist(bins, bins, weights=n, histtype = 'step', color='k')

#n,bins,patches=plt.hist(calib_bins, bins, weights=N, histtype = 'step', color='k')

from scipy.signal import find_peaks
peaks,properties = find_peaks(n,width=3,height=600, rel_height=0.7)
print("index of peaks:", peaks)
print("X Values of peaks:", x[peaks] )
print(properties)

from scipy.optimize import curve_fit

def gauss (x,a,sigma,mu):
    return a*np.exp(-(x-mu)**2/(2*sigma**2)) 



for peakNr in range(1, len(peaks) ):
    
    #peakNr = 1 
    fit_width = 0.5
    
    x_values = bins[np.int(peaks[peakNr] - (fit_width) * properties["widths"][peakNr]) :np.int( peaks[peakNr]+ (fit_width+0.15) * (properties["widths"][peakNr])) ]
    y_values = n[np.int(peaks[peakNr] - (fit_width) * properties["widths"][peakNr]) :np.int( peaks[peakNr]+ (fit_width+0.15) * (properties["widths"][peakNr]))  ]
  
    a_bound = [0, 1.1*np.max(y_values)]
    sigma_bound = [0, 1]
    mu_bound = [np.min(x_values), np.max(x_values)]


    popt, pcov = curve_fit(gauss, x_values , y_values, bounds= ([a_bound[0],sigma_bound[0],mu_bound[0]], [a_bound[1],sigma_bound[1],mu_bound[1]] ) ) 


    #plt.plot(x, gauss(x, *popt), label='gauss Fit', color='g')
    #plt.plot((popt[2], popt[2]), (0, 1500), 'r-')
    
    #statistic error of curve fit 
    perr = np.sqrt(np.diag(pcov))
  
    #FWHM
    FWHM = 2* np.sqrt(2*np.log(2)) * popt[1]
    delta_FWHM = 2* np.sqrt(2*np.log(2)) * perr[1]

  
    print("Peak Nr:" , peakNr) 
    print("Fit Bereich:" , mu_bound)
    print("Energy:" , popt[2], "statistic error of energy:" , perr[2], "FWHM:", FWHM, "statistic error of FWHM", delta_FWHM )
  

for peakNr in range(1, len(peaks) ):
    
    #peakNr = 0 
    fit_width = 0.1
    
    x_values = bins[np.int(peaks[peakNr] - (fit_width) * properties["widths"][peakNr]) :np.int( peaks[peakNr]+ (fit_width+0.2) * (properties["widths"][peakNr])) ]
    y_values = n[np.int(peaks[peakNr] - (fit_width) * properties["widths"][peakNr]) :np.int( peaks[peakNr]+ (fit_width+0.2) * (properties["widths"][peakNr]))  ]
  
    a_bound = [0, 1.1*np.max(y_values)]
    sigma_bound = [0, 1]
    mu_bound = [np.min(x_values), np.max(x_values)]


    popt, pcov = curve_fit(gauss, x_values , y_values, bounds= ([a_bound[0],sigma_bound[0],mu_bound[0]], [a_bound[1],sigma_bound[1],mu_bound[1]] ) ) 


    plt.plot(x, gauss(x, *popt), label='gauss Fit', color='g')
    plt.plot((popt[2], popt[2]), (0, 1500), 'r-')
    
    #statistic error of curve fit 
    perr = np.sqrt(np.diag(pcov))
  
    #FWHM
    FWHM = 2* np.sqrt(2*np.log(2)) * popt[1]
    delta_FWHM = 2* np.sqrt(2*np.log(2)) * perr[1]

  
    print("Peak Nr:" , peakNr) 
    print("Fit Bereich:" , mu_bound)
    print("Energy:" , popt[2], "statistic error of energy:" , perr[2], "FWHM:", FWHM, "statistic error of FWHM", delta_FWHM )
  





xticks = np.arange(0,110,10)
plt.xticks(xticks)

plt.xlim(0,100)
#plt.legend()
plt.grid()
plt.title("Histogramm der ToT")
plt.xlabel("ToT (ns)")
plt.ylabel("Anzahl N")

plt.xlim(0,30)
#plt.ylim(0,1250)

#import tikzplotlib
#tikzplotlib.save("GSI_generator.tex")




plt.show()




# FÜr den fi der mit dem Frequenzgernator erstellten Histogramme ergibt sich eine
# FWHM von 0.25 bis 0.35 bei 1 ns ungenauigkeit des genrator 
# Trailing und falling von je  2,5ns beachten 
#
#
#
#
#
#
#
#
#








