import numpy as np 
import matplotlib.pyplot as plt 
from scipy import stats 

#zylinder ohne  shaper
#E_real = np.array([80,186,242,295,352,32,662,1173,1332,1063,569])
#E_fit = np.array([0.97844,2.4133,3.137,3.7597,4.3779,0.3457,7.2341,10.7222,11.7675,10.2165,6.485])

#E_real = np.array([569, 662,1173,1332,1063])
#E_fit = np.array([6.472, 7.19,10.703,11.733,10.19])
#FWHM = np.array([0.51, 0.44, 0.87, 0.64, 0.48]) 

# zylinder mit shaper
E_real = np.array([662,1173,1332,569,1063])
E_fit = np.array([0.753,1.163,1.271,0.7018,1.1080])
 
FWHM_E = np.array([0.052, 0.069,0.0645, 0.065, 0.055])


#mit shaper ToT
#ToT_real = np.array([662, 1063,1252])
#ToT_fit = np.array([91.00116606316232,95.21278835637344,107.05266288124014])

#FWHM_T = np.array([9.798621093190247,15.77351053281185,22.63164384355587])
#x_err = np.array([0,0,79.5])



#tot mit shaper final

ToT_real = np.array([662,727,1069,1252])
ToT_fit = np.array([77.29,80.00,88.00,91.91]) 

FWHM_T = np.array([4.26,5,4,11]) 

#print(ToT_real)
#print(ToT_fit)

#GSI ToT

ToT_GSI = np.array([18.53,18.91,20.73,21.19,22.76])
E_real_GSI = np.array([569,662,1063,1173,1332])
#E_fit_GSI = np.array([0.701803960984309,0.7535521180502402,1.1657715187599953,1.273437535904615])
#ToT_GSI = np.array([18.53,18.95,21.2,22.75])

FWHM_GSI = np.array([1.88,1.68,2.81,3.31,3.46])



#bi207 1069keV 20,2

#plt.plot(E_fit,E_real,ls='None',marker='.')
#plt.errorbar(E_fit,E_real, xerr= 0.5*FWHM_E , fmt='.',capsize = 3,color='k', label='zugeordnete Energie') 

#plt.plot(ToT_fit,ToT_real, ls='None', marker='.')
#plt.errorbar(ToT_fit,ToT_real, xerr=0.5*FWHM_T,fmt='.',capsize = 3, label='zugeordnete Energie', color ='k') 

#plt.plot(test_tot,test_E, ls='None', marker='.')

plt.errorbar(ToT_GSI,E_real_GSI,xerr=0.5*FWHM_GSI,fmt='.',capsize = 3,color='k', label='zugeordnete Energie') 


#Energieauflösung

#E_fit_input = E_fit
#E_real_input = E_real
#FWHM_input = FWHM

#delta_E_fit = np.zeros([len(E_fit_input)])
#delta_E_real = np.zeros(len(E_fit_input))
#for i in range(0, len(E_fit_input), 1):
 #   print(i)
  #  delta_E_fit[i] = FWHM_input[i]/E_fit_input[i] 
   # delta_E_real[i] = delta_E_fit[i]*E_real_input[i]

#plt.plot(E_real, delta_E_fit, ls='None', marker = '.')


#print('%:',delta_E_fit)
#print('absolut:', delta_E_real)
#print(delta_E_fit*100)

from scipy.optimize import curve_fit

def fit_exp(x,a,b):
    return a*np.exp(b*x)

def fit_log(x,a,b):
    return a*np.log(x)+b

popt, pcov = curve_fit(fit_exp, ToT_GSI , E_real_GSI, bounds=([15,0],[35,1]) ,method ='trf')

#popt,pcov = curve_fit(fit, E_fit, E_real , bounds=([0,0,0],[1000,10,100]), method='trf'  )

perr = np.sqrt(np.diag(pcov))

print(popt)
print(perr)



#slope, intercept, r_value, p_value, std_err = stats.linregress(E_fit, E_real)

slope, intercept, r_value, p_value, std_err = stats.linregress(ToT_GSI, E_real_GSI)


fit_data = np.array([slope, intercept])
#print(slope)
#print(intercept)
#print(r_value)
#print(p_value)
#print(std_err)


x_fit_E = np.linspace(17,24,10)

plt.plot(x_fit_E, intercept + slope*x_fit_E, 'r',label='linearer Fit')

x_fit = np.linspace(17,24,1000)

#plt.plot(x_fit, fit_exp(x_fit, *popt), label='exp. Fit',color='r')

plt.legend()
plt.grid()
#plt.text(470,1.2, 'Q(E) = %5.5fE + %5.5f' %tuple(fit_data))
plt.text(17,1200, '$E(T) = m \cdot T + n$')
#plt.text(17,1200, '$E(T)=a \cdot e^{(b \cdot x)}$')
plt.text(17,1125, '$m= 187.69{keV/ns} $')
plt.text(17,1050, '$n= -2873.54 keV $')



plt.title("")
plt.ylabel("Energie E (keV)")
plt.xlabel("ToT T (ns)")

import tikzplotlib

#tikzplotlib.save("calib_ToT_GSI_linear.tex")



plt.show()
