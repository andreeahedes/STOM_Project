import STOM_higgs_tools as stom
import matplotlib.pyplot as plt 
import numpy as np 
import scipy.optimize as spo 
from scipy.stats import chisquare

from scipy.special import gamma 
def middles(x):
    return (x[:(len(x)-1)]+x[1:])/2

def exp_fit(data, lmda, a): 
    return a*np.exp(-data/lmda)

def getchisquared(binedges, binheights, a, lamda):
    binmiddle = (np.array(binedges[:(len(binedges) - 1)]) + np.array(binedges[1:])) / 2
    expected = exp_fit(binmiddle, lamda, a)
    chinom = sum((binheights-expected)**2/expected)
    return chinom


def fitting(data):#function for exponential fitting
    initial=[30,40000]
    bin_heights, bin_edges=np.histogram(data, range=(104, 155), bins=30 )
    bin_middles=(bin_edges[:(len(bin_edges)-1)]+bin_edges[1:])/2
    fit, fit_cov= spo.curve_fit(exp_fit, bin_middles, bin_heights, initial, maxfev=10**6)
    #plt.plot(bin_middles, exp_fit(bin_middles,*fit) )
    #plt.show()
    chi= getchisquared(bin_edges, bin_heights, fit[1], fit[0])
    return chi

data=stom.generate_data(400, A=0.1)
c = fitting(data)
    
    
chi2=[]
for i in range(0,500):
     back=stom.generate_background(100000, 30)
     ch=fitting(back)
     chi2.append(ch)
    
bin_heights2, bin_edges2, patches2 = plt.hist(chi2, bins=30, density=True)


def chi_fit(x,k):
    f= 1/(2**(k/2)*gamma(k/2))*x**(k/2 -1)*np.exp(-x/2)
    return f 
def areas(bin_heights, bin_edges):
    area=sum(bin_heights*(-bin_edges[:(len(bin_edges)-1)]+bin_edges[1:]))
    return area

middles2=middles(bin_edges2)
Area2=areas(bin_heights2, bin_edges2)
'''
def normalisation(heights, edges):
    width=edges[1]-edges[0]
    length=edges[len(edges)-1]-edges[0]
    h=np.zeros_like(heights)
    for i in range(0, 30):
        h[i]=heights[i]*length/width
    return h
'''    
initial=[30]
#bin_norm1=normalisation(bin_heights1, bin_edges1)
bin_norm= bin_heights2/Area2
fit2, cov_fit2=spo.curve_fit(chi_fit, middles2 ,bin_norm, initial, maxfev=10**6) 
plt.plot(middles2, chi_fit(middles2,fit2[0]))
print(fit2)
print(Area2)
plt.show()


t = np.linspace(0, 60, 1000)
fit_curve= chi_fit(t, fit2)

plt.plot(t, fit_curve)
plt.show()
#%%
from scipy.stats import chi2

print(1-chi2.cdf(c, fit2))
#%%
""" 
signal estimation thingy 

"""

def signal_back(x,A, lamb, mu, sig, signal_amp):
    return A*np.exp(-x/lamb) + signal_gaus(x, mu, sig, signal_amp)

def fitting(data):#function for exponential fitting
    initial=[30,40000,125, 1.5, 700]
    bin_heights, bin_edges=np.histogram(data, range=(104, 155), bins=30 )
    bin_middles=(bin_edges[:(len(bin_edges)-1)]+bin_edges[1:])/2
    fit, fit_cov= spo.curve_fit(signal_back, bin_middles, bin_heights, initial, maxfev=10**6)
    #plt.plot(bin_middles, exp_fit(bin_middles,*fit) )
    #plt.show()
    chi= getchisquared(bin_edges, bin_heights, fit[1], fit[0])
    return chi, fit
data=stom.generate_data(400)
fit, chi= fitting(data)
bin_heights, bin_edges=plt.hist(data, range=(104, 155), bins=30 )
plt.plot(middles(bin_edges), signal_back(middles, *fit))
