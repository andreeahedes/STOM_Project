#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 13 11:43:33 2021

@author: matekoszta
"""


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

def getpvalue(observed, expected, ddof):
   chi, p= chisquare(observed, expected, ddof)
   return p

chi1=[]
for i in range(0,1000):
    data=stom.generate_data(1000)
    c = fitting(data)
    chi1.append(c)
    
# chi2=[]
# for i in range(0,500):
#     back=stom.generate_data(1000,0)
#     ch=fitting(back)
#     chi2.append(ch)
    
bin_heights1, bin_edges1, patches1 = plt.hist(chi1, bins=30, density=True)
# bin_heights2, bin_edges2, patches2 = plt.hist(chi2, bins=30)


def chi_fit(x,k):
    f= 1/(2**(k/2)*gamma(k/2))*x**(k/2 -1)*np.exp(-x/2)
    return f 
def areas(bin_heights, bin_edges):
    area=sum(bin_heights*(bin_edges[:(len(bin_edges)-1)]+bin_edges[1:]))
    return area

middles1=middles(bin_edges1)
#middles2=middles(bin_edges2)
Area1=areas(bin_heights1,bin_edges1)
#Area2=areas(bin_heights2, bin_edges2)

initial=[29]
fit1, cov_fit1=spo.curve_fit(chi_fit, middles1, bin_heights1, initial, maxfev=10**6) 
plt.plot(middles1, chi_fit(middles1,fit1[0]))
print(fit1)






