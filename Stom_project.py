# -*- coding: utf-8 -*-
"""
Created on Mon May 31 19:51:11 2021

@author: ah920
"""

import STOM_higgs_tools as stom
import matplotlib.pyplot as plt 
import numpy as np 
import scipy.optimize as spo 
from scipy.stats import chi2_contingency

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

def getpvalue(table):
    stat, p, dof, expected = chi2_contingency(table)
    return p

# print(fitting(vals))#use function above to get both fitting parameters and profile 
# #for both data and background
# print(fitting(background))
chi1=[]
for i in range(0,101):
    data=stom.generate_data(1000)
    c = fitting(data)
    chi1.append(c)
    
chi2=[]
for i in range(0,101):
    back=stom.generate_background(1000000, 30)
    ch=fitting(data)
    chi2.append(ch)
    
bin_heigths1, bin_edges1=np.histogram(chi1, bins=30)
bin_heights2, bin_edges2=np.histogram(chi2, bins=30)

print(getpvalue([bin_heigths1,bin_heights2]))
   
