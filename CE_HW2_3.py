# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 15:12:59 2017

@author: frazeranderson1
"""

#packages
import sys
sys.path.append('10_mcs')
import math
import numpy as np
np.set_printoptions(suppress=True,
                    formatter={'all': lambda x: '%7.6f' % x})
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['font.family'] = 'serif'
import scipy.interpolate as sci
from scipy.optimize import fmin

#Import Data
import os
import pandas as pd

os.chdir("/Users/frazeranderson1/Desktop/School/Econ 647 - Applied Computational/Problem Sets/PS2/Data")

df = pd.read_csv("Data_2.csv")
df.index = pd.DatetimeIndex(df['Date'])
df = df.drop("Date", axis = 1)


#Create The inputs

t_list = np.array((0.25, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 30.0))/1.

for i in range(0,5):
	r_list = df.values[i]/100
	
	
	#Formatting the inputs
	factors = (1 + t_list * r_list)
	zero_rates = 1 / t_list * np.log(factors)
	
	r0 = r_list[0] # 0.0  # set to zero 
	
	#
	# Interpolation of Market Data
	#
	
	tck = sci.splrep(t_list, zero_rates, k=3)  # cublen(tic splines
	tn_list = np.linspace(0.0, 1.0, 24)
	ts_list = sci.splev(tn_list, tck, der=0)
	de_list = sci.splev(tn_list, tck, der=1)
	
	f = ts_list + de_list * tn_list
	
	
	
	#The CIR Forward Rate Function
	def CIR_forward_rate(opt):
	    ''' Function for forward rates in CIR85 model.
	
	    Parameters
	    ==========
	    kappa_r: float
	        mean-reversion factor
	    theta_r: float 
	        long-run mean
	    sigma_r: float
	        volatility factor
	
	    Returns
	    =======
	    forward_rate: float
	        forward rate
	    '''
	    kappa_r, theta_r, sigma_r = opt
	    t = tn_list
	    g = np.sqrt(kappa_r ** 2 + 2 * sigma_r ** 2)
	    sum1 = ((kappa_r * theta_r * (np.exp(g * t) - 1)) /
	          (2 * g + (kappa_r + g) * (np.exp(g * t) - 1)))
	    sum2 = r0 * ((4 * g ** 2 * np.exp(g * t)) /
	            (2 * g + (kappa_r + g) * (np.exp(g * t) - 1)) ** 2)
	    forward_rate = sum1 + sum2
	    return forward_rate
	
	
	#The MSE Function
	def CIR_error_function(opt):
	    ''' Error function for CIR85 model calibration. '''
	    kappa_r, theta_r, sigma_r = opt
	    if 2 * kappa_r * theta_r < sigma_r ** 2:
	        return 100
	    if kappa_r < 0 or theta_r < 0 or sigma_r < 0.001:
	        return 100
	    forward_rates = CIR_forward_rate(opt)
	    MSE = np.sum((f - forward_rates) ** 2) / len(f)
	    # print opt, MSE
	    return MSE
	
	
	#
	# Calibration Procedure
	#
	def CIR_calibration():
	    opt = fmin(CIR_error_function, [1.0, 0.02, 0.1],
	            xtol=0.00001, ftol=0.00001,
	            maxiter=300, maxfun=500)
	    return opt
	
	
	print(CIR_forward_rate(CIR_calibration()))
	
	print(CIR_calibration())