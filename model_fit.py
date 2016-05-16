"""
Model-data overlay with fit score from Three-State Melanoma model and mixing experiment.

The current parameter estimation only includes data from parental replicates from 
mixing experiment.

Expectation: Need to add data from additional mixing of subclones.

Based on code from: http://adventuresinpython.blogspot.com/2012/08/fitting-differential-equation-system-to.html

Implemented by: Corey E. Hayford (Vito Quaranta Lab, Vanderbilt University)

"""

import pandas as pd
from pylab import *
from numpy import *
from numpy.random import rand
import numpy as np
from scipy.optimize import fmin
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy import integrate

### Data - need nl2 and time ###
data = pd.read_csv('/Users/Corey/Documents/QuarantaLab/mixing_all', sep = "\t")
print data

data_parental = pd.read_csv('/Users/Corey/Documents/QuarantaLab/mixing_parental', sep = "\t")
print data_parental

data_parental_treated = data_parental.loc[data_parental['conc'] == 8]
print data_parental_treated

time = data_parental_treated["Time"]
print time

nl2 = data_parental_treated["nl2"]
print nl2

cell_count = data_parental_treated["Cell.Nucleus"]
print cell_count

### Function to define model and parameters ###
def eq(par,initial_cond,start_t,end_t,incr):
     # Create time space
     t  = np.linspace(start_t, end_t,incr)
     # Setup ODE system
     def funct(y,t):
        A=y[0]
        B=y[1]
        C=y[2]
        par = k_div_A, k_div_B, k_div_C, k_death_A, k_death_B, k_death_C, k_AB, k_BA, k_CB, k_BC
        # the model equations from PySB model
        f0 = -A*k_AB - A*k_death_A + A*k_div_A + B*k_BA
        f1 = A*k_AB - B*k_BA - B*k_BC - B*k_death_B + B*k_div_B + C*k_CB
        f2 = B*k_BC - C*k_CB -C*k_death_C + C*k_div_C
        return [f0, f1, f2]
     # Integrate over the equations
     ds = integrate.odeint(funct,initial_cond,t)
     return (ds[:,0],ds[:,1],ds[:,2],t)
#=======================================================

### Define Parameters ###
# May want to just use DIP rate instead of div and death
k_div_A = 0.033
k_div_B = 0
k_div_C =  0.016
    
k_death_A = 0.066
k_death_B = 0
k_death_C = 0
    
k_AB = 0.025
k_BA = 0.00004
k_CB = 0.025
k_BC = 0.00004

rates = (k_div_A, k_div_B, k_div_C, k_death_A, k_death_B, k_death_C, k_AB, k_BA, k_CB, k_BC)


### Define Initial Conditions ###
A_0 = 1000
B_0 = 1000
C_0 = 1000
y_0 = [A_0, B_0, C_0]

### Redefine model steps ###
start_time=0.0
end_time=120.0
intervals=10
mt=np.linspace(start_time,end_time,intervals)

### Model index to compare to data ###
findindex=lambda x:np.where(mt>=x)[0][0]
mindex=map(findindex,time)
print mindex

### Score Fit of System ###
def score(params):
    #a.Get Solution to system
    F0,F1,F2,T=eq(params,y_0,start_time,end_time,intervals)
    #b.Pick of Model Points to Compare
    Zm=np.log2((F0+F1+F2)/(F0[0]+F1[0]+F2[0]))[mindex]
    #c.Score Difference between model and data points
    ss=lambda data,model:((data-model)**2).sum()
    return ss(nl2,Zm)
 
fit_score=score(rates)
print fit_score

### Generate Solution to System ###
F0,F1,F2,T=eq(rates,y_0,start_time,end_time,intervals)
Zm=np.log2((F0+F1+F2)/(F0[0]+F1[0]+F2[0]))[mindex]
Tm=T[mindex]

### Plot Solution to System with experimental data ###
plt.figure()
plt.plot(T,np.log2((F0+F1+F2)/(F0[0]+F1[0]+F2[0])),'-b', Tm, Zm, "ro" ,time,nl2,'go')
#plt.legend(('A','B','C','Experimental'),'upper center')
plt.xlabel('Time')
plt.ylabel('Population Doublings')
title='Model_Data + Fit Score: '+str(fit_score)
plt.title(title)

plt.show()