"""
Model-data overlay from Three-State Melanoma model and mixing experiment.

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



# Data - need nl2 and time

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

plt.figure()
plt.plot(time, nl2)
plt.legend(loc=0, prop={'size': 16})
plt.xlabel("Time", fontsize=22)
plt.ylabel("Population Doublings", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.title("Mixing", fontsize=22)

plt.show()


### Function to define model and parameters ###
def eq(par,initial_cond,start_t,end_t,incr):
     #-time-grid
     t  = np.linspace(start_t, end_t,incr)
     #differential-eq-system----------------------
     def funct(y,t):
        A=y[0]
        B=y[1]
        C=y[2]
        par = k_div_A, k_div_B, k_div_C, k_death_A, k_death_B, k_death_C, k_AB, k_BA, k_CB, k_BC
        # the model ODEs - from PySB model
        f0 = -A*k_AB - A*k_death_A + A*k_div_A + B*k_BA
        f1 = A*k_AB - B*k_BA - B*k_BC - B*k_death_B + B*k_div_B + C*k_CB
        f2 = B*k_BC - C*k_CB -C*k_death_C + C*k_div_C
        return [f0, f1, f2]
     #integrate
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
#=======================================================

### Define Initial Conditions ###
A_0 = 1000
B_0 = 1000
C_0 = 1000
y_0 = [A_0, B_0, C_0]
#=======================================================

### Use eq() function defined above to get population dynamics for each of three states ###
F0,F1,F2,T= eq(rates,y_0,0,120,10)

### Plot the output (total) of the model and the data ###
plt.figure()
plt.plot(T,np.log2((F0+F1+F2)/(F0[0]+F1[0]+F2[0])),'-b',time,nl2,'go')
#plt.legend(('A','B','C','Experimental'),'upper center')
plt.xlabel('Time')
plt.ylabel('Population Doublings')
plt.title('Model_Data')

plt.show()    