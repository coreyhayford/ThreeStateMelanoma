from pysb import *
from pysb.util import *
from pysb.macros import *
from pysb.bng import *
from pysb.core import *
from pysb.integrate import odesolve
from pylab import *
import pylab as pl
import numpy as np
from numpy import linspace
from sympy import sympify
from scipy import constants 
import matplotlib.pyplot as plt


 
# Import Model
import barcoding_SCFD_states

plt.figure("Barcoding Stochastic Simulations")
####
for i in range(1,2):
    
    Model()  
        
    barcoding_SCFD_states.declare_monomers()
    barcoding_SCFD_states.declare_parameters()
    barcoding_SCFD_states.declare_initial_conditions()
    barcoding_SCFD_states.declare_observables()
    barcoding_SCFD_states.declare_functions()
    
    
    for m in model.monomers:
        print m
         
    for p in model.parameters:
        print p
         
    for ic in model.initial_conditions:
        print ic
     
    for obs in model.observables:
        print obs
     
    for rules in model.rules:
        print rules
        
    for exp in model.expressions:
        print exp
        
    
    t = linspace(0, 200, 200)
    
    #y = odesolve(model, t, verbose=True)
    y = run_ssa(model,t_end = t[-1], n_steps = len(t)-1, verbose=True)
    
    
    plt.subplot(3,3,i)
    plt.plot(t, np.log2(y["Obs_A"]/y["Obs_A"][0]), 'r-', lw=2)
    plt.plot(t, np.log2(y["Obs_B"]/y["Obs_B"][0]), 'b-', lw=2)
    plt.plot(t, np.log2(y["Obs_C"]/y["Obs_C"][0]), 'g-', lw=2)
    plt.plot(t, np.log2(y["Obs_D"]/y["Obs_D"][0]), 'm-', lw=2)
    plt.plot(t, np.log2(y["Obs_E"]/y["Obs_E"][0]), 'c-', lw=2)
    plt.plot(t, np.log2(y["Obs_F"]/y["Obs_F"][0]), 'r:', lw=2)
    plt.plot(t, np.log2(y["Obs_G"]/y["Obs_G"][0]), 'b:', lw=2)
    plt.plot(t, np.log2(y["Obs_H"]/y["Obs_H"][0]), 'g:', lw=2)
    plt.plot(t, np.log2(y["Obs_I"]/y["Obs_I"][0]), 'm:', lw=2)
    plt.plot(t, np.log2(y["Obs_J"]/y["Obs_J"][0]), 'c:', lw=2)
    plt.plot(t, np.log2(y["Obs_All"]/y["Obs_All"][0]), 'k-', lw=4)
    plt.xlabel("Time")
    plt.ylabel("Population Doublings")
    plt.title("Replicate " + str(i))

plt.show()