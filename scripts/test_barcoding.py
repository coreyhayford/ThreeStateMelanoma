from pysb import *
from pysb.util import *
from pysb.macros import *
from pysb.bng import *
from pysb.integrate import odesolve
import pylab as pl
from numpy import linspace
from sympy import sympify
from scipy.stats import multivariate_normal
import numpy as np
import matplotlib.pyplot as plt

Model()

def declare_monomers():
    [Monomer("Barcode_%d" %i) for i in range(1,101)]
    
def declare_parameters():
    [Parameter("Barcode_%d_0" %i, 1) for i in range(1, 101)]
    
    Parameter("k_div_noDrug", 0.011)
    Parameter("k_death_noDrug", 0.01)
    
    alias_model_components()

def declare_ICs():
    [Initial("Barcode_%d", "Barcode_%d_0" %i) for i in range(1, 101)]
    
def declare_observables():
#     Observable("Obs_1", Barcode_1)
#     Observable("Obs_2", Barcode_2)
#     Observable("Obs_3", Barcode_3)
#     Observable("Obs_4", Barcode_4)
#     Observable("Obs_5", Barcode_5)
#     Observable("Obs_6", Barcode_6)
#     Observable("Obs_7", Barcode_7)
#     Observable("Obs_8", Barcode_8)
#     Observable("Obs_9", Barcode_9)
#     Observable("Obs_10", Barcode_10)
    #Observable("Obs_All", A()+B()+C()+D()+E()+F()+G()+H()+I()+J()
    
    [Observable("Obs_%d", "Barcode_%d" %i) for i in range(1,101)]

    Observable("Obs_All", np.sum(a = np.array("Barcode_%d" %i))) for i in range(1,101)

def declare_rules():
    [Rule("div_%d", "Barcode_%d" >> "Barcode_%d" + "Barcode_%d" %i, k_div_noDrug)]

    [Rule("death_%d", "Barcode_%d" %i >> None, k_death_noDrug)]
    


plt.figure("Pre-drug Stochastic Simulations")
####
for i in range(1,2):
    
    Model()  
    declare_monomers()
    declare_parameters()
    declare_ICs()
    declare_observables()
    declare_rules()

print(model.monomers)
print(model.parameters)
print(model.initial_conditions)
print(model.observables)
print(model.rules)

quit()
#     t = linspace(0, 200, 200)
#     
#     #y = odesolve(model, t, verbose=True)
#     y = run_ssa(model,t_end = t[-1], n_steps = len(t)-1, verbose=True)
#     
#     
#     #plt.subplot(3,3,i)
#     plt.plot(t, y["Obs_A"], 'r-', lw=2)
#     plt.plot(t, y["Obs_B"], 'b-', lw=2)
#     plt.plot(t, y["Obs_C"], 'g-', lw=2)
#     plt.plot(t, y["Obs_D"], 'm-', lw=2)
#     plt.plot(t, y["Obs_E"], 'c-', lw=2)
#     plt.plot(t, y["Obs_F"], 'r:', lw=2)
#     plt.plot(t, y["Obs_G"], 'b:', lw=2)
#     plt.plot(t, y["Obs_H"], 'g:', lw=2)
#     plt.plot(t, y["Obs_I"], 'm:', lw=2)
#     plt.plot(t, y["Obs_J"], 'c:', lw=2)
#     plt.plot(t, y["Obs_All"], 'k-', lw=4)
#     plt.xlabel("Time")
#     plt.ylabel("Number of Barcodes")
#     plt.title("Replicate " + str(i))
# 
# plt.show()
# 
# print(y["Obs_A"][199])
# print(y["Obs_B"][199])
# print(y["Obs_C"][199])
# print(y["Obs_D"][199])
# print(y["Obs_E"][199])
# print(y["Obs_F"][199])
# print(y["Obs_G"][199])
# print(y["Obs_H"][199])
# print(y["Obs_I"][199])
# print(y["Obs_J"][199])
# 
# print(y["Obs_A"][199]/y["Obs_All"][199])
# print(y["Obs_B"][199]/y["Obs_All"][199])
# print(y["Obs_C"][199]/y["Obs_All"][199])
# print(y["Obs_D"][199]/y["Obs_All"][199])
# print(y["Obs_E"][199]/y["Obs_All"][199])
# print(y["Obs_F"][199]/y["Obs_All"][199])
# print(y["Obs_G"][199]/y["Obs_All"][199])
# print(y["Obs_H"][199]/y["Obs_All"][199])
# print(y["Obs_I"][199]/y["Obs_All"][199])
# print(y["Obs_J"][199]/y["Obs_All"][199])
# 
# # These will be ICs for Post Drug...
