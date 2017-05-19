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
    Monomer("A")
    Monomer("B")
    Monomer("C")
    Monomer("D")
    Monomer("E")
    Monomer("F")
    Monomer("G")
    Monomer("H")
    Monomer("I")
    Monomer("J")
    

def declare_parameters():
    Parameter("A_0", 1)
    Parameter("B_0", 1)
    Parameter("C_0", 1)
    Parameter("D_0", 1)
    Parameter("E_0", 1)
    Parameter("F_0", 1)
    Parameter("G_0", 1)
    Parameter("H_0", 1)
    Parameter("I_0", 1)
    Parameter("J_0", 1)

    Parameter("k_div_noDrug", 0.05)
    Parameter("k_death_noDrug", 0.025)
    
    alias_model_components()

def declare_ICs():
    Initial(A, A_0)
    Initial(B, B_0)
    Initial(C, C_0)
    Initial(D, D_0)
    Initial(E, E_0)
    Initial(F, F_0)
    Initial(G, G_0)
    Initial(H, H_0)
    Initial(I, I_0)
    Initial(J, J_0)
    
def declare_observables():
    Observable("Obs_A", A)
    Observable("Obs_B", B)
    Observable("Obs_C", C)
    Observable("Obs_D", D)
    Observable("Obs_E", E)
    Observable("Obs_F", F)
    Observable("Obs_G", G)
    Observable("Obs_H", H)
    Observable("Obs_I", I)
    Observable("Obs_J", J)
    Observable("Obs_All", A()+B()+C()+D()+E()+F()+G()+H()+I()+J())

def declare_rules():
    Rule('div_A', A() >> A() + A(), k_div_noDrug)
    Rule('div_B', B() >> B() + B(), k_div_noDrug)
    Rule('div_C', C() >> C() + C(), k_div_noDrug)
    Rule('div_D', D() >> D() + D(), k_div_noDrug)
    Rule('div_E', E() >> E() + E(), k_div_noDrug)
    Rule('div_F', F() >> F() + F(), k_div_noDrug)
    Rule('div_G', G() >> G() + G(), k_div_noDrug)
    Rule('div_H', H() >> H() + H(), k_div_noDrug)
    Rule('div_I', I() >> I() + I(), k_div_noDrug)
    Rule('div_J', J() >> J() + J(), k_div_noDrug)
    
    Rule('death_A', A() >> None, k_death_noDrug)
    Rule('death_B', B() >> None, k_death_noDrug)
    Rule('death_C', C() >> None, k_death_noDrug)
    Rule('death_D', D() >> None, k_death_noDrug)
    Rule('death_E', E() >> None, k_death_noDrug)
    Rule('death_F', F() >> None, k_death_noDrug)
    Rule('death_G', G() >> None, k_death_noDrug)
    Rule('death_H', H() >> None, k_death_noDrug)
    Rule('death_I', I() >> None, k_death_noDrug)
    Rule('death_J', J() >> None, k_death_noDrug)
    

# plt.figure("Pre-drug Stochastic Simulations")
####

def plot(t, t_start=0, label=False):
    plt.plot(t+t_start, np.log2(y["Obs_A"]), 'c-', lw=2, 
             label= ('A DIP = %g' %(k_div_noDrug.value-k_death_noDrug.value)))
    plt.plot(t+t_start, np.log2(y["Obs_B"]), 'b-', lw=2, 
             label= ('B DIP = %g' %(k_div_noDrug.value-k_death_noDrug.value)))
    plt.plot(t+t_start, np.log2(y["Obs_C"]), 'g-', lw=2, 
             label= ('C DIP = %g' %(k_div_noDrug.value-k_death_noDrug.value)))
    plt.plot(t+t_start, np.log2(y["Obs_D"]), 'r-', lw=2, 
             label= ('D DIP = %g' %(k_div_noDrug.value-k_death_noDrug.value)))
    plt.plot(t+t_start, np.log2(y["Obs_E"]), 'm-', lw=2, 
             label= ('E DIP = %g' %(k_div_noDrug.value-k_death_noDrug.value)))
    plt.plot(t+t_start, np.log2(y["Obs_F"]), 'c:', lw=2, 
             label= ('F DIP = %g' %(k_div_noDrug.value-k_death_noDrug.value)))
    plt.plot(t+t_start, np.log2(y["Obs_G"]), 'b:', lw=2, 
             label= ('G DIP = %g' %(k_div_noDrug.value-k_death_noDrug.value)))
    plt.plot(t+t_start, np.log2(y["Obs_H"]), 'g:', lw=2, 
             label= ('H DIP = %g' %(k_div_noDrug.value-k_death_noDrug.value)))
    plt.plot(t+t_start, np.log2(y["Obs_I"]), 'r:', lw=2, 
             label= ('I DIP = %g' %(k_div_noDrug.value-k_death_noDrug.value)))
    plt.plot(t+t_start, np.log2(y["Obs_J"]), 'm:', lw=2, 
             label= ('J DIP = %g' %(k_div_noDrug.value-k_death_noDrug.value)))
    
    plt.plot(t+t_start, np.log2(y["Obs_All"]), 'k-', lw=4, 
             label = ('All'))
    plt.xlabel("Time")
    plt.ylabel("Population Doublings")
    plt.legend(loc = 0) #, fontsize = 6)

for i in range(1,2):
    
    Model()  
    declare_monomers()
    declare_parameters()
    declare_ICs()
    declare_observables()
    declare_rules()
    
#     for m in model.monomers:
#         print m
#     quit()
    
    t = linspace(0, 200, 201)
    
    #y = odesolve(model, t, verbose=True)
    y = run_ssa(model,t_end = t[-1], n_steps = len(t)-1, verbose=True)
    
    
    plot(t, label = True)
    #plt.subplot(2,2,i)
#     plt.plot(t, y["Obs_A"], 'r-', lw=2, label = "A")
#     plt.plot(t, y["Obs_B"], 'b-', lw=2, label = "B")
#     plt.plot(t, y["Obs_C"], 'g-', lw=2, label = "C")
#     plt.plot(t, y["Obs_D"], 'm-', lw=2, label = "D")
#     plt.plot(t, y["Obs_E"], 'c-', lw=2, label = "E")
#     plt.plot(t, y["Obs_F"], 'r:', lw=2, label = "F")
#     plt.plot(t, y["Obs_G"], 'b:', lw=2, label = "G")
#     plt.plot(t, y["Obs_H"], 'g:', lw=2, label = "H")
#     plt.plot(t, y["Obs_I"], 'm:', lw=2, label = "I")
#     plt.plot(t, y["Obs_J"], 'c:', lw=2, label = "J")
#     plt.plot(t, y["Obs_All"], 'k-', lw=4, label = "All")
#     plt.xlabel("Time")
#     plt.ylabel("Number of Barcodes")
#     plt.legend(loc = 0)
#     plt.title("Replicate " + str(i))

#plt.tight_layout()
# plt.show()

# print(y["Obs_A"][-1])
# print(y["Obs_B"][-1])
# print(y["Obs_C"][-1])
# print(y["Obs_D"][-1])
# print(y["Obs_E"][-1])
# print(y["Obs_F"][-1])
# print(y["Obs_G"][-1])
# print(y["Obs_H"][-1])
# print(y["Obs_I"][-1])
# print(y["Obs_J"][-1])
# 
# print(y["Obs_A"][-1]/y["Obs_All"][-1])
# print(y["Obs_B"][-1]/y["Obs_All"][-1])
# print(y["Obs_C"][-1]/y["Obs_All"][-1])
# print(y["Obs_D"][-1]/y["Obs_All"][-1])
# print(y["Obs_E"][-1]/y["Obs_All"][-1])
# print(y["Obs_F"][-1]/y["Obs_All"][-1])
# print(y["Obs_G"][-1]/y["Obs_All"][-1])
# print(y["Obs_H"][-1]/y["Obs_All"][-1])
# print(y["Obs_I"][-1]/y["Obs_All"][-1])
# print(y["Obs_J"][-1]/y["Obs_All"][-1])

# These will be ICs for Post Drug...
def declare_newparameters():
    Parameter("postDrugIC_A", y["Obs_A"][-1])
    Parameter("postDrugIC_B", y["Obs_B"][-1])
    Parameter("postDrugIC_C", y["Obs_C"][-1])
    Parameter("postDrugIC_D", y["Obs_D"][-1])
    Parameter("postDrugIC_E", y["Obs_E"][-1])
    Parameter("postDrugIC_F", y["Obs_F"][-1])
    Parameter("postDrugIC_G", y["Obs_G"][-1])
    Parameter("postDrugIC_H", y["Obs_H"][-1])
    Parameter("postDrugIC_I", y["Obs_I"][-1])
    Parameter("postDrugIC_J", y["Obs_J"][-1])
 
    Parameter("postDrugIC_All", y["Obs_A"][-1] + y["Obs_B"][-1] + y["Obs_C"][-1] + y["Obs_D"][-1] + y["Obs_E"][-1] + y["Obs_F"][-1] + y["Obs_G"][-1] + y["Obs_H"][-1] + y["Obs_I"][-1] + y["Obs_J"][-1])
    
    alias_model_components()

# declare_newparameters()    
# print(model.parameters)