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


def declare_drug_parameters():
    
    Parameter("k_divide_A", 0.25)
    Parameter("k_divide_B", 0.25)
    Parameter("k_divide_C", 0.25)
    Parameter("k_divide_D", 0.25)
    Parameter("k_divide_E", 0.25)
    Parameter("k_divide_F", 0.25)
    Parameter("k_divide_G", 0.25)
    Parameter("k_divide_H", 0.25)
    Parameter("k_divide_I", 0.25)
    Parameter("k_divide_J", 0.25)
    
    Parameter("k_death_A", -bins_avg_array[0]+0.25)
    Parameter("k_death_B", -bins_avg_array[1]+0.25)
    Parameter("k_death_C", -bins_avg_array[2]+0.25)
    Parameter("k_death_D", -bins_avg_array[3]+0.25)
    Parameter("k_death_E", -bins_avg_array[4]+0.25)
    Parameter("k_death_F", -bins_avg_array[5]+0.25)
    Parameter("k_death_G", -bins_avg_array[6]+0.25)
    Parameter("k_death_H", -bins_avg_array[7]+0.25)
    Parameter("k_death_I", -bins_avg_array[8]+0.25)
    Parameter("k_death_J", -bins_avg_array[9]+0.25)
    
def declare_drug_rules():
    Rule('divide_A', A() >> A() + A(), k_divide_A)
    Rule('divide_B', B() >> B() + B(), k_divide_B)
    Rule('divide_C', C() >> C() + C(), k_divide_C)
    Rule('divide_D', D() >> D() + D(), k_divide_D)
    Rule('divide_E', E() >> E() + E(), k_divide_E)
    Rule('divide_F', F() >> F() + F(), k_divide_F)
    Rule('divide_G', G() >> G() + G(), k_divide_G)
    Rule('divide_H', H() >> H() + H(), k_divide_H)
    Rule('divide_I', I() >> I() + I(), k_divide_I)
    Rule('divide_J', J() >> J() + J(), k_divide_J)
    
    Rule('die_A', A() >> None, k_death_A)
    Rule('die_B', B() >> None, k_death_B)
    Rule('die_C', C() >> None, k_death_C)
    Rule('die_D', D() >> None, k_death_D)
    Rule('die_E', E() >> None, k_death_E)
    Rule('die_F', F() >> None, k_death_F)
    Rule('die_G', G() >> None, k_death_G)
    Rule('die_H', H() >> None, k_death_H)
    Rule('die_I', I() >> None, k_death_I)
    Rule('die_J', J() >> None, k_death_J)
    
# plt.figure("Pre-drug Stochastic Simulations")
####

def plot(t, t_start=0, label=False):
    plt.plot(t+t_start, y["Obs_A"], 'c-', lw=2, 
             label= ('A DIP = %g' %(k_divide_A.value-k_death_A.value)) if label else None)
    plt.plot(t+t_start, y["Obs_B"], 'b-', lw=2, 
             label= ('B DIP = %g' %(k_divide_B.value-k_death_B.value)) if label else None)
    plt.plot(t+t_start, y["Obs_C"], 'g-', lw=2, 
             label= ('C DIP = %g' %(k_divide_C.value-k_death_C.value)) if label else None)
    plt.plot(t+t_start, y["Obs_D"], 'r-', lw=2, 
             label= ('D DIP = %g' %(k_divide_D.value-k_death_D.value)) if label else None)
    plt.plot(t+t_start, y["Obs_E"], 'm-', lw=2, 
             label= ('E DIP = %g' %(k_divide_E.value-k_death_E.value)) if label else None)
    plt.plot(t+t_start, y["Obs_F"], 'c:', lw=2, 
             label= ('F DIP = %g' %(k_divide_F.value-k_death_F.value)) if label else None)
    plt.plot(t+t_start, y["Obs_G"], 'b:', lw=2, 
             label= ('G DIP = %g' %(k_divide_G.value-k_death_G.value)) if label else None)
    plt.plot(t+t_start, y["Obs_H"], 'g:', lw=2, 
             label= ('H DIP = %g' %(k_divide_H.value-k_death_H.value)) if label else None)
    plt.plot(t+t_start, y["Obs_I"], 'r:', lw=2, 
             label= ('I DIP = %g' %(k_divide_I.value-k_death_I.value)) if label else None)
    plt.plot(t+t_start, y["Obs_J"], 'm:', lw=2, 
             label= ('J DIP = %g' %(k_divide_J.value-k_death_J.value)) if label else None)
    
#     plt.plot(t+t_start, y["Obs_All"], 'k-', lw=4, 
#              label = ('All') if label else None)

    plt.xlabel("Time")
    plt.ylabel("Cell Count")
    plt.legend(loc = 0) #, fontsize = 6)
#     plt.ylim([0,1000])


for i in range(1,10):
    
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
    
    plt.subplot(3,3,i)
    plot(t, label = False)
    
    plt.axvline(x=t[-1], color='red', ls='--', lw=4)
    
    dips = np.random.normal(-0.033, 0.01, 100)
# print(dips)
    count, bins, ignored = plt.hist(dips, 10, normed=True)
    # plt.show()
    # print(count)
    # print(bins)
    
    # print(sum(count))
    
    count_normalized = count/sum(count)
    # print(count_normalized)
    # print(sum(count_normalized))
    
    bins_avg = []
    for i in range(1,len(bins+1)):
    #     print(i)
    #     print(a[i])
    #     print(a[i-1])
    #     print((bins[i]+bins[i-1])/2.)
        bins_avg.append((bins[i]+bins[i-1])/2.)
    bins_avg_array = np.array(bins_avg)
        # print(bins_avg)
    # print(bins_avg_array)
    
    
    ## Number of cells in each state - NOT necessary is explicitly modeled
    picks = np.random.multinomial(10, count_normalized)
    # print(picks)

# These will be ICs for Post Drug...

    A_0.value = y["Obs_A"][-1]
    B_0.value = y["Obs_B"][-1]
    C_0.value = y["Obs_C"][-1]
    D_0.value = y["Obs_D"][-1]
    E_0.value = y["Obs_E"][-1]
    F_0.value = y["Obs_F"][-1]
    G_0.value = y["Obs_G"][-1]
    H_0.value = y["Obs_H"][-1]
    I_0.value = y["Obs_I"][-1]
    J_0.value = y["Obs_J"][-1]



    declare_drug_parameters()
    declare_newparameters()
    declare_drug_rules()

    y = run_ssa(model,t_end = t[-1], n_steps = len(t)-1, verbose=True)
    plot(t, t_start=t[-1], label=False)

#     print(A_0.value, B_0.value, C_0.value, D_0.value, E_0.value, F_0.value, G_0.value, H_0.value, I_0.value, J_0.value)

plt.tight_layout()
plt.show()
## DONE: established initial conditions as a function of gillespie sims
## DONE: Pick DIP rates from multinomial distribution
## DONE: Plot before and after drug simulated results
## DONE: Apply new plot function to stitch together results
## DONE: Put in line for drug introduction
## DONE: Stitch together subplots and make into figure - made in a loop


## GOAL: Incorporate James' GPU simulator with only 1 cell and run 1000s of sims
## GOAL: Plot the clonal distributions like Hata
## GOAL: Make it for 1 million barcoded cells rather than 10
## GOAL: Play with different DIP rates/distributions to model diversification phenotype
## GOAL: Only pick time points similar to experiment??
## GOAL: Include state transitions in the model

## GOAL: Add more barcodes and plot like Hata
