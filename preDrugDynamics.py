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

    Parameter("k_div_noDrug", 0.031)
    Parameter("k_death_noDrug", 0.03)
    
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
    
    Rule('death_A', A() >> A() + A(), k_death_noDrug)
    Rule('death_B', B() >> B() + B(), k_death_noDrug)
    Rule('death_C', C() >> C() + C(), k_death_noDrug)
    Rule('death_D', D() >> D() + D(), k_death_noDrug)
    Rule('death_E', E() >> E() + E(), k_death_noDrug)
    Rule('death_F', F() >> F() + F(), k_death_noDrug)
    Rule('death_G', G() >> G() + G(), k_death_noDrug)
    Rule('death_H', H() >> H() + H(), k_death_noDrug)
    Rule('death_I', I() >> I() + I(), k_death_noDrug)
    Rule('death_J', J() >> J() + J(), k_death_noDrug)

plt.figure("Pre-drug Stochastic Simulations")
####
for i in range(1,2):
    
    Model()  
    declare_monomers()
    declare_parameters()
    declare_ICs()
    declare_observables()
    declare_rules()

    t = linspace(0, 200, 200)
    
    #y = odesolve(model, t, verbose=True)
    y = run_ssa(model,t_end = t[-1], n_steps = len(t)-1, verbose=True)
    
    
    #plt.subplot(2,2,i)
    plt.plot(t, y["Obs_A"], 'r-', lw=2)
    plt.plot(t, y["Obs_B"], 'b-', lw=2)
    plt.plot(t, y["Obs_C"], 'g-', lw=2)
    plt.plot(t, y["Obs_D"], 'm-', lw=2)
    plt.plot(t, y["Obs_E"], 'c-', lw=2)
    plt.plot(t, y["Obs_F"], 'r:', lw=2)
    plt.plot(t, y["Obs_G"], 'b:', lw=2)
    plt.plot(t, y["Obs_H"], 'g:', lw=2)
    plt.plot(t, y["Obs_I"], 'm:', lw=2)
    plt.plot(t, y["Obs_J"], 'c:', lw=2)
    plt.plot(t, y["Obs_All"], 'k-', lw=4)
    plt.xlabel("Time")
    plt.ylabel("Number of Barcodes")
    plt.title("Replicate " + str(i))

plt.tight_layout()
plt.show()

print(y["Obs_A"][199])
print(y["Obs_B"][199])
print(y["Obs_C"][199])
print(y["Obs_D"][199])
print(y["Obs_E"][199])
print(y["Obs_F"][199])
print(y["Obs_G"][199])
print(y["Obs_H"][199])
print(y["Obs_I"][199])
print(y["Obs_J"][199])

print(y["Obs_A"][199]/y["Obs_All"][199])
print(y["Obs_B"][199]/y["Obs_All"][199])
print(y["Obs_C"][199]/y["Obs_All"][199])
print(y["Obs_D"][199]/y["Obs_All"][199])
print(y["Obs_E"][199]/y["Obs_All"][199])
print(y["Obs_F"][199]/y["Obs_All"][199])
print(y["Obs_G"][199]/y["Obs_All"][199])
print(y["Obs_H"][199]/y["Obs_All"][199])
print(y["Obs_I"][199]/y["Obs_All"][199])
print(y["Obs_J"][199]/y["Obs_All"][199])

# These will be ICs for Post Drug...
def declare_newparameters():
    postDrugIC_A = y["Obs_A"][199]
    postDrugIC_B = y["Obs_B"][199]
    postDrugIC_C = y["Obs_C"][199]
    postDrugIC_D = y["Obs_D"][199]
    postDrugIC_E = y["Obs_E"][199]
    postDrugIC_F = y["Obs_F"][199]
    postDrugIC_G = y["Obs_G"][199]
    postDrugIC_H = y["Obs_H"][199]
    postDrugIC_I = y["Obs_I"][199]
    postDrugIC_J = y["Obs_J"][199]
    postDrugIC_All = y["Obs_A"][199] + y["Obs_B"][199] + y["Obs_C"][199] + y["Obs_D"][199] + y["Obs_E"][199] + y["Obs_F"][199] + y["Obs_G"][199] + y["Obs_H"][199] + y["Obs_I"][199] + y["Obs_J"][199]
    
    alias_model_components()