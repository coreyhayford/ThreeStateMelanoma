from pysb import *
from pysb.util import *
from pysb.macros import *
from pysb.bng import *
from pysb.integrate import odesolve
import pylab as pl
from numpy import linspace
from sympy import sympify
import numpy as np


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
#     for i in xrange(1,10):
#         
#         print(i)
#         Monomer(i)
#         
#         print(model.monomers())
#         
# quit()
def declare_parameters():
#     
#     # Initial Conditions

    Parameter("A_0", 1000)
    Parameter("B_0", 1000)
    Parameter("C_0", 1000)
    Parameter("D_0", 1000)
    Parameter("E_0", 1000)
    Parameter("F_0", 1000)
    Parameter("G_0", 1000)
    Parameter("H_0", 1000)
    Parameter("I_0", 1000)
    Parameter("J_0", 1000)
   
#     for j in xrange(1,10):
#         Parameter(j + "_0", 10)
#     
#     # Rate Parameters

    Parameter("k_div_A", 0.02 * np.random.randn() + 0.06)
    Parameter("k_div_B", 0.02 * np.random.randn() + 0.06)
    Parameter("k_div_C", 0.02 * np.random.randn() + 0.06)
    Parameter("k_div_D", 0.02 * np.random.randn() + 0.06)
    Parameter("k_div_E", 0.02 * np.random.randn() + 0.06)
    Parameter("k_div_F", 0.02 * np.random.randn() + 0.06)
    Parameter("k_div_G", 0.02 * np.random.randn() + 0.06)
    Parameter("k_div_H", 0.02 * np.random.randn() + 0.06)
    Parameter("k_div_I", 0.02 * np.random.randn() + 0.06)
    Parameter("k_div_J", 0.02 * np.random.randn() + 0.06)
    
    Parameter("k_death_A", 0.02 * np.random.randn() + 0.09)
    Parameter("k_death_B", 0.02 * np.random.randn() + 0.09)
    Parameter("k_death_C", 0.02 * np.random.randn() + 0.09)
    Parameter("k_death_D", 0.02 * np.random.randn() + 0.09)
    Parameter("k_death_E", 0.02 * np.random.randn() + 0.09)
    Parameter("k_death_F", 0.02 * np.random.randn() + 0.09)
    Parameter("k_death_G", 0.02 * np.random.randn() + 0.09)
    Parameter("k_death_H", 0.02 * np.random.randn() + 0.09)
    Parameter("k_death_I", 0.02 * np.random.randn() + 0.09)
    Parameter("k_death_J", 0.02 * np.random.randn() + 0.09)
    
#     for k in xrange(1,10):
#         Parameter("k_dip_" + k, 0.02 * np.random.randn(1,10) + 0)
#         
    alias_model_components()

def declare_initial_conditions():
    
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
#     
#     for l in xrange(1,10):
#         Initial(l, l + "_0")

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
#     
#     for m in xrange(1,10):
#         Observable("Obs_" + m, m)
#         
#     #Observable("Obs_All", )
    
def declare_functions():
    
    Rule('Divide_A', A() >> A() + A(), k_div_A)
    Rule('Divide_B', B() >> B() + B(), k_div_B)
    Rule('Divide_C', C() >> C() + C(), k_div_C)
    Rule('Divide_D', D() >> D() + D(), k_div_D)
    Rule('Divide_E', E() >> E() + E(), k_div_E)
    Rule('Divide_F', F() >> F() + F(), k_div_F)
    Rule('Divide_G', G() >> G() + G(), k_div_G)
    Rule('Divide_H', H() >> H() + H(), k_div_H)
    Rule('Divide_I', I() >> I() + I(), k_div_I)
    Rule('Divide_J', J() >> J() + J(), k_div_J)
    
    Rule('Death_A', A() >> None, k_death_A)
    Rule('Death_B', B() >> None, k_death_B)
    Rule('Death_C', C() >> None, k_death_C)
    Rule('Death_D', D() >> None, k_death_D)
    Rule('Death_E', E() >> None, k_death_E)
    Rule('Death_F', F() >> None, k_death_F)
    Rule('Death_G', G() >> None, k_death_G)
    Rule('Death_H', H() >> None, k_death_H)
    Rule('Death_I', I() >> None, k_death_I)
    Rule('Death_J', J() >> None, k_death_J)

#     
#     for n in xrange(1,10):
#         Rule("DIP_" + n, n >> n + n, "k_dip_" + n)
    
