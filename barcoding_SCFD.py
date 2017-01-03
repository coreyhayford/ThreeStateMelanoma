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

    Parameter("A_0", 10)
    Parameter("B_0", 10)
    Parameter("C_0", 10)
    Parameter("D_0", 10)
    Parameter("E_0", 10)
    Parameter("F_0", 10)
    Parameter("G_0", 10)
    Parameter("H_0", 10)
    Parameter("I_0", 10)
    Parameter("J_0", 10)
   
#     for j in xrange(1,10):
#         Parameter(j + "_0", 10)
#     
#     # Rate Parameters

    Parameter("k_dip_A", 0.02 * np.random.randn() + 0)
    Parameter("k_dip_B", 0.02 * np.random.randn() + 0)
    Parameter("k_dip_C", 0.02 * np.random.randn() + 0)
    Parameter("k_dip_D", 0.02 * np.random.randn() + 0)
    Parameter("k_dip_E", 0.02 * np.random.randn() + 0)
    Parameter("k_dip_F", 0.02 * np.random.randn() + 0)
    Parameter("k_dip_G", 0.02 * np.random.randn() + 0)
    Parameter("k_dip_H", 0.02 * np.random.randn() + 0)
    Parameter("k_dip_I", 0.02 * np.random.randn() + 0)
    Parameter("k_dip_J", 0.02 * np.random.randn() + 0)
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
    
    Rule('DIP_A', A() >> A() + A(), k_dip_A)
    Rule('DIP_B', B() >> B() + B(), k_dip_B)
    Rule('DIP_C', C() >> C() + C(), k_dip_C)
    Rule('DIP_D', D() >> D() + D(), k_dip_D)
    Rule('DIP_E', E() >> E() + E(), k_dip_E)
    Rule('DIP_F', F() >> F() + F(), k_dip_F)
    Rule('DIP_G', G() >> G() + G(), k_dip_G)
    Rule('DIP_H', H() >> H() + H(), k_dip_H)
    Rule('DIP_I', I() >> I() + I(), k_dip_I)
    Rule('DIP_J', J() >> J() + J(), k_dip_J)

#     
#     for n in xrange(1,10):
#         Rule("DIP_" + n, n >> n + n, "k_dip_" + n)
    
