from pysb import *
from pysb.util import *
from pysb.macros import *
from pysb.bng import *
from pysb.integrate import odesolve
import pylab as pl
from numpy import linspace
from sympy import sympify

Model()

def declare_monomers():   
    Monomer("A")
    Monomer("B")
    Monomer("C")
    
def declare_parameters():
    
    ### Initial Conditions ###
    
    Parameter("A_0", 1500)
    Parameter("B_0", 1500)
    Parameter("C_0", 0)
    
    ### Growth Parameters ###
    Parameter("k_dip_A", -0.033)
    Parameter("k_dip_B",  0.000)
    Parameter("k_dip_C",  0.016)
     
#     Parameter("k_div_A", 0.033)
#     Parameter("k_div_B",  0)
#     Parameter("k_div_C",  0.016)
#     
#     Parameter("k_death_A", 0.066)
#     Parameter("k_death_B",  0)
#     Parameter("k_death_C",  0)
    
    Parameter("k_AB",     0.025)
    Parameter("k_BA",     0.00004)
    Parameter("k_CB",     0.025)
    Parameter("k_BC",     0.00004)
    
    alias_model_components()
    
def declare_initial_conditions():
    
    Initial(A, A_0)
    Initial(B, B_0)
    Initial(C, C_0)
    
def declare_observables():
    
    Observable("Obs_A", A)
    Observable("Obs_B", B)
    Observable("Obs_C", C)
    Observable("Obs_AB", A() + B())
    Observable("Obs_BC", B() + C())
    Observable("Obs_AC", A() + C())
    Observable("Obs_All", A() + B() + C())

def declare_functions():
    
#     Rule('Divide_A', A() >> A() + A(), k_div_A)
#     Rule('Death_A',  A() >> None,  k_death_A)
#     Rule('Divide_B', B() >> B() + B(), k_div_B)
#     Rule('Death_B',  B() >> None,  k_death_B)
#     Rule('Divide_C', C() >> C() + C(), k_div_C)
#     Rule('Death_C',  C() >> None,  k_death_C)
    
    Rule('DIP_A', A() >> A() + A(), k_dip_A)
    Rule('DIP_B', B() >> B() + B(), k_dip_B)
    Rule('DIP_C', C() >> C() + C(), k_dip_C)
    
    Rule("A_to_B", A() >> B(),       k_AB)
    Rule("B_to_A", B() >> A(),       k_BA)
    Rule("B_to_C", B() >> C(),       k_BC)
    Rule("C_to_B", C() >> B(),       k_CB)
       
 
