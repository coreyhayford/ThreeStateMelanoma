from pysb import *
from pysb.util import *
from pysb.macros import *
from pysb.bng import *
from pysb.core import *
from pysb.integrate import odesolve
#import pylab as pl
from numpy import linspace
from sympy import sympify
from scipy import constants 

 
import three_state

Model()  
    
three_state.declare_monomers()
three_state.declare_parameters()
three_state.declare_initial_conditions()
three_state.declare_observables()
three_state.declare_functions()


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
    
generate_equations(model, verbose=True)

print len(model.rules)
print len(model.initial_conditions)
print len(model.reactions)
print len(model.species)
quit()