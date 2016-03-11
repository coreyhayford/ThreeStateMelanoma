from pysb import *
from pysb.util import *
from pysb.macros import *
from pysb.bng import *
from pysb.core import *
from pysb.integrate import odesolve
from pylab import *
import pylab as pl
from numpy import linspace
from sympy import sympify
from scipy import constants 
import matplotlib.pyplot as plt


 
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

#quit()


t = linspace(0,8000, 800)

y = odesolve(model,t,verbose=True)

plt.figure()
for obs in ["Obs_A", "Obs_B", "Obs_C"]:
    plt.plot(t, y[obs], label=re.match(r"Obs_(\w+)", obs).group(1), linewidth=3)
plt.legend(loc=0, prop={'size': 16})
plt.xlabel("Time", fontsize=22)
plt.ylabel("Cell Population", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.title("Three-State Model", fontsize=22)

plt.show()
quit() 
 
plt.figure()

for obs in ["Obs_A", "Obs_B", "Obs_C"]:
    plt.plot(t, y[obs], label=re.match(r"Obs_(\w+)", obs).group(1), linewidth=3)
plt.legend(loc=0, prop={'size': 16})
plt.xlabel("Time", fontsize=22)
plt.ylabel("Cell Population", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.title("Three-State Model", fontsize=22)

plt.show()

#pl.savefig("three_state_model", format= "png")