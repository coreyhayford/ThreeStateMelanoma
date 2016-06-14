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
    
for exp in model.expressions:
    print exp
    
generate_equations(model, verbose=True)
from pysb.generator.bng import BngGenerator
print BngGenerator(model).get_content()

for i in range(len(model.odes)):
    print str(i)+":", model.odes[i]





print len(model.rules)
print len(model.initial_conditions)
print len(model.reactions)
print model.species

#quit()


t = linspace(0, 100, 100)

y = odesolve(model,t,verbose=True)
#y = run_ssa(model,t_end = t[-1], n_steps = len(t)-1, verbose=True)
 
plt.figure()

#for obs in ["Obs_A", "Obs_B", "Obs_C"]:
#    plt.plot(t, y[obs], label=re.match(r"Obs_(\w+)", obs).group(1), linewidth=3)
#for obs in ["Obs_AB", "Obs_BC", "Obs_AC"]:
#    plt.plot(t, y[obs], '--', label=re.match(r"Obs_(\w+)", obs).group(1), linewidth=3)
plt.plot(t, y["Obs_AC"], 'k--', lw=3, label="Total")
plt.legend(loc=0, prop={'size': 16})
plt.xlabel("Time", fontsize=22)
plt.ylabel("Cell Population", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.title("Three-State Model", fontsize=22)

# plt.show()
plt.savefig("three_state_model_mix.pdf", format= "pdf")


plt.figure()

# fig, ax = plt.subplots()
# ax.set_yscale('log', basey=2)
# ax.plot(range(100))

#for obs in ["Obs_A", "Obs_B", "Obs_C"]:
#    plt.plot(t, y[obs], label=re.match(r"Obs_(\w+)", obs).group(1), linewidth=3)
#for obs in ["Obs_AB", "Obs_BC", "Obs_AC"]:
#    plt.plot(t, y[obs], '--', label=re.match(r"Obs_(\w+)", obs).group(1), linewidth=3)
plt.plot(t, np.log2(y["Obs_All"]/y["Obs_All"][0]), 'b--', lw=3, label="1:1:1")
#plt.yscale('log', basey=2)

plt.xlabel("Time", fontsize=22)
plt.ylabel("Cell Population log2", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.title("Three-State Model", fontsize=22)

#########

model.parameters["A_0"].value = 600
model.parameters["B_0"].value = 1800
model.parameters["C_0"].value = 600

y = odesolve(model, t, verbose=True)
plt.plot(t, np.log2(y["Obs_All"]/y["Obs_All"][0]), 'r--', lw=3, label="1:3:1")
plt.legend(loc=0, prop={'size': 16})

#plt.savefig("three_state_model_mix_log2.pdf", format= "pdf")

plt.show()