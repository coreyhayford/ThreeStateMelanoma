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
import scipy as sp
import pandas as pd

#import seaborn as sns; sns.set(color_codes=True)
# x = np.linspace(0, 15, 31)
# data = np.sin(x) + np.random.rand(10, 31) + np.random.randn(10, 1)
# print(data)
# ax = sns.tsplot(data=data)
# plt.show()
# quit()

import three_state_dip


Model()  
    
three_state_dip.declare_monomers()
three_state_dip.declare_parameters()
three_state_dip.declare_initial_conditions()
three_state_dip.declare_observables()
three_state_dip.declare_functions()
  
generate_equations(model, verbose=True)
from pysb.generator.bng import BngGenerator
print BngGenerator(model).get_content()

for i in range(len(model.odes)):
    print str(i)+":", model.odes[i]

t = linspace(0, 200, 200)

import matplotlib.cm as cm
import matplotlib.colors as mplcolors
logCC = np.arange(-5.0, 2.0, 0.1)
colors = [cm.spectral(i) for i in np.linspace(0.2, 0.9, len(logCC)+1)]

plt.figure()
ax = plt.subplot(111)

y = odesolve(model,t,verbose=True)
# df = pd.DataFrame(columns = range(1,201))
#y = odesolve(model,t_end = t[-1], n_steps = len(t)-1, verbose=True)
    #print type(y)
    #ax = sns.tsplot(data=np.log2(y['Obs_All']/y['Obs_All'][0]), time = y['time'],condition = "100% SC01", err_style = "ci_band", ci=95, color="m")
#plt.plot(y['time'], np.log2(y['Obs_All']/y['Obs_All'][0]), 'r-', lw=3, label="100% SC07" if sims == 0 else "") #colors[int(i)])

plt.plot(t, np.log2(y['Obs_All']/y['Obs_All'][0]), 'r--', lw=3, label="100% Subclone 01") #colors[int(i)])





model.parameters["A_0"].value = 0
model.parameters["B_0"].value = 3000
model.parameters["C_0"].value = 0


y = odesolve(model,t,verbose=True)
plt.plot(t, np.log2(y['Obs_All']/y['Obs_All'][0]), 'g--', lw=3, label="100% Subclone 07") #colors[int(i)])


model.parameters["A_0"].value = 0  
model.parameters["B_0"].value = 0    
model.parameters["C_0"].value = 3000

y = odesolve(model,t,verbose=True)
plt.plot(t, np.log2(y['Obs_All']/y['Obs_All'][0]), 'b--', lw=3, label="100% Subclone 10") #colors[int(i)]
#  dashes = (10,5/10/15)

model.parameters["A_0"].value = 1000   
model.parameters["B_0"].value = 1000    
model.parameters["C_0"].value = 1000

y = odesolve(model,t,verbose=True)
plt.plot(t, np.log2(y['Obs_All']/y['Obs_All'][0]), 'k-', lw=5, label="Parental") #colors[int(i)]

    
plt.xlabel("Time (hours)", fontsize=22)
plt.ylabel("Population Doublings (nl2)", fontsize=18)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
#plt.ylim(-1.25, 1.25)
plt.xlim(0, 200)

plt.title("Model Simulations", fontsize=20)
# Shrink current axis's height by 10% on the bottom
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])

# Put a legend below current axis
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
          fancybox=True, shadow=True, ncol=5)
#plt.legend(loc=0, prop={'size': 16})
plt.show()
quit()

#y = odesolve(model,t,verbose=True)

# y1 = run_ssa(model,t_end = t[-1], n_steps = len(t)-1, verbose=True)
# y2 = run_ssa(model,t_end = t[-1], n_steps = len(t)-1, verbose=True)
# y3 = run_ssa(model,t_end = t[-1], n_steps = len(t)-1, verbose=True)


# plt.figure()
# 
# #for obs in ["Obs_A", "Obs_B", "Obs_C"]:
# #    plt.plot(t, y[obs], label=re.match(r"Obs_(\w+)", obs).group(1), linewidth=3)
# #for obs in ["Obs_AB", "Obs_BC", "Obs_AC"]:
# #    plt.plot(t, y[obs], '--', label=re.match(r"Obs_(\w+)", obs).group(1), linewidth=3)
# plt.plot(t, y1["Obs_All"], 'r-', lw=3, label="Total")
# plt.plot(t, y2["Obs_All"], 'g-', lw=3, label="Total")
# plt.plot(t, y3["Obs_All"], 'b-', lw=3, label="Total")
# plt.legend(loc=0, prop={'size': 16})
# plt.xlabel("Time", fontsize=22)
# plt.ylabel("Cell Population", fontsize=22)
# plt.xticks(fontsize=18)
# plt.yticks(fontsize=18)
# plt.title("Three-State Model: Stochastic", fontsize=22)
# 
# # plt.show()
# plt.savefig("three_state_model_nonnormal_stochastic_131.pdf", format= "pdf")


# fig, ax = plt.subplots()
# ax.set_yscale('log', basey=2)
# ax.plot(range(100))

#for obs in ["Obs_A", "Obs_B", "Obs_C"]:
#    plt.plot(t, y[obs], label=re.match(r"Obs_(\w+)", obs).group(1), linewidth=3)
#for obs in ["Obs_AB", "Obs_BC", "Obs_AC"]:
#    plt.plot(t, y[obs], '--', label=re.match(r"Obs_(\w+)", obs).group(1), linewidth=3)
#ssa_list = []
#print ssa_list
# SSA simulations (run each in triplicate)
#ssa_list.append(3*[None])
#print ssa_list

        
# fig = plt.figure('% Drugged Cells (SSA)')
# plt.plot(y['time'][1:], y['Cell_drug'][1:]/y['Cell_total'][1:]*100., lw=3, color=colors[i])

#plt.legend(loc=0, prop={'size': 16})

#plt.savefig("three_state_model_mix_log2.pdf", format= "pdf")

# plt.show()