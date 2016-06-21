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
import matplotlib.cm as cm
import matplotlib.colors as mplcolors



 
import three_state

Model()  
    
three_state.declare_monomers()
three_state.declare_parameters()
three_state.declare_initial_conditions()
three_state.declare_observables()
three_state.declare_functions()
  
generate_equations(model, verbose=True)
from pysb.generator.bng import BngGenerator
print BngGenerator(model).get_content()

for i in range(len(model.odes)):
    print str(i)+":", model.odes[i]

t = linspace(0, 100, 100)

logCC = np.arange(-5.0, 2.0, 0.1)
colors = [cm.spectral(i) for i in np.linspace(0.2, 0.9, len(logCC)+1)]
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
ssa_list = []
#print ssa_list
# SSA simulations (run each in triplicate)
ssa_list.append(3*[None])
#print ssa_list
 
for sims in range(20):
    y = run_ssa(model,t_end = t[-1], n_steps = len(t)-1, verbose=True)
    #y = run_ssa(model, t_span, verbose=False)
#     try:
#         # throw out the first 10% of points (just to be safe)
#         slope, intercept, r_value, p_value, std_err = ss.linregress(y['time'][100:], np.log(y['Cell_total'][100:])) # Note: natural log
#         except RuntimeWarning:
#             pass
#         else:
#             slopes_ssa[-1][j] = slope
                
plt.figure('Total Cells (Stochastic)')
plt.plot(y['time'], np.log2(y['Obs_All']/y['Obs_All'][0]), lw=3, color=colors[int(i)])
plt.xlabel("Time", fontsize=22)
plt.ylabel("Cell Population (nl2)", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.title("Three-State Model: Stochastic", fontsize=22)
plt.show()
quit()
        
# fig = plt.figure('% Drugged Cells (SSA)')
# plt.plot(y['time'][1:], y['Cell_drug'][1:]/y['Cell_total'][1:]*100., lw=3, color=colors[i])

plt.plot(t, np.log2(y1["Obs_All"]/y1["Obs_All"][0]), 'b-', lw=3, label="1:1:1")
plt.plot(t, np.log2(y2["Obs_All"]/y2["Obs_All"][0]), 'b-', lw=3, label="1:1:1")
plt.plot(t, np.log2(y3["Obs_All"]/y3["Obs_All"][0]), 'b-', lw=3, label="1:1:1")
#plt.yscale('log', basey=2)

plt.xlabel("Time", fontsize=22)
plt.ylabel("Cell Population (nl2)", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.title("Three-State Model: Stochastic", fontsize=22)

#########

model.parameters["A_0"].value = 600
model.parameters["B_0"].value = 1800
model.parameters["C_0"].value = 600

z1 = run_ssa(model,t_end = t[-1], n_steps = len(t)-1, verbose=True)
z2 = run_ssa(model,t_end = t[-1], n_steps = len(t)-1, verbose=True)
z3 = run_ssa(model,t_end = t[-1], n_steps = len(t)-1, verbose=True)
plt.plot(t, np.log2(z1["Obs_All"]/z1["Obs_All"][0]), 'r--', lw=3, label="1:3:1")
plt.plot(t, np.log2(z2["Obs_All"]/z2["Obs_All"][0]), 'r--', lw=3, label="1:3:1")
plt.plot(t, np.log2(z3["Obs_All"]/z3["Obs_All"][0]), 'r--', lw=3, label="1:3:1")
#plt.legend(loc=0, prop={'size': 16})

#plt.savefig("three_state_model_mix_log2.pdf", format= "pdf")

plt.show()