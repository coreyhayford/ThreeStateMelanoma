from pysb import *
from pysb.util import *
from pysb.macros import *
from pysb.bng import run_ssa
from pysb.core import *
from pysb.integrate import odesolve
from pysb.simulator.bng import BngSimulator
from pylab import *
import pylab as pl
import numpy as np
import scipy as sp
from numpy import linspace
from sympy import sympify
from scipy import constants 
import matplotlib.pyplot as plt
from ggplot import *

from scipy.stats import sem

def plot_mean_min_max(name, title = None, color = "red"):
    trajectories = y1
    x = np.array([tr[:][name] for tr in trajectories]).T
    tout = x1.tout.T
    # plt.figure()
    plt.title(title)
    plt.subplot(111)
    mean = np.log2(x.mean(1)/x.mean(1)[0])
    ci_95 = 1.96 * sem(np.log2(x/x[0]).T)
    ci_min = mean - ci_95
    ci_max = mean + ci_95
    # print(mean)
    # print(ci_95)
    # print(mean.shape)
    # print(np.log2(x/x[0]))
    # print(ci_95)
    # print(ci_95.shape)
    # print(tout.shape)
    # print(ci_min.shape)
    # print(ci_max.shape)
    # min = np.log2(x.min(1)/x.min(1)[0])
    # max = np.log2(x.max(1)/x.max(1)[0])
    plt.fill_between(t, ci_min, ci_max, facecolor = color, alpha = 0.25)
    plt.plot(tout, np.log2(x/x[0]), '0.5', color = color, lw=1, alpha=0.25)  # individual trajectories
    # plt.plot(tout, np.log2(x.mean(1)/x.mean(1)[0]), 'k-*', lw=3)#, label="Mean")
    # plt.plot(tout, np.log2(x.min(1)/x.min(1)[0]), 'b--', lw=3)#, label="Minimum")
    # plt.plot(tout, np.log2(x.max(1)/x.max(1)[0]), 'r--', lw=3)#, label="Maximum")
    # if model is not None:
    #     y = odesolve(model, tspan = t)
    #     plt.plot(tout, y[title], 'g--', lw=3, label="ODE")
    # plt.tight_layout()
    # if savename is not None:
    #     plt.savefig("%s.png" % savename)
    # if show:
    #     plt.show()
    # plt.close()

# def mean_confidence_interval(data, confidence=0.95):
#     a = 1.0*np.array(data)
#     n = len(a)
#     m, se = np.mean(a), sp.stats.sem(a)
#     h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
#     return m, m-h, m+h

import three_state
# import three_state_dip

# def plot_mean_min_max(tout, trajectories, param_values=None, title=None,
#                       savename=None, model=None, tspan=None):


Model()  
    
three_state.declare_monomers()
three_state.declare_parameters()
three_state.declare_initial_conditions()
three_state.declare_observables()
three_state.declare_functions()

#
# for m in model.monomers:
#     print m
#
# for p in model.parameters:
#     print p
#
# for ic in model.initial_conditions:
#     print ic
#
# for obs in model.observables:
#     print obs
#
# for rules in model.rules:
#     print rules
#
# for exp in model.expressions:
#     print exp
#
# generate_equations(model, verbose=True)
# from pysb.generator.bng import BngGenerator
# print BngGenerator(model).get_content()
#
# for i in range(len(model.odes)):
#     print str(i)+":", model.odes[i]
#
#
#
#
# print len(model.rules)
# print len(model.initial_conditions)
# print len(model.reactions)
# print model.species

#quit()

plt.figure()
t = linspace(0, 200, 201)
#
sims = 5
sim = BngSimulator(model, tspan=t, verbose=False)
x1 = sim.run(tspan=t, verbose=False, n_runs=sims, cleanup=False, method='ssa')
# tout = x1.tout
y1 = np.array(x1.observables)
print(y1['Obs_All'])
print(y1['Obs_All'].T)
print((np.log2(y1['Obs_All'])).T)
print((np.log2(y1['Obs_All']).T/np.log2(y1['Obs_All']).T[0]))
print(x1.all[-1]["Obs_All"])




plot_mean_min_max('Obs_All', color = 'red')
y = odesolve(model = model,tspan = t, verbose = True)
plt.plot(t, np.log2(y["Obs_All"]/y["Obs_All"][0]), 'r-', lw=4, label = "100% SC01")

# [plt.plot(t, np.log2(x1.all[i]['Obs_All']/x1.all[i]['Obs_All'][0]), lw = 1, color = 'r') for i in range(sims)]#label = obs.name, color = colors[i])
# plt.plot(tout.T, (y1['Obs_All']).T, color = 'r')
# plt.plot(tout, y1['Obs_All'], color = 'b')

# y = odesolve(model = model,tspan = t, verbose = True)
# plt.plot(t, np.log2(y["Obs_All"]/y["Obs_All"][0]), 'r-', lw=3)


model.parameters["A_0"].value = 2550
model.parameters["B_0"].value = 300
model.parameters["C_0"].value = 150

sim = BngSimulator(model, tspan=t, verbose=False)
x1 = sim.run(tspan=t, verbose=False, n_runs=sims, cleanup=False, method='ssa')
tout = x1.tout
y1 = np.array(x1.observables)
# plt.plot(tout.T, y1['Obs_All'].T, color = 'g')

plot_mean_min_max('Obs_All', color = "green")
y = odesolve(model, t, verbose=True)
plt.plot(t, np.log2(y["Obs_All"]/y["Obs_All"][0]), 'g-', lw=3, dashes = (8,2), label = "85% SC01 : 10% SC07 : 5% SC10")
# plt.plot(t, y["Obs_All"], 'g--', lw=3)#, label="2:1:1")

model.parameters["A_0"].value = 2100
model.parameters["B_0"].value = 600
model.parameters["C_0"].value = 300

sim = BngSimulator(model, tspan=t, verbose=False)
x1 = sim.run(tspan=t, verbose=False, n_runs=sims, cleanup=False, method='ssa')
tout = x1.tout
y1 = np.array(x1.observables)
# plt.plot(tout.T, y1['Obs_All'].T, color = 'b')

plot_mean_min_max('Obs_All', color = "blue")
y = odesolve(model, t, verbose=True)
plt.plot(t, np.log2(y["Obs_All"]/y["Obs_All"][0]), 'b-', lw=3, dashes = (10,4), label = "70% SC01 : 20% SC07 : 10% SC10")
# plt.plot(t, y["Obs_All"], 'b--', lw=3, dashes=(10, 1))#, label="3:1:1")

model.parameters["A_0"].value = 750
model.parameters["B_0"].value = 1500
model.parameters["C_0"].value = 750

sim = BngSimulator(model, tspan=t, verbose=False)
x1 = sim.run(tspan=t, verbose=False, n_runs=sims, cleanup=False, method='ssa')
tout = x1.tout
y1 = np.array(x1.observables)
# plt.plot(tout.T, y1['Obs_All'].T, color = 'm')

plot_mean_min_max('Obs_All', color = "magenta")
y = odesolve(model, t, verbose=True)
plt.plot(t, np.log2(y["Obs_All"]/y["Obs_All"][0]), 'm-', lw=3, dashes = (12,6), label = "25% SC01 : 50% SC07 : 25% SC10")
# plt.plot(t, y["Obs_All"], 'm--', lw=3, dashes=(10, 1.5))#, label="3:1:1")


#plt.yscale('log', basey=2)

plt.xlabel("Time", fontsize=22)
plt.ylabel("Population Doublings (nl2)", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.ylim(-1.25,1.25)
plt.xlim(-5,205)


plt.title("Three-State Model", fontsize=22)

# plt.legend(loc=9, bbox_to_anchor=(0.5, -0.1), ncol=4)


plt.legend(loc=0, prop={'size': 16})
plt.show()
quit()
# y2 = run_ssa(model,t_end = t[-1], n_steps = len(t)-1, verbose=True)

# plt.figure()
# 
# #for obs in ["Obs_A", "Obs_B", "Obs_C"]:
# #    plt.plot(t, y[obs], label=re.match(r"Obs_(\w+)", obs).group(1), linewidth=3)
# #for obs in ["Obs_AB", "Obs_BC", "Obs_AC"]:
# #    plt.plot(t, y[obs], '--', label=re.match(r"Obs_(\w+)", obs).group(1), linewidth=3)
# plt.plot(t, y["Obs_AC"], 'k--', lw=3, label="Total")
# plt.legend(loc=0, prop={'size': 16})
# plt.xlabel("Time", fontsize=22)
# plt.ylabel("Cell Population", fontsize=22)
# plt.xticks(fontsize=18)
# plt.yticks(fontsize=18)
# 
# plt.title("Three-State Model", fontsize=22)
# 
# # plt.show()
# plt.savefig("three_state_model_mix.pdf", format= "pdf")




# fig, ax = plt.subplots()
# ax.set_yscale('log', basey=2)
# ax.plot(range(100))

#for obs in ["Obs_A", "Obs_B", "Obs_C"]:
#    plt.plot(t, y[obs], label=re.match(r"Obs_(\w+)", obs).group(1), linewidth=3)
#for obs in ["Obs_AB", "Obs_BC", "Obs_AC"]:
#    plt.plot(t, y[obs], '--', label=re.match(r"Obs_(\w+)", obs).group(1), linewidth=3)
plt.plot(t, np.log2(y["Obs_All"]/y["Obs_All"][0]), 'r-', lw=3)#, label="1:1:1")
# plt.plot(t, np.log2(y2["Obs_All"]/y2["Obs_All"][0]), 'r:', lw=3, label="1:1:1")


model.parameters["A_0"].value = 2550
model.parameters["B_0"].value = 300
model.parameters["C_0"].value = 150

y = odesolve(model, t, verbose=True)
plt.plot(t, np.log2(y["Obs_All"]/y["Obs_All"][0]), 'g--', lw=3)#, label="2:1:1")

model.parameters["A_0"].value = 2100
model.parameters["B_0"].value = 600
model.parameters["C_0"].value = 300

y = odesolve(model, t, verbose=True)
plt.plot(t, np.log2(y["Obs_All"]/y["Obs_All"][0]), 'b--', lw=3, dashes=(10, 1))#, label="3:1:1")

model.parameters["A_0"].value = 750
model.parameters["B_0"].value = 1500
model.parameters["C_0"].value = 750

y = odesolve(model, t, verbose=True)
plt.plot(t, np.log2(y["Obs_All"]/y["Obs_All"][0]), 'm--', lw=3, dashes=(10, 1.5))#, label="3:1:1")


#plt.yscale('log', basey=2)

plt.xlabel("Time", fontsize=22)
plt.ylabel("Population Doublings (nl2)", fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.ylim(-1.25,1.25)


plt.title("Three-State Model", fontsize=22)
plt.legend(loc=0, prop={'size': 16})
plt.show()
quit()
#########

model.parameters["A_0"].value = 600
model.parameters["B_0"].value = 1800
model.parameters["C_0"].value = 600

y = odesolve(model, t, verbose=True)
plt.plot(t, np.log2(y["Obs_All"]/y["Obs_All"][0]), 'r--', lw=3, label="1:3:1")

#plt.savefig("three_state_model_mix_log2.pdf", format= "pdf")

plt.show()