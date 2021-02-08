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

t = linspace(0, 200, 200)

import matplotlib.cm as cm
import matplotlib.colors as mplcolors
logCC = np.arange(-5.0, 2.0, 0.1)
colors = [cm.spectral(i) for i in np.linspace(0.2, 0.9, len(logCC)+1)]



plt.figure()
# samples_time = []
# 
# 
# # df = pd.DataFrame(columns = range(1,201))
# for sims in range(5):
#     y = run_ssa(model,t_end = t[-1], n_steps = len(t)-1, verbose=True)
#     #print type(y)
#     #ax = sns.tsplot(data=np.log2(y['Obs_All']/y['Obs_All'][0]), time = y['time'],condition = "100% SC01", err_style = "ci_band", ci=95, color="m")
#     plt.plot(y['time'], np.log2(y['Obs_All']/y['Obs_All'][0]), 'r-', lw=3, label="100% SC07" if sims == 0 else "") #colors[int(i)])
#     samples_time.append(y['Obs_All']/y['Obs_All'][0])
#     #df.loc[1] = pd.Series(y['time'])
# #     df.loc[sims+1] = pd.Series(y['Obs_All']/y['Obs_All'][0])
#     #sample_mean = mean(y['Obs_All']/y['Obs_All'][0])
#     #quit()
# print (samples_time)
# 
# mean_time = []
# sem_time = []
# 
# for time_point_samples in samples_time:
#     mean_time.append(np.mean(time_point_samples))
#     sem_time.append(np.sem(time_point_samples))
# 
# # print df
# # df_t = pd.DataFrame.transpose(df)
# # print df_t
# # sns.tsplot(data = df_t, color = "indianred")
# # data = {'time': y['time']}
# # quit()
# #print(np.mean(np.log2(y['Obs_All']/y['Obs_All'][0])))
# df = pd.DataFrame()
# #ax = sns.tsplot(data=np.log2(y['Obs_All']/y['Obs_All'][0]), ci=[95], color="m")
# #sample_mean = mean()
# #sample_95ci = 1.96*


for sims in range(5):
    y = run_ssa(model,t_end = t[-1], n_steps = len(t)-1, verbose=True)
    plt.plot(y['time'], np.log2(y['Obs_All']/y['Obs_All'][0]), 'r-', lw=3, 
              label="100% A" if sims == 0 else "") #colors[int(i)])


model.parameters["A_0"].value = 2250
model.parameters["B_0"].value = 450
model.parameters["C_0"].value = 300


for sims in range(5):
    y = run_ssa(model,t_end = t[-1], n_steps = len(t)-1, verbose=True)
    plt.plot(y['time'], np.log2(y['Obs_All']/y['Obs_All'][0]), 'g--', lw=3, dashes = (10,5),
             label="75% A, 15% B, 10%C" if sims == 0 else "") #colors[int(i)]


model.parameters["A_0"].value = 1500  
model.parameters["B_0"].value = 900    
model.parameters["C_0"].value = 600

for sims in range(5):
    y = run_ssa(model,t_end = t[-1], n_steps = len(t)-1, verbose=True)
    plt.plot(y['time'], np.log2(y['Obs_All']/y['Obs_All'][0]), 'b--', lw=3, dashes = (10,10),
              label="50% A, 30% B, 20% C" if sims == 0 else "") #colors[int(i)]


model.parameters["A_0"].value = 750   
model.parameters["B_0"].value = 1500    
model.parameters["C_0"].value = 750
    
for sims in range(5):
    y = run_ssa(model,t_end = t[-1], n_steps = len(t)-1, verbose=True)
    plt.plot(y['time'], np.log2(y['Obs_All']/y['Obs_All'][0]), 'm--', lw=3, dashes = (10,15),
              label="25% A, 50% B, 25% C" if sims == 0 else "") #colors[int(i)]

    
plt.xlabel("Time (hours)", fontsize=22)
plt.ylabel("Population Doublings (nl2)", fontsize=18)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.ylim(-1.5, 1.0)
#plt.title("Model Simulations", fontsize=20)
plt.legend(loc=2)
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