from pysb import *
from pysb.util import *
from pysb.macros import *
from pysb.simulator.bng_ssa import BngSimulator
from pysb.integrate import odesolve
import pylab as pl
from numpy import linspace
from sympy import sympify
import scipy.stats
import numpy as np
import matplotlib.pyplot as plt

# def plot(t, t_start=0, label=False):
#     plt.plot(tout.T+t_start, y_sim, lw=2, color = 'grey', label = ('DIP = %g' %(k_divide_Cell.value-k_death_Cell.value)) if label else None)

# def plot_mean_min_max(tout, trajectories, param_values=None, title=None,
#                       savename=None, model=None, tspan=None):
#
#     x = np.array([tr[:] for tr in trajectories]).T
#     tout = tout.T
#     plt.figure()
#     plt.title(title)
#     plt.subplot(111)
#     plt.plot(tout, x, '0.5', lw=2, alpha=0.25)  # individual trajectories
#     plt.plot(tout, x.mean(1), 'k-*', lw=3, label="Mean")
#     plt.plot(tout, x.min(1), 'b--', lw=3, label="Minimum")
#     plt.plot(tout, x.max(1), 'r--', lw=3, label="Maximum")
#     if model is not None:
#         y = odesolve(model, tspan, param_values)
#         plt.plot(tout, y[title], 'g--', lw=3, label="ODE")
#     plt.tight_layout()
#     if savename is not None:
#         plt.savefig("%s.png" % savename)
#     if show:
#         plt.show()
#     plt.close()

def declare_monomers():
    Monomer("Cell")

def declare_parameters():
    Parameter("Cell_0", 1)
    # Parameter("Cell_0_1", y_sim["Obs_Cell"][-1])

    Parameter("k_div_noDrug", 0.05)
    Parameter("k_death_noDrug", 0.025)

def declare_ICs():
    Initial(Cell, Cell_0)
    # Initial(Cell, Cell_0_1)

def declare_observables():
    Observable("Obs_Cell", Cell)

def declare_rules():
    Rule('div_Cell', Cell() >> Cell() + Cell(), k_div_noDrug)
    Rule('death_Cell', Cell() >> None, k_death_noDrug)

# def declare_drug_ICs():
#     Parameter("postDrugIC_Cell", y_sim[-1][0])

def declare_drug_parameters():
    Parameter("k_divide_Cell", 0.05)
    Parameter("k_death_Cell", -dip_choice[0] + 0.05)

def declare_drug_rules():
    Rule('divide_Cell', Cell() >> Cell() + Cell(), k_divide_Cell)
    Rule('die_Cell', Cell() >> None, k_death_Cell)

# alias_model_components()

# plt.figure("Pre-drug Stochastic Simulations")
####
#
# def plot(t, t_start=0, label=False):
#     plt.plot(t + t_start, y["Obs_Cell"], 'c-', lw=2,
#              label=('Cell DIP = %g' % (k_divide_Cell.value - k_death_Cell.value)) if label else None)
#
#     plt.xlabel("Time")
#     plt.ylabel("Cell Count")
#     plt.legend(loc=0)  # , fontsize = 6)
#     plt.ylim([0,1000])


Model()
declare_monomers()
declare_parameters()
declare_ICs()
declare_observables()
declare_rules()

# y = odesolve(model, t, verbose=True)
# y = run_ssa(model, t_end=t[-1], n_steps=len(t) - 1, verbose=True, n_sims = 1000)
# plot(t, label=False)

# We will integrate from t=0 to t=200
t = np.linspace(0, 200, 201)
# Simulate the model
print("Simulating Pre Drug...")
sim = BngSimulator(model, cleanup = False)
x = sim.run(tspan=t, verbose=False, n_sim=3, cleanup = False,)
# x_sim1 = sim1.run(tspan=t, t_end=t[-1], n_steps=len(t) - 1, verbose=False, n_sim=3)
y1 = np.array(x.observables)
print(y1)
# tout1 = x_sim1.tout
initials = np.array(x.species)
# print(initials)
print(initials[:,-1, :])
quit()
# print(initials)
# print(y_sim1)
# print(initials[:,-1,0])

# print(y_sim[2])
# print(type(y_sim))
# print(y_sim[-1])
# print(type(y_sim[-1]))
# print(y_sim[-1][0])

# plot_mean_min_max(tout, y_sim,
#                   model=model, tspan=t)

# plt.plot(tout1.T, y_sim1.T)
#
# plt.show()
#
# quit()
# plot(t, label=False)

plt.axvline(x=t[-1], color='red', ls='--', lw=4)


# Multinomial Choice of DIP (death now) rate from normal distribution around SC01-like DIP rate
dips = np.random.normal(-0.033, 0.01, 100)
# print(dips)
count, bins, ignored = plt.hist(dips, 10, normed=True)
# print(count)
# print(bins)

# print(sum(count))

count_normalized = (count / sum(count))
# print(count_normalized)
# print(sum(count_normalized))

bins_avg = []
for i in range(1, len(bins + 1)):
    # print(i)
    # print(a[i])
    # print(a[i-1])
    # print((bins[i]+bins[i-1])/2.)
    bins_avg.append((bins[i] + bins[i - 1]) / 2.)
bins_avg_array = np.array(bins_avg)
# print(bins_avg)
# print(bins_avg_array)

dip_state = np.random.multinomial(1, count_normalized)
# print(dip_state)
dip_state_index = [i for i,x in enumerate(dip_state) if x == 1]
# print(dip_state_index)
dip_choice = bins_avg_array[dip_state_index]
# print(dip_choice[0])

# count_list = np.ndarray.tolist(count_normalized)
# print(count_list)
# dip = scipy.stats.multinomial(n=10, p=count_list)
# dip_dist = dip.pmf(bins_avg)
# print(dip_dist)

# These will be ICs for Post Drug...

print(model.initial_conditions)

# Cell_0_ICs = []
# for j,y in enumerate(y1):
#     Cell_0.value = y1[j][-1]
#     # print(Cell_0.value)
#     Cell_0_ICs.append(Cell_0.value)
#     # print(Cell_0_ICs)
# # print(Cell_0.value)
# Cell_0_ICs = np.array(Cell_0_ICs)
# print(Cell_0_ICs)
# declare_drug_ICs()
declare_drug_parameters()
declare_drug_rules()

print(model.monomers)
print(model.parameters)
print(model.initial_conditions)
print(model.observables)
print(model.rules)

# Simulate the model
print("Simulating Post Drug...")
# sim2 = BngSimulator(model)
# x_sim = sim.run(tspan=t, initials=Cell_0_ICs, t_end=t[-1], n_steps=len(t) - 1, verbose=False, n_sim=3)

# print(self._model.param_values)
# print(self._model.param_values.shape)
print(model.parameters)
# quit()
x = sim.run(tspan=t, verbose = False, n_sim=3, cleanup = False, initials=initials[:,-1,:])
initials = np.array(x.species)
tout = x.tout
y2 = np.array(x_sim2.observables)
# print(y_sim2)
# plot(t, t_start=t[-1], label=True)
plt.plot(tout.T, y1["Obs_Cell"].T, 'r')
plt.plot(tout.T+200, y2["Obs_Cell"].T, 'b')
# plt.plot(tout.T[-1] + tout.T, y2.T, 'b')
# y = run_ssa(model, t_end=t[-1], n_steps=len(t) - 1, verbose=True, n_sims = 1000)
# plot(t, t_start=t[-1], label=False)

# print(A_0.value, B_0.value, C_0.value, D_0.value, E_0.value, F_0.value, G_0.value, H_0.value, I_0.value, J_0.value)

# plt.tight_layout()
plt.show()
plt.show()
## DONE: established initial conditions as a function of gillespie sims
## DONE: Pick DIP rates from multinomial distribution
## DONE: Plot before and after drug simulated results
## DONE: Apply new plot function to stitch together results
## DONE: Put in line for drug introduction
## DONE: Stitch together subplots and make into figure - made in a loop
## DONE: Incorporate James' new bng simulator with only 1 cell and run 1000s of sims - not stitched yet


## GOAL:
## GOAL: Plot the clonal distributions like Hata
## GOAL: Make it for 1 million barcoded cells rather than 10
## GOAL: Play with different DIP rates/distributions to model diversification phenotype
## GOAL: Only pick time points similar to experiment??
## GOAL: Include state transitions in the model

## GOAL: Add more barcodes and plot like Hata
