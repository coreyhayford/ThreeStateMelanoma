"""
Simulations of barcoded cell dynamics in the absense and presence of drug
Implemented By: Corey Hayford
Support By: James Pino (BNG SSA Simulator Branch)

"""
from pysb import *

from pysb.simulator.bng_ssa import BngSimulator
from pysb.integrate import odesolve
import numpy as np
import matplotlib.pyplot as plt

import barcoding_model

Model()
barcoding_model.declare_monomers()
barcoding_model.declare_parameters()
barcoding_model.declare_ICs()
barcoding_model.declare_observables()
barcoding_model.declare_rules()


# We will integrate from t=0 to t=200
t = np.linspace(0, 200, 201)
# Simulate the model
print("Simulating Pre Drug...")
sim = BngSimulator(model, cleanup = False)
num_sims = 2
for i in range(num_sims):
    x = sim.run(tspan=t, verbose=False, n_sim=1, cleanup = False,)

    y1 = np.array(x.observables)
    initials = np.array(x.species)
    print(initials)
    tout = x.tout
    plt.plot(tout.T, y1['Obs_Cell'].T, '0.5', lw=2, alpha=0.25, label=None)  # individual trajectories

    plt.axvline(x=t[-1], color='red', ls='--', lw=4)

    # Cell_0.value = y1["Obs_Cell"][-1]

    print(model.parameters)

    dips = np.random.normal(-0.033, 0.03, 100)
    count, bins, ignored = plt.hist(dips, 10, normed=True)
    count_normalized = (count / sum(count))

    bins_avg = []
    for i in range(1, len(bins + 1)):
        bins_avg.append((bins[i] + bins[i - 1]) / 2.)
    bins_avg_array = np.array(bins_avg)

    dip_state = np.random.multinomial(1, count_normalized)
    dip_state_index = [i for i, x in enumerate(dip_state) if x == 1]
    dip_choice = bins_avg_array[dip_state_index]

    x = sim.run(tspan=t, verbose=False, n_sim=1, cleanup=False,
                param_values={"k_death": -dip_choice[0] + 0.05, "Cell_0": y1["Obs_Cell"][-1]})

    y2 = np.array(x.observables)
    print(model.parameters)
    print(y2)
    plt.plot(tout.T + 200, y2['Obs_Cell'].T, '0.5', lw=2, alpha=0.25, label=None)  # individual trajectories

plt.show()
quit()
# print(initials[:,-1,0])
# print(type(y1[-1]))
# print(y1[-1][0])

# Multinomial Choice of DIP (death now) rate from normal distribution around SC01-like DIP rate
dips = np.random.normal(-0.033, 0.03, 100)
count, bins, ignored = plt.hist(dips, 10, normed=True)
count_normalized = (count / sum(count))

bins_avg = []
for i in range(1, len(bins + 1)):
    bins_avg.append((bins[i] + bins[i - 1]) / 2.)
bins_avg_array = np.array(bins_avg)

dip_state = np.random.multinomial(1, count_normalized)
dip_state_index = [i for i,x in enumerate(dip_state) if x == 1]
dip_choice = bins_avg_array[dip_state_index]

x = sim.run(tspan=t, verbose = False, n_sim=100, cleanup = False,
            initials=initials[:,-1,:],
            param_values={"k_death":-dip_choice[0] + 0.05})

y2 = np.array(x.observables)
tout = x.tout



# y3 = odesolve(model = model, tspan = t)
# print(initials)
# print(model.species)
# quit()
# y4 = odesolve(model = model, tspan = t, #y0=y3[-1]["Obs_Cell"],
#               param_values={"k_death":-dip_choice[0] + 0.05})
#
# # print(y3["Obs_Cell"])
# print(initials)
# print(initials[:,-1,:])
# print(y3[-1]["Obs_Cell"])
# print(np.array(y3))
# # print(y3["ObsCell"][-1])
# quit()
# print(y3)
# quit()
# print(y4["Obs_Cell"])
# print(y3)
# print(y3["Obs_Cell"])
# print(y3.shape)
# print(y3[0][0])
# quit()
# print(y4)
# print(y4[2])
# quit()

plt.figure()
plt.title("Barcoding Sims")
plt.plot(tout.T, y1['Obs_Cell'].T, '0.5', lw=2, alpha=0.25, label=None)  # individual trajectories
plt.plot(tout.T, y1['Obs_Cell'].T.mean(1), 'k-*', lw=3, label="Mean")
plt.plot(tout.T, y1['Obs_Cell'].T.min(1), 'b--', lw=3, label="Minimum")
plt.plot(tout.T, y1['Obs_Cell'].T.max(1), 'm--', lw=3, label="Maximum")
# plt.plot(t, y3['Obs_Cell'], 'g--', lw=3, label="ODE")

plt.axvline(x=t[-1], color='red', ls='--', lw=4)


plt.plot(tout.T+200, y2['Obs_Cell'].T, '0.5', lw=2, alpha=0.25, label=None)  # individual trajectories
plt.plot(tout.T+200, y2['Obs_Cell'].T.mean(1), 'k-*', lw=3, label="Mean")
plt.plot(tout.T+200, y2['Obs_Cell'].T.min(1), 'b--', lw=3, label="Minimum")
plt.plot(tout.T+200, y2['Obs_Cell'].T.max(1), 'm--', lw=3, label="Maximum")
# plt.plot(t+200, y4['Obs_Cell'], 'g--', lw=3, label="ODE")

plt.yscale('log', base=2)
plt.xlabel("Time")
plt.ylabel("Log2(Cell Count)")
# plt.legend(loc = 0)
plt.show()
