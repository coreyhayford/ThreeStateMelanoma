"""
Simulation of Barcoded Cell Dynamics in the absence and presence of drug treatment.
Implemented By: Corey Hayford
Support By: James Pino (BNG SSA Simulator Branch)

"""

from pysb import *
from pysb.simulator.bng_ssa import BngSimulator
from pysb.integrate import odesolve
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sp

n_clones = 17

def declare_monomers():
    Monomer('Cell', ['dip'], {'dip': [str(i) for i in range(1, n_clones)]})
    # Monomer("Cell", ['state'], {'state': ['dip_1', 'dip_2', 'dip_3', 'dip_4', 'dip_5', 'dip_6', 'dip_7', 'dip_8',
    #                                       'dip_9', 'dip_10', 'dip_11', 'dip_12', 'dip_13', 'dip_14', 'dip_15',
    #                                       'dip_16']})


def declare_parameters():
    # Parameter("Cell_dip_1_0", 1)
    # Parameter("Cell_dip_2_0", 0)
    # Parameter("Cell_dip_3_0", 0)
    # Parameter("Cell_dip_4_0", 0)
    # Parameter("Cell_dip_5_0", 0)
    # Parameter("Cell_dip_6_0", 0)
    # Parameter("Cell_dip_7_0", 0)
    # Parameter("Cell_dip_8_0", 0)
    # Parameter("Cell_dip_9_0", 0)
    # Parameter("Cell_dip_10_0", 0)
    # Parameter("Cell_dip_11_0", 0)
    # Parameter("Cell_dip_12_0", 0)
    # Parameter("Cell_dip_13_0", 0)
    # Parameter("Cell_dip_14_0", 0)
    # Parameter("Cell_dip_15_0", 0)
    # Parameter("Cell_dip_16_0", 0)
    Parameter("dip_0", 0)
    Parameter("dip_0_1", 1.0)

    Parameter("k_div", 0.03)
    Parameter("k_death", 0.005)
    # Parameter("k_dip_1_death", 0.005)
    # Parameter("k_dip_2_death", 0.005)
    # Parameter("k_dip_3_death", 0.005)
    # Parameter("k_dip_4_death", 0.005)
    # Parameter("k_dip_5_death", 0.005)
    # Parameter("k_dip_6_death", 0.005)
    # Parameter("k_dip_7_death", 0.005)
    # Parameter("k_dip_8_death", 0.005)
    # Parameter("k_dip_9_death", 0.005)
    # Parameter("k_dip_10_death", 0.005)
    # Parameter("k_dip_11_death", 0.005)
    # Parameter("k_dip_12_death", 0.005)
    # Parameter("k_dip_13_death", 0.005)
    # Parameter("k_dip_14_death", 0.005)
    # Parameter("k_dip_15_death", 0.005)
    # Parameter("k_dip_16_death", 0.005)

    # Parameter("k_dip_1_div", 0.03)
    # Parameter("k_dip_1_death", 0.005)
    # Parameter("k_dip_2_div", 0.03)
    # Parameter("k_dip_2_death", 0.005)
    # Parameter("k_dip_3_div", 0.03)
    # Parameter("k_dip_3_death", 0.005)
    # Parameter("k_dip_4_div", 0.03)
    # Parameter("k_dip_4_death", 0.005)
    # Parameter("k_dip_5_div", 0.03)
    # Parameter("k_dip_5_death", 0.005)
    # Parameter("k_dip_6_div", 0.03)
    # Parameter("k_dip_6_death", 0.005)
    # Parameter("k_dip_7_div", 0.03)
    # Parameter("k_dip_7_death", 0.005)
    # Parameter("k_dip_8_div", 0.03)
    # Parameter("k_dip_8_death", 0.005)
    # Parameter("k_dip_9_div", 0.03)
    # Parameter("k_dip_9_death", 0.005)
    # Parameter("k_dip_10_div", 0.03)
    # Parameter("k_dip_10_death", 0.005)
    # Parameter("k_dip_11_div", 0.03)
    # Parameter("k_dip_11_death", 0.005)
    # Parameter("k_dip_12_div", 0.03)
    # Parameter("k_dip_12_death", 0.005)
    # Parameter("k_dip_13_div", 0.03)
    # Parameter("k_dip_13_death", 0.005)
    # Parameter("k_dip_14_div", 0.03)
    # Parameter("k_dip_14_death", 0.005)
    # Parameter("k_dip_15_div", 0.03)
    # Parameter("k_dip_15_death", 0.005)
    # Parameter("k_dip_16_div", 0.03)
    # Parameter("k_dip_16_death", 0.005)


def declare_ICs():
    Initial(Cell(dip=str(1)), dip_0_1)
    for i in range(2, n_clones):
        Initial(Cell(dip=str(i)), dip_0)
    # Initial(Cell(dip=str(i)) for i in range(1, n_clones), 0)
#     Initial(Cell(dip=str(1)), dip_0)
    # Initial(Cell(state='dip_1'), Cell_dip_1_0)
    # Initial(Cell(state='dip_2'), Cell_dip_2_0)
    # Initial(Cell(state='dip_3'), Cell_dip_3_0)
    # Initial(Cell(state='dip_4'), Cell_dip_4_0)
    # Initial(Cell(state='dip_5'), Cell_dip_5_0)
    # Initial(Cell(state='dip_6'), Cell_dip_6_0)
    # Initial(Cell(state='dip_7'), Cell_dip_7_0)
    # Initial(Cell(state='dip_8'), Cell_dip_8_0)
    # Initial(Cell(state='dip_9'), Cell_dip_9_0)
    # Initial(Cell(state='dip_10'), Cell_dip_10_0)
    # Initial(Cell(state='dip_11'), Cell_dip_11_0)
    # Initial(Cell(state='dip_12'), Cell_dip_12_0)
    # Initial(Cell(state='dip_13'), Cell_dip_13_0)
    # Initial(Cell(state='dip_14'), Cell_dip_14_0)
    # Initial(Cell(state='dip_15'), Cell_dip_15_0)
    # Initial(Cell(state='dip_16'), Cell_dip_16_0)

    # Initial(Cell, Cell_0_1)

def declare_observables():
    Observable("Obs_Cell", Cell(dip=str(1)))

def declare_rules():
    Rule("Divide", Cell(str(1)) >> Cell() + Cell(), k_div)
    Rule("Death", Cell() >> None, k_death)
    # Rule('div_Cell', Cell(state = 'dip_1') >> Cell(state = 'dip_1') + Cell(state = 'dip_1'), k_div)
    # Rule('death_Cell', Cell(state = 'dip_1') >> None, k_dip_1_death)

Model()
declare_monomers()
declare_parameters()
declare_ICs()
declare_observables()
declare_rules()

print(model.monomers)
print(model.parameters)
print(model.observables)
print(model.rules)
for ic in model.initial_conditions:
    print ic


# We will integrate from t=0 to t=200
t = np.linspace(0, 200, 201)
# Simulate the model


print("Simulating Pre Drug...")
sim = BngSimulator(model, cleanup = False)
x = sim.run(tspan=t, verbose=False, n_sim=100, cleanup = False)

y1 = np.array(x.observables)
initials = np.array(x.species)
print(np.shape(initials))
print(initials[0,-1,:])

tout = x.tout
y1_final = np.array([np.log2(y) for y in y1["Obs_Cell"].T[-1] if y != 0])
print(len(y1_final))
print(y1_final)
# quit()
plt.figure("Final Pre-Drug Cell Counts")
plt.hist(y1_final)

plt.figure("Time Courses")
plt.plot(tout.T, np.log2(y1['Obs_Cell'].T), '0.5', lw=4, alpha=0.25)  # individual trajectories
# plt.yscale('log', base = 2)

quit()
plt.axvline(x=t[-1], color='red', ls='--', lw=4)


# Multinomial Choice of DIP (death now) rate from normal distribution around SC01-like DIP rate
# plt.figure("DIP Rate Histogram")
dips = np.linspace(-0.062, -0.002, 16)
print(dips)
normal = sp.norm.pdf(dips, -0.033, 0.01)
print(normal)
# plt.plot(dips, normal, lw = 2)

sum = 0
for i in range(1, len(dips)):
    sum += normal[i] * (dips[i]-dips[i-1])
print(sum)
normal_hist = normal * (dips[1]-dips[0])
print(normal_hist)
print(2**y1_final)
for y1_f in y1_final:
    # plt.figure()
    dip_state = np.random.multinomial(int(2**y1_f), normal_hist) # *100 for true normal
    print(dip_state)
    # plt.bar(dips-0.5*(dips[1]-dips[0]), 1.*dip_state/int(2**y1_f), dips[1]-dips[0], label = "n_cells = %d" %int(2**y1_f)) # *100 for true normal
    # plt.plot(dips, normal_hist, lw = 2, color = "r")
    # plt.legend(loc = 0)

    # dips_avg = dips-0.5*(dips[1]-dips[0])
    # dips_avg = []
    # for i in range(1, len(dips)):
    #     print(dips[i])
    #     # dips_avg.append(dips[i]+ (0.5*(dips[i+1] - dips[i])))
    #     dips_avg.append((dips[i] + dips[i - 1]) / 2.)
    # dips_avg_array = np.array(dips_avg)
    # print(dips_avg_array)
    print(dips)
    print(dip_state[0], dip_state[15])

    x1 = sim.run(tspan=t, verbose=False, n_sim=1, cleanup=False,
                    initials=initials[y1_f, -1, :],
                    param_values={"Cell_dip_1_0": dip_state[0],
                                  "Cell_dip_2_0": dip_state[1],
                                  "Cell_dip_3_0": dip_state[2],
                                  "Cell_dip_4_0": dip_state[3],
                                  "Cell_dip_5_0": dip_state[4],
                                  "Cell_dip_6_0": dip_state[5],
                                  "Cell_dip_7_0": dip_state[6],
                                  "Cell_dip_8_0": dip_state[7],
                                  "Cell_dip_9_0": dip_state[8],
                                  "Cell_dip_10_0": dip_state[9],
                                  "Cell_dip_11_0": dip_state[10],
                                  "Cell_dip_12_0": dip_state[11],
                                  "Cell_dip_13_0": dip_state[12],
                                  "Cell_dip_14_0": dip_state[13],
                                  "Cell_dip_15_0": dip_state[14],
                                  "Cell_dip_16_0": dip_state[15],
                                  "k_dip_1_death": -dips[0],
                                  "k_dip_2_death": -dips[1],
                                  "k_dip_3_death": -dips[2],
                                  "k_dip_4_death": -dips[3],
                                  "k_dip_5_death": -dips[4],
                                  "k_dip_6_death": -dips[5],
                                  "k_dip_7_death": -dips[6],
                                  "k_dip_8_death": -dips[7],
                                  "k_dip_9_death": -dips[8],
                                  "k_dip_10_death": -dips[9],
                                  "k_dip_11_death": -dips[10],
                                  "k_dip_12_death": -dips[11],
                                  "k_dip_13_death": -dips[12],
                                  "k_dip_14_death": -dips[13],
                                  "k_dip_15_death": -dips[14],
                                  "k_dip_16_death": -dips[15],
                                  })
    y2 = np.array(x1.observables)
    tout = x1.tout
    plt.plot(tout.T+200, np.log2(y2['Obs_Cell'].T), '0.5', lw=2, alpha=0.25)  # individual trajectories
# print(np.shape(initials))
# print(initials[:, -1, :])
plt.show()
quit()
###################

#     dip_state_index = [i for i, x in enumerate(dip_state) if x == 1]
#     dip_choice = dips_avg_array[dip_state_index]
#     print(dip_choice)
#     quit()
#     # post drug sims
#     x1 = sim.run(tspan=t, verbose=False, n_sim=1, cleanup=False,
#                 initials=initials[:, -1, :],
#                 param_values={"k_death": -dip_choice[0] + 0.05})
#     # dip states = 16 = bins
#     # monomers = # states
#     # n_clones = 16
#     # count up tot cells
#
# plt.show()
# quit()