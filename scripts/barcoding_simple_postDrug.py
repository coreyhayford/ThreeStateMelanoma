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

import barcoding_simple_preDrug as preDrugModel
print("New File...")
print(preDrugModel.y1_final)

cellCounts = preDrugModel.y1_final

# Multinomial Choice of DIP (death now) rate from normal distribution around SC01-like DIP rate
plt.figure("DIP Rate Histogram")
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
print(cellCounts)
for y1_f in cellCounts:
    plt.figure()
    dip_state = np.random.multinomial(int(y1_f), normal_hist) # *100 for true normal
    print(dip_state)
    plt.bar(dips-0.5*(dips[1]-dips[0]), 1.*dip_state/int(y1_f), dips[1]-dips[0], label = "n_cells = %d" %int(y1_f)) # *100 for true normal
    plt.plot(dips, normal_hist, lw = 2, color = "r")
    plt.legend(loc = 0)

n_clones = len(dips)

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

    # Parameter("k_div", 0.03)
    # Parameter("k_death", 0.005)
    for index,count in enumerate(cellCounts):
        print(index, int(count))
        print("Initial_Count_%d" %index, int(cellCounts[index]))
        Parameter("Initial_Count_%d" %index, int(cellCounts[index]))


Model()
declare_parameters()
print(model.parameters)
quit()
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



# def declare_ICs():
#     Initial(Cell(dip=str(1)), dip_0_1)
#     for i in range(2, n_clones):
#         Initial(Cell(dip=str(i)), dip_0)
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

quit()

plt.show()