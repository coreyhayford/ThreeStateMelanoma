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


def declare_monomers():
    Monomer("Cell")

def declare_parameters():
    Parameter("Cell_0", 1)
    Parameter("k_div", 0.03)
    Parameter("k_death", 0.005)

def declare_ICs():
    Initial(Cell, Cell_0)
    # Initial(Cell, Cell_0_1)

def declare_observables():
    Observable("Obs_Cell", Cell)

def declare_rules():
    Rule('div_Cell', Cell() >> Cell() + Cell(), k_div)
    Rule('death_Cell', Cell() >> None, k_death)

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
x = sim.run(tspan=t, verbose=False, n_sim=100, cleanup = False,)

y1 = np.array(x.observables)
initials = np.array(x.species)
tout = x.tout
y1_final = np.array([y for y in y1["Obs_Cell"].T[-1] if y != 0])
# y1_final = np.array([np.log2(y) for y in y1["Obs_Cell"].T[-1] if y != 0])
print(len(y1_final))
print(y1_final)

plt.figure("Final Pre-Drug Cell Counts")
plt.hist(y1_final)

plt.figure("Time Courses")
plt.plot(tout.T, y1['Obs_Cell'].T, '0.5', lw=4, alpha=0.25)  # individual trajectories
# plt.yscale('log', base = 2)
