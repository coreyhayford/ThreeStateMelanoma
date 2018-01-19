from pysb import *
from pysb.integrate import odesolve
from pysb.macros import equilibrate, bind
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sp
from pysb.bng import run_ssa
# from pysb.simulator.bng_ssa import BngSimulator
from pysb.bng import generate_equations
import pandas as pd

alpha = 10

Model()

# Three Component Model
Monomer("Cell", ['a', 'b'])
Monomer("Drug_A", ['a'])
Monomer("Drug_B", ['b'])

print(model.monomers)

# Initial Params
Parameter("Cell_0", 1000)
Parameter("CellDrug_A0", 0)
Parameter("CellDrug_B0", 0)
Parameter("CellDrug_A0B0", 0)

# DIP Rates
Parameter("k_dip_C", 0.1)
Parameter("k_dip_CDa", 0.1)
Parameter("k_dip_CDb", 0.1)
Parameter("k_dip_CDaDb", 0.1)

# Forward and Reverse Rate Constants
Parameter("kf_trans_C_CDa", 0.1)
Parameter("kf_trans_C_CDb", 0.1)
Parameter("kf_trans_CDa_CDaDb", alpha * 0.1)
Parameter("kf_trans_CDb_CDaDb", alpha * 0.1)

Parameter("kr_trans_C_CDa", 0.1)
Parameter("kr_trans_C_CDb", 0.1)
Parameter("kr_trans_CDa_CDaDb", 0.1)
Parameter("kr_trans_CDb_CDaDb", 0.1)

# Can I put alpha as a parameter? Did not work when trying to multiply in rule...
# Parameter("alpha", 10)

print(model.parameters)

# Setting ICs
Initial(Cell(a=None, b=None), Cell_0)
Initial(Cell(a=1, b=None) % Drug_A(a=1), CellDrug_A0)
Initial(Cell(a=None, b=1) % Drug_B(b=1), CellDrug_B0)
Initial(Cell(a=1, b=1) % Drug_A(a=1) % Drug_B(b=1), CellDrug_A0B0)

print(model.initial_conditions)

# Setting Observables
Observable("OBS_Cell", Cell(a=None, b=None))
Observable("OBS_Cell_DrugA", Cell(a=1, b=None) % Drug_A(a=1))
Observable("OBS_Cell_DrugB", Cell(a=None, b=1) % Drug_B(b=1))
Observable("OBS_Cell_DrugAB", Cell(a=1, b=1) % Drug_A(a=1) % Drug_B(b=1))

print(model.observables)

# DIP Rate Rules (div-death, can break apart into div and death rules)
# Here, cells inherit drugged state from parent - not sure if that's true
Rule("DIP_Cell",Cell(a=None,b=None) >> Cell(a=None,b=None) + Cell(a=None,b=None), k_dip_C)
Rule("DIP_CellA",Cell(a=1,b=None) >> Cell(a=1,b=None) + Cell(a=1,b=None), k_dip_CDa)
Rule("DIP_CellB",Cell(a=None,b=1) >> Cell(a=None,b=1) + Cell(a=None,b=1), k_dip_CDb)
Rule("DIP_CellAB",Cell(a=1,b=1) >> Cell(a=1,b=1) + Cell(a=1,b=1), k_dip_CDaDb)

# Binding and Unbinding of Cell to Drugs A and B Rules
Rule("bind_Cell_DrugA",Cell(a=None,b=None) + Drug_A(a=None) >> Cell(a=1,b=None) % Drug_A(a=1), kf_trans_C_CDa)
Rule("bind_Cell_DrugB",Cell(a=None,b=None) + Drug_B(b=None) >> Cell(a=None,b=1) % Drug_B(b=1), kf_trans_C_CDb)
Rule("bind_CellA_DrugAB",Cell(a=1,b=None) % Drug_A(a=1) + Drug_B(b=None) >> Cell(a=1,b=1) % Drug_A(a=1) % Drug_B(b=1), kf_trans_CDa_CDaDb)
Rule("bind_CellB_DrugAB",Cell(a=None,b=1) % Drug_B(b=1) + Drug_A(a=None) >> Cell(a=1,b=1) % Drug_A(a=1) % Drug_B(b=1), kf_trans_CDb_CDaDb)

Rule("unbind_Cell_DrugA",Cell(a=1,b=None) % Drug_A(a=1) >> Cell(a=None, b=None) + Drug_A(a=None), kr_trans_C_CDa)
Rule("unbind_Cell_DrugB",Cell(a=None,b=1) % Drug_B(b=1) >> Cell(a=None, b=None) + Drug_B(b=None), kr_trans_C_CDb)
Rule("unbind_CellA_DrugAB",Cell(a=1,b=1) % Drug_A(a=1) % Drug_B(b=1) >> Cell(a=1,b=None) % Drug_A(a=1) + Drug_B(b=None), kr_trans_CDa_CDaDb)
Rule("unbind_CellB_DrugAB",Cell(a=1,b=1) % Drug_A(a=1) % Drug_B(b=1) >> Cell(a=None,b=1) % Drug_B(b=1) + Drug_A(a=None), kr_trans_CDb_CDaDb)
#
# bind(Cell(a=None,b=None),'a',Drug_A(a=1), 'a', [kf_trans_C_CDa, kr_trans_C_CDa])
# bind(Cell(a=None,b=None),'b',Drug_B(b=1), 'b', [kf_trans_C_CDb, kr_trans_C_CDb])
# bind(Cell(a=1,b=None),'b',Drug_B(b=1), 'b', [kf_trans_CDa_CDaDb, kr_trans_CDa_CDaDb])
# bind(Cell(a=None,b=1),'a',Drug_A(a=1), 'a', [kf_trans_CDb_CDaDb, kf_trans_CDb_CDaDb])

print(model.rules)
# quit()
plt.figure()
t = np.linspace(0, 200, 201)
y = odesolve(model = model, tspan = t, verbose=True)
plt.plot(t, "OBS_Cell", 'k', lw = 4)
plt.plot(t, "OBS_Cell_DrugA", 'r', lw = 2)
plt.plot(t, "OBS_Cell_DrugB", 'g', lw = 2)
plt.plot(t, "OBS_Cell_DrugAB", 'b', lw = 2)

plt.show()
