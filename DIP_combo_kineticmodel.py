from pysb import *
from pysb.integrate import odesolve
from pysb.simulator import ScipyOdeSimulator
from pysb.macros import equilibrate, bind
from sympy import sympify
import numpy as np
import matplotlib.pyplot as plt

Model()

Monomer("Cell", ['state'], {'state':['0','1','2','3']})
Monomer("Drug1")
Monomer("Drug2")

# print(model.monomers)

# Initial Params
Parameter("Cell_S0_0", 1000)
Parameter("Cell_S1_0", 0)
Parameter("Cell_S2_0", 0)
Parameter("Cell_S3_0", 0)

# DIP Rates
Parameter("k_dip_Cell_S0", 0.1)
Parameter("k_dip_Cell_S1", 0.1)
Parameter("k_dip_Cell_S2", 0.1)
Parameter("k_dip_Cell_S3", 0.1)

# Forward and Reverse Rate Constants
Parameter("kf_Cell_S0S1", 0.1)
Parameter("kf_Cell_S0S2", 0.1)
Parameter("kf_Cell_S1S3", 0.1)
Parameter("kf_Cell_S2S3", 0.1)

Parameter("kr_Cell_S1S0", 0.1)
Parameter("kr_Cell_S2S0", 0.1)
Parameter("kr_Cell_S3S1", 0.1)
Parameter("kr_Cell_S3S2", 0.1)

# Drug Concentrations
Parameter("Drug1_conc", 0.1)
Parameter("Drug2_conc", 0.1)

Parameter("alpha", 10)

# State Transition Params
Parameter("statetrans_C0C1", kf_Cell_S0S1.value * Drug1_conc.value)
Parameter("statetrans_C0C2", kf_Cell_S0S2.value * Drug2_conc.value)
Parameter("statetrans_C1C3", kf_Cell_S1S3.value * Drug2_conc.value * alpha.value)
Parameter("statetrans_C2C3", kf_Cell_S2S3.value * Drug1_conc.value * alpha.value)

print(model.parameters)

# Setting ICs
Initial(Cell(state="0"), Cell_S0_0)
Initial(Cell(state="1"), Cell_S1_0)
Initial(Cell(state="2"), Cell_S2_0)
Initial(Cell(state="3"), Cell_S3_0)

print(model.initial_conditions)

# Setting Observables
Observable("OBS_Cell_S0", Cell(state="0"))
Observable("OBS_Cell_S1", Cell(state="1"))
Observable("OBS_Cell_S2", Cell(state="2"))
Observable("OBS_Cell_S3", Cell(state="3"))
# Observable("Drug1", Drug1())
# Observable("Drug2", Drug2())

print(model.observables)

# Expressions if necessary
# Expression("statetrans_C0C1", sympify("kf_Cell_S0S1 * Drug1_conc"))
# Expression("statetrans_C0C2", sympify("kf_Cell_S0S2 * Drug2_conc"))
# Expression("statetrans_C1C3", sympify("kf_Cell_S1S3 * Drug2_conc * alpha"))
# Expression("statetrans_C2C3", sympify("kf_Cell_S2S3 * Drug1_conc * alpha"))

# DIP Rate Rules (div-death, can break apart into div and death rules)
# Here, cells inherit drugged state from parent - not sure if that's true
Rule("DIP_Cell_S0",Cell(state="0") >> Cell(state="0") + Cell(state="0"), k_dip_Cell_S0)
Rule("DIP_Cell_S1",Cell(state="1") >> Cell(state="1") + Cell(state="1"), k_dip_Cell_S1)
Rule("DIP_Cell_S2",Cell(state="2") >> Cell(state="2") + Cell(state="2"), k_dip_Cell_S2)
Rule("DIP_Cell_S3",Cell(state="3") >> Cell(state="3") + Cell(state="3"), k_dip_Cell_S3)

# State Transitions
Rule("trans_C0C1", Cell(state="0") + Drug1() >> Cell(state="1") + Drug1(), statetrans_C0C1)
Rule("trans_C0C2", Cell(state="0") + Drug2() >> Cell(state="2") + Drug2(), statetrans_C0C2)
Rule("trans_C1C3", Cell(state="1") + Drug2() >> Cell(state="3") + Drug2(), statetrans_C1C3)
Rule("trans_C2C3", Cell(state="2") + Drug1() >> Cell(state="3") + Drug1(), statetrans_C2C3)

Rule("trans_C1C0", Cell(state="1") >> Cell(state="0"), kr_Cell_S1S0)
Rule("trans_C2C0", Cell(state="2") >> Cell(state="1"), kr_Cell_S2S0)
Rule("trans_C3C1", Cell(state="3") >> Cell(state="1"), kr_Cell_S3S1)
Rule("trans_C3C2", Cell(state="3") >> Cell(state="2"), kr_Cell_S3S2)


print(model.rules)

plt.figure()
t = np.linspace(0, 200, 201)
y = ScipyOdeSimulator(model).run(tspan=t).all
# y = odesolve(model = model, tspan = t, verbose=True)
plt.plot(t, y["OBS_Cell_S0"], 'k', lw = 4)
plt.plot(t, y["OBS_Cell_S1"], 'r', lw = 2)
plt.plot(t, y["OBS_Cell_S2"], 'g', lw = 2)
plt.plot(t, y["OBS_Cell_S3"], 'b', lw = 2)

plt.show()