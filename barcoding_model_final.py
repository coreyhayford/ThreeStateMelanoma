from pysb import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sp
from pysb.bng import run_ssa
from pysb.simulator.bng_ssa import BngSimulator
from pysb.bng import generate_equations

n_cell_types = 16

Model()
[Monomer("Cell", ['dip'], {'dip': ["%d" %i for i in range(n_cell_types)]})]
print(model.monomers)

# Monomer("Cell", ['dip'], {"dip": ["0","1"]})


Parameter('cellInit_0', 1)
[Parameter('cellInit_%d' %i)for i in range(1,n_cell_types)]
# Parameter('Cell_1init')
print(model.parameters)

Initial(Cell(dip = "0"), cellInit_0)
Initial(Cell(dip = "1"), cellInit_1)
Initial(Cell(dip = "2"), cellInit_2)
Initial(Cell(dip = "3"), cellInit_3)
Initial(Cell(dip = "4"), cellInit_4)
Initial(Cell(dip = "5"), cellInit_5)
Initial(Cell(dip = "6"), cellInit_6)
Initial(Cell(dip = "7"), cellInit_7)
Initial(Cell(dip = "8"), cellInit_8)
Initial(Cell(dip = "9"), cellInit_9)
Initial(Cell(dip = "10"), cellInit_10)
Initial(Cell(dip = "11"), cellInit_11)
Initial(Cell(dip = "12"), cellInit_12)
Initial(Cell(dip = "13"), cellInit_13)
Initial(Cell(dip = "14"), cellInit_14)
Initial(Cell(dip = "15"), cellInit_15)

# [Initial(Cell(dip = "%d" %i), "cellInit_%d" %i) for i in range(1,n_cell_types)]
print(model.initial_conditions)

# Initial(Cell(dip = "0"), Cell_0init)
# Initial(Cell(dip = "1"), Cell_1init)

[Observable("Obs_Cell%d" %i, Cell(dip = "%d" %i)) for i in range(n_cell_types)]
print(model.observables)

# Observable("Obs_Cell0", Cell(dip = "0"))
# Observable("Obs_Cell1", Cell(dip = "1"))

low_dips = -0.062
high_dips = -0.002
dips = np.linspace(low_dips, high_dips, n_cell_types)
print(dips)
[Parameter("k_div_%d" % i, 0.03) for i in range(n_cell_types)]
k_death = -dips+0.03
[Parameter("k_death_%d" % i, k_death[i]) for i in range(n_cell_types)]

print(model.parameters)

# [Parameter("k_div_%d" % i, 0.03) for i in range(n_cell_types)]
# k_death = [0.04, 0.05]
# [Parameter("k_death_%d" % i, k_death[i]) for i in range(n_cell_types)]

Rule("Cell0_Div", Cell(dip = "0") >> Cell(dip = "0") + Cell(dip = "0"), k_div_0)
Rule("Cell0_Death", Cell(dip = "0") >> None, k_death_0)

Rule("Cell1_Div", Cell(dip = "1") >> Cell(dip = "1") + Cell(dip = "1"), k_div_1)
Rule("Cell1_Death", Cell(dip = "1") >> None, k_death_1)

Rule("Cell2_Div", Cell(dip = "2") >> Cell(dip = "2") + Cell(dip = "2"), k_div_2)
Rule("Cell2_Death", Cell(dip = "2") >> None, k_death_2)

Rule("Cell3_Div", Cell(dip = "3") >> Cell(dip = "3") + Cell(dip = "3"), k_div_3)
Rule("Cell3_Death", Cell(dip = "3") >> None, k_death_3)

Rule("Cell4_Div", Cell(dip = "4") >> Cell(dip = "4") + Cell(dip = "4"), k_div_4)
Rule("Cell4_Death", Cell(dip = "4") >> None, k_death_4)

Rule("Cell5_Div", Cell(dip = "5") >> Cell(dip = "5") + Cell(dip = "5"), k_div_5)
Rule("Cell5_Death", Cell(dip = "5") >> None, k_death_5)

Rule("Cell6_Div", Cell(dip = "6") >> Cell(dip = "6") + Cell(dip = "6"), k_div_6)
Rule("Cell6_Death", Cell(dip = "6") >> None, k_death_6)

Rule("Cell7_Div", Cell(dip = "7") >> Cell(dip = "7") + Cell(dip = "7"), k_div_7)
Rule("Cell7_Death", Cell(dip = "7") >> None, k_death_7)

Rule("Cell8_Div", Cell(dip = "8") >> Cell(dip = "8") + Cell(dip = "8"), k_div_8)
Rule("Cell8_Death", Cell(dip = "8") >> None, k_death_8)

Rule("Cell9_Div", Cell(dip = "9") >> Cell(dip = "9") + Cell(dip = "9"), k_div_9)
Rule("Cell9_Death", Cell(dip = "9") >> None, k_death_9)

Rule("Cell10_Div", Cell(dip = "10") >> Cell(dip = "10") + Cell(dip = "10"), k_div_10)
Rule("Cell10_Death", Cell(dip = "10") >> None, k_death_10)

Rule("Cell11_Div", Cell(dip = "11") >> Cell(dip = "11") + Cell(dip = "11"), k_div_11)
Rule("Cell11_Death", Cell(dip = "11") >> None, k_death_11)

Rule("Cell12_Div", Cell(dip = "12") >> Cell(dip = "12") + Cell(dip = "12"), k_div_12)
Rule("Cell12_Death", Cell(dip = "12") >> None, k_death_12)

Rule("Cell13_Div", Cell(dip = "13") >> Cell(dip = "13") + Cell(dip = "13"), k_div_13)
Rule("Cell13_Death", Cell(dip = "13") >> None, k_death_13)

Rule("Cell14_Div", Cell(dip = "14") >> Cell(dip = "14") + Cell(dip = "14"), k_div_14)
Rule("Cell14_Death", Cell(dip = "14") >> None, k_death_14)

Rule("Cell15_Div", Cell(dip = "15") >> Cell(dip = "15") + Cell(dip = "15"), k_div_15)
Rule("Cell15_Death", Cell(dip = "15") >> None, k_death_15)
# generate_equations(model, verbose=True)
# for sp in model.species:
#     print sp

print(model.rules)
quit()
t1 = np.linspace(0,200,201)
sim = BngSimulator(model, tspan=t1, verbose=False)#, seed = 1095205711)
x1 = sim.run(seed = 1095205711, param_values={"k_death_0": 0.005}) # returns np.array with species and obs
cell_tot = x1.all["Obs_Cell0"][-1]
cell_pop = [] #multinom here
cell_pop.append(sp.binom.rvs(n = cell_tot, p = 0.5))
cell_pop.append(cell_tot-cell_pop[0])
print("Sim1 finished.")
t2 = np.linspace(0,100,101)

x2 = sim.run(tspan = t2, param_values={"k_death_0": k_death_0.value},
             initials={model.species[i]:cell_pop[i] for i in range(len(cell_pop))})

# colors = ['r','b']
for i,obs in enumerate(model.observables):
    # plt.plot(t1, x1.all[obs.name], lw = 2, label = obs.name, color = colors[i])
    plt.plot(t1[-1]+t2, x2.all[obs.name], lw = 2, label = obs.name)
plt.annotate(s = "n_cells = %d" % cell_tot, xy = (0.1,0.1), xycoords = "axes fraction", fontsize = 24)
# plt.text(-0.04, 2800, r'$\mu=%g$' % (np.mean(distr_all)))
plt.legend(loc = 0)
plt.show()

