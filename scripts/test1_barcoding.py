from pysb import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sp
from pysb.bng import run_ssa
from pysb.simulator.bng_ssa import BngSimulator
from pysb.bng import generate_equations

Model()

Monomer("Cell", ['dip'], {"dip": ["0","1"]})

Parameter('Cell_0init', 1)
Parameter('Cell_1init')
Initial(Cell(dip = "0"), Cell_0init)
Initial(Cell(dip = "1"), Cell_1init)

Observable("Obs_Cell0", Cell(dip = "0"))
Observable("Obs_Cell1", Cell(dip = "1"))

n_cell_types = 2
[Parameter("k_div_%d" % i, 0.03) for i in range(n_cell_types)]
k_death = [0.04, 0.05]
[Parameter("k_death_%d" % i, k_death[i]) for i in range(n_cell_types)]

Rule("Cell0_Div", Cell(dip = "0") >> Cell(dip = "0") + Cell(dip = "0"), k_div_0)
Rule("Cell0_Death", Cell(dip = "0") >> None, k_death_0)

Rule("Cell1_Div", Cell(dip = "1") >> Cell(dip = "1") + Cell(dip = "1"), k_div_1)
Rule("Cell1_Death", Cell(dip = "1") >> None, k_death_1)

# generate_equations(model, verbose=True)
# for sp in model.species:
#     print sp

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

