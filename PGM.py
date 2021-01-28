from pysb import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sp
from pysb.simulator.bng import BngSimulator
import pandas as pd
import seaborn as sns
import math

# Import DIP rate distribution data
cFP_rates = pd.read_csv("cFP_rates_VUlines.csv")

# Initalize division and death rates (2 each - 2 states)
div = [0.032, 0.033]
dth = [0.0311, 0.0265]

# Function to run PySB monoclonal growth model (MGM)
## Remember stochasticity...not all simulations will look the same
def dist_compare(div, dth):
    ## Pull out only rates from DS8
    dat = cFP_rates[cFP_rates['Cell_Line'] == 'PC9-DS8']
    ## Initiate division and death rates
    kdiv = div
    kdth = dth
    ## cFP assays always start with 1 cell
    num_cells = 1
    ## Initiate model with components
    ### Model = cells divide and die according to rates associated with state
    Model()
    Monomer('Cell', ['type'], {'type': ["x"+str(i) for i in range(len(kdiv))]})
    [Initial(Cell(type="x"+str(i)), Parameter('Cell%d_Init' % i, num_cells)) for i in range(len(kdiv))]
    [Rule('divide_Cell%d' % i, Cell(type="x"+str(i)) >> Cell(type="x"+str(i)) + Cell(type="x"+str(i)),
          Parameter('kdiv_%d' % i, kdiv[i])) for i in range(len(kdiv))]
    [Rule('die_Cell%d' % i, Cell(type="x"+str(i)) >> None,
          Parameter('kdth_%d' % i, kdth[i])) for i in range(len(kdiv))]
    Observable("Cell_total", Cell())
    [Observable("Cell_t_%s" % i, Cell(type="x"+str(i))) for i in range(len(kdiv))]
    ## Collect simulated DIP rate distributions in list
    sim_dist = []
    ## Run first part of simulation (pre-drug: set division and death rates)
    t1 = np.linspace(0, 200, 201) #hrs
    sim1 = BngSimulator(model, tspan=t1, verbose=False)  # , seed = 1095205711)
    x1 = sim1.run(tspan=t1, param_values={'kdiv_%d' % i: 0.04 * np.log(2),
                                          'kdth_%d' % i: 0.005 * np.log(2)},
                 n_runs=len(dat), verbose=False)
    ## Capture response trajectories from first part of simulation
    trajs = np.array(np.array([tr[:]["Cell_total"] for tr in np.array(x1.observables)]).T)
    ## Initiate second part of simulation (drug-treated) with last time point from first
    ## part of simulation
    cell_tot = trajs[-1]
    ## Split cells into two states (1 has 75% of cells, 2 has 25% cells)
    cell_tot_num = len(cell_tot)
    cell_tot_s = np.random.multinomial(cell_tot_num, [3/4.,1/4.])
    ### Two helpful objects - not used later in code though
    prob_list = ([0]*int(cell_tot_s[0])) + ([1]*int(cell_tot_s[1]))
    cell_pop_n = [np.random.multinomial(int(round(cell)), [1,0]) for i,cell in enumerate(cell_tot)]
    ### Randomly distribute cells according to proportional split between states above
    cell_pop1 = [np.random.multinomial(int(round(cell)), [1,0]) for i,cell in enumerate(cell_tot)][:int(cell_tot_s[0])]
    cell_pop2 = [np.random.multinomial(int(round(cell)), [0,1]) for i,cell in enumerate(cell_tot)][int(cell_tot_s[0]):]
    ## Run second part of simulation and capture trajectories - state 1
    trajects = []
    for ind, cell in enumerate(cell_tot[:int(cell_tot_s[0])]):
        print("%d" % ind)
        cell_pop_1n = cell_pop1[ind]
        t2 = np.linspace(0, 225, 226)  # in drug
        x2 = sim1.run(tspan=t2,
                      initials={model.species[i]: pop for i, pop in enumerate(cell_pop_1n)},
                      verbose=False)
        ## Put DIP rates into list (if not zero)
        ### Keep DIP rate if not zero, keep trajectory if end of first simulation > 50 cells
        if cell != 0:
            slope, intercept, r_value, p_value, std_err = sp.stats.linregress(t2, np.log2(
                x2.observables["Cell_total"] / x2.observables["Cell_total"][0]))
            if math.isnan(slope) == False:
                sim_dist = sim_dist + [slope]
        if cell > 50:
            trajects.append(np.log2(x2.observables["Cell_total"] / x2.observables["Cell_total"][0]))
    ## Run second part of simulation and capture trajectories - state 2
    for ind, cell in enumerate(cell_tot[int(cell_tot_s[0]):]):
        print("%d" % ind)
        cell_pop_2n = cell_pop2[ind]
        t3 = np.linspace(0, 225, 226)  # in drug
        x3 = sim1.run(tspan=t3,
                      initials={model.species[i]: pop for i, pop in enumerate(cell_pop_2n)},
                      verbose=False)
        ## Put DIP rates into list (if not zero)
        ### Keep DIP rate if not zero, keep trajectory if end of first simulation > 50 cells
        if cell != 0:
            slope, intercept, r_value, p_value, std_err = sp.stats.linregress(t3, np.log2(
                x3.observables["Cell_total"] / x3.observables["Cell_total"][0]))
            if math.isnan(slope) == False:
                sim_dist = sim_dist + [slope]
        if cell > 50:
            trajects.append(np.log2(x3.observables["Cell_total"] / x3.observables["Cell_total"][0]))
    ## Capture p-value from K-S test comparison of experimental and simulated data
    D_stat, p_val = sp.ks_2samp(dat['DIP_Rate'], sim_dist)
    ## Return key elements
    return trajects, sim_dist, p_val

# Run function for DS8 - PGM) at specified example rates and capture key data
data_DS8 = dist_compare([0.032, 0.033],[0.0311, 0.0265])

## Output model trajectories
trajectories_DS8 = pd.DataFrame(data_DS8[0])
trajectories_DS8 = trajectories_DS8.transpose()
trajectories_DS8.to_csv('trajectories_DS8_G50.csv')

## Output model distributions
distributions_DS8 = pd.DataFrame({'DS8': data_DS8[1]})
distributions_DS8.to_csv('distributions_DS8_G50.csv')

## Output model p-value (won't plot if just one - manually inserted it)
# data_DS8[2].to_csv('pval_DS8.csv')