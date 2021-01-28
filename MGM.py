from pysb import *
import numpy as np
import scipy.stats as sp
from pysb.simulator.bng import BngSimulator
import pandas as pd
import math

# Import DIP rate distribution data
cFP_rates = pd.read_csv("cFP_rates_VUlines.csv")

# Function to run PySB monoclonal growth model (MGM)
## Remember stochasticity...not all simulations will look the same
def dist_compare(samp, div, dth):
    ## Pull out only rates from that subline
    dat = cFP_rates[cFP_rates['Cell_Line'] == samp]
    ## Initiate division and death rates
    kdiv = div
    kdth = dth
    ## cFP assays always start with 1 cell
    num_cells = 1
    ## Initiate model with components
    ### Model = 1 cell divides and dies with specified rate
    Model()
    Monomer('Cell')
    Parameter('Cell_init', num_cells)
    Initial(Cell, Cell_init)
    Observable('Obs_Cell', Cell())
    Parameter('k_div', kdiv)
    Parameter('k_dth', kdth)
    Rule('Division', Cell() >> Cell() + Cell(), k_div)
    Rule('Death', Cell() >> None, k_dth)
    ## Collect simulated DIP rate distributions in list
    sim_dist = []
    ## Run first part of simulation (pre-drug: set division and death rates)
    t1 = np.linspace(0, 200, 201) #hrs
    sim1 = BngSimulator(model, tspan=t1, verbose=False)
    x1 = sim1.run(tspan=t1, param_values={'k_div': 0.04 * np.log(2),
                                          'k_dth': 0.005 * np.log(2)},
                  n_runs=len(dat), verbose=False)
    ## Capture response trajectories from first part of simulation
    trajs = np.array(np.array([tr[:]["Obs_Cell"] for tr in np.array(x1.observables)]).T)
    ## Initiate second part of simulation (drug-treated) with last time point from first
    ## part of simulation
    cell_pop = trajs[-1]
    ## Run second part of simulation and capture trajectories
    trajects = []
    for ind, cell in enumerate(cell_pop):
        print("%d" % ind)
        print("%d" % cell)
        t2 = np.linspace(0, 225, 226)  # in drug
        x2 = sim1.run(tspan=t2,
                      initials={model.species[0]: cell},
                      verbose=False)
        ## Put DIP rates into list (if not zero)
        ### Keep DIP rate if not zero, keep trajectory if end of first simulation > 50 cells
        if cell != 0:
            slope, intercept, r_value, p_value, std_err = sp.stats.linregress(t2, np.log2(
                x2.observables["Obs_Cell"] / x2.observables["Obs_Cell"][0]))
            if math.isnan(slope) == False:
                sim_dist = sim_dist + [slope]
        if cell > 50:
            trajects.append(np.log2(x2.observables["Obs_Cell"] / x2.observables["Obs_Cell"][0]))
    ## Capture p-value from K-S test comparison of experimental and simulated data
    D_stat, p_val = sp.ks_2samp(dat['DIP_Rate'], sim_dist)
    ## Return key elements
    return trajects, sim_dist, p_val

# Run function for sublines (except DS8 - PGM) at specified example rates
# and capture key data
DSs = ['PC9-DS1', 'PC9-DS3', 'PC9-DS4', 'PC9-DS6', 'PC9-DS7', 'PC9-DS9']
divs_a = [0.032, 0.030, 0.035, 0.030, 0.028, 0.025]
dths_a = [0.03050, 0.03075, 0.03175, 0.02940, 0.02625, 0.02375]
data_a = []
for i in range(len(DSs)):
    data_a.append(dist_compare(DSs[i], divs_a[i], dths_a[i]))

## Output model trajectories
trajectories_DS1 = pd.DataFrame(data_a[0][0])
trajectories_DS1 = trajectories_DS1.transpose()
trajectories_DS1.to_csv('trajectories_DS1_G50.csv')

trajectories_DS3 = pd.DataFrame(data_a[1][0])
trajectories_DS3 = trajectories_DS3.transpose()
trajectories_DS3.to_csv('trajectories_DS3_G50.csv')

trajectories_DS4 = pd.DataFrame(data_a[2][0])
trajectories_DS4 = trajectories_DS4.transpose()
trajectories_DS4.to_csv('trajectories_DS4_G50.csv')

trajectories_DS6 = pd.DataFrame(data_a[3][0])
trajectories_DS6 = trajectories_DS6.transpose()
trajectories_DS6.to_csv('trajectories_DS6_G50.csv')

trajectories_DS7 = pd.DataFrame(data_a[4][0])
trajectories_DS7 = trajectories_DS7.transpose()
trajectories_DS7.to_csv('trajectories_DS7_G50.csv')

trajectories_DS9 = pd.DataFrame(data_a[5][0])
trajectories_DS9 = trajectories_DS9.transpose()
trajectories_DS9.to_csv('trajectories_DS9_G50.csv')

## Output model distributions
distributions_0 = pd.DataFrame({'DS1': data_a[0][1]})
distributions_1 = pd.DataFrame({'DS3': data_a[1][1]})
distributions_2 = pd.DataFrame({'DS4': data_a[2][1]})
distributions_3 = pd.DataFrame({'DS6': data_a[3][1]})
distributions_4 = pd.DataFrame({'DS7': data_a[4][1]})
distributions_5 = pd.DataFrame({'DS9': data_a[5][1]})
distributions = pd.concat([distributions_0, distributions_1,
                           distributions_2, distributions_3,
                           distributions_4, distributions_5],
                          ignore_index=True, axis = 1)
distributions.columns = ['DS1', 'DS3', 'DS4', 'DS6', 'DS7', 'DS9']
distributions.to_csv('distributions_G50.csv')

## Output model p-values
pvalues = pd.DataFrame({'DS1': data_a[0][2], 'DS3': data_a[1][2],
                       'DS4': data_a[2][2], 'DS6': data_a[3][2],
                       'DS7': data_a[4][2], 'DS9': data_a[5][2]}, index=[0])
pvalues.to_csv('pvalues_G50.csv')
quit()