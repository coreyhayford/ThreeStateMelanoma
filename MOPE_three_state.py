"""
Multi-objective Parameter Estimation (MOPE) for Three-State Melanoma model and subclone mixing experiments.
The parameters are estimated using the minimize function in SciPy.optimize.

The current script includes 8 different experimental conditions: Parental, SC01 only, SC07 only, 
SC10 only, SC01/SC07 mix, SC01/SC10 mix, SC07/SC10 mix, and SC01/SC07/SC10 mix. These experimental
data are computationally fit to deterministic ordinary differential equations (ODEs) using a gradient
descent method to find the optimal parameter set that fits all 8 experimental conditions 
simultaneously. The best fit in local parameter space is achieved by minimizing over the sum of
a sum-of-squares cost function for all experimental conditions between the average of 3 experimental
replicates and the ODE values at experimental time points.

Expectation: Need higher time resolution in early time points and fewer at later time points to 
to balance out the weighting of behaviors in experiment - early dynamics take place much faster.

Next Steps: To create a multi-objective parameter estimation function in PySB

Implemented by: Corey E. Hayford and Leonard A. Harris (Vito Quaranta Lab, Vanderbilt University)

"""
 
### Import necessary python packages ###

import pandas as pd
from pylab import *
from numpy import *
from numpy.random import rand
import numpy as np
from scipy.optimize import minimize, fmin
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy import integrate
from pysb import *
from pysb.integrate import Solver

### Import Model ###
import three_state_dip as three_state

### Create new model -- call functions in previous model ###
Model()  
    
three_state.declare_monomers()
three_state.declare_parameters()
three_state.declare_initial_conditions()
three_state.declare_observables()
three_state.declare_functions()


####################################################################################################
### Experimental Data Cleaning and Plotting
####################################################################################################

# Import experimental data - need normalized log2 cell count and time
mixing_all = pd.read_csv('/Users/Corey/Documents/QuarantaLab/Mixing/mixing_all', sep = "\t")

# Open an empty list which will be appended with necessary experimental data
mixing = []

mixing.append(mixing_all.loc[mixing_all['CellLine'] == "SKMEL5 Parental"])
mixing.append(mixing_all.loc[mixing_all['CellLine'] == "SKMEL5 Subclone01"])
mixing.append(mixing_all.loc[mixing_all['CellLine'] == "SKMEL5 Subclone07"])
mixing.append(mixing_all.loc[mixing_all['CellLine'] == "SKMEL5 Subclone10"])
mixing.append(mixing_all.loc[mixing_all['CellLine'] == "Subclone01+Subclone07"])
mixing.append(mixing_all.loc[mixing_all['CellLine'] == "Subclone01+Subclone10"])
mixing.append(mixing_all.loc[mixing_all['CellLine'] == "Subclone07+Subclone10"])
mixing.append(mixing_all.loc[mixing_all['CellLine'] == "All"])


# Open an empty list of averaged (from 3 replicates) normalized log2 cell counts
nl2_average = []
time_points = 10 # Replace with however many time points in experimental data

# Loop over the mixing list and create sublists corresponding to the experimental replicates, and
# take the average, which will be appended to the empty list above. Loop for all conditions, and 
# plot the experimental replicates and average.
for i,m in enumerate(mixing):
    sublist = [m[x:x+time_points] for x in range(0, len(m), time_points)]
    rep1 = sublist[0] # first replicate data points
    print rep1
    rep2 = sublist[1]
    print rep2
    rep3 = sublist[2]
    print rep3
    # Create a time list for experimental and computational results -- same for all
    time = [x for x in rep1["Time"]] 
    # Make lists numpy arrays so can easily do math on them
    nl2_rep1 = np.array([x for x in rep1["nl2"]])
    nl2_rep2 = np.array([x for x in rep2["nl2"]])
    nl2_rep3 = np.array([x for x in rep3["nl2"]])

    nl2_average.append((nl2_rep1 + nl2_rep2 + nl2_rep3) / 3.0)
    
    plt.figure(i)
    
    plt.plot(time, nl2_average[-1], ms = 0.50, mfc = "0.25")
    
    plt.plot(time, nl2_rep1, '*', ms = 12, mfc = "0.75")
    plt.plot(time, nl2_rep2, '*', ms = 12, mfc = "0.75")
    plt.plot(time, nl2_rep3, '*', ms = 12, mfc = "0.75")  
plt.show()
    
nl2_average = np.array(nl2_average) # Make nl2_average list an array - easier to work with

### Define other initial conditions ###
# Function that creates an empty list that is appended with initial cell counts and returns
# that list. This list will be used to continuously change simulation ICs.
def ICs():
    ics = []
    ics.append([600, 1800, 600]) #parental
  #  ics.append([3000, 0, 0]) # all A 
    ics.append([0, 3000, 0]) # all B
    ics.append([0, 0, 3000]) # all C
  #  ics.append([1500, 1500, 0]) # AB
    ics.append([1500, 0, 1500]) # AC
    ics.append([0, 1500, 1500]) # BC
    ics.append([1000, 1000, 1000]) #ABC
    
    return ics

####################################################################################################
### Simulation
####################################################################################################

# Define solver object
solver = Solver(model, time, verbose = True)

# Loop over initial conditions to simulate ODEs with inital parameters and initial
# conditions and plot simulation results to experimental data plots
for i, ic in enumerate(ICs()):
    solver.run(y0 = ic) # ic is 3 element list
    # Making data normalized log2
    nl2_sim = np.log2(solver.yobs["Obs_All"]/solver.yobs["Obs_All"][0]) 
    plt.figure(i)
    plt.plot(time, nl2_sim, 'ro-', ms = 12, mfc = "None", mew = 2, mec = "r")
#plt.show()

####################################################################################################
### Parameter Estimation
####################################################################################################

### Score Function ###
param_values = np.array([p.value for p in model.parameters]) # Do first so can manipulate in function
def score(params):
    # k_value parameters (7) + Initial conditions parameters (3 * #conditions)
    # ex: params = [kdipA, kdipB, kdipB, kAB, kBA, kBC, kCB, 
    #          A_0[0], B_0[0], C_0[0], A_0[1], B_0[1], C_0[1], ...]
    params = np.array(params) # change list of params to array
    param_values = np.zeros(len(model.parameters)) # fill param_values object with an array of zeros (lenght of number of parameters)
    # ex: param_values = [A_0, B_0, C_0, kdipA, kdipB, kdipB, kAB, kBA, kBC, kCB]
    k_params = params[:7] # the first 7 parameters (index 0-6) are rate parameters
    param_values[3:] = k_params # replace 0s in param_values list with rate paramters (placed after ICs)
    init_params =  params[7:] # the final n*3 parameters - different ICs
    n_sims = len(init_params) / 3 # define number of simulations by number of ICs divided by 3 (ICs per simulation) - integer
    sse = 0 # set initial sum-of-squared errors to 0
    # Loop over each simulation to overwrite IC param_values with different ICs, olve ODE for those params, normalize, and return the 
    # SSE cost (which is additive for each sequential (new IC) simulation
    for n in range(n_sims):
        start = n*3
        end = start+3
        param_values[:3] = init_params[start:end]
        solver.run(param_values = param_values)
        nl2_sim = np.log2(solver.yobs["Obs_All"] / solver.yobs["Obs_All"][0]) 
        sse += ((nl2_average[n] - nl2_sim)**2).sum() # Sum of squared errors - cost function
    return sse

#fit_score=score([p.value for p in model.parameters[:]]) # No indecies reported so all parameters can float (including initial conditions)
#print fit_score

### Making parameter initial guess list ###
# Store rate parameter values in a list that is appended (at the end) with chosen ICs (in groups of 3) additively
x0 = [p.value for p in model.parameters[3:]]
for ic in ICs():
    x0 += ic

### Setting parameter bounds ###
# Store bounds in a list of rate parameters that is appended by initial IC guesses that can vary +/- bound_variance
# Some parameters are set to 0 - so the variance does not affect those
bounds = [(None, 0), (-0.5, 0.5), (0, None), (0, None), (0, None), (0, None), (0, None)]
bound_variance = 0.1 # 10% variance
for i in range(len(bounds), len(x0)):
    bounds.append((x0[i]-bound_variance*x0[i], x0[i]+bound_variance*x0[i]))
 
### Minimization of parameters ###
## Parameters (param_values) are minimized over the cost function (sse) using 
## an inital parameter guess (x0) with bounded constraints (bounds) and a max
## number of iterations -- FOR CONSTRAINED OPTIMIZATION
# Method L-BFGS-B works by identifying fixed and free variables at every step 
# (using a simple gradient method), and then using the L-BFGS method on the free
# variables only to get higher accuracy, and then repeating the process.
# See https://en.wikipedia.org/wiki/Limited-memory_BFGS#L-BFGS-B for more details
result = minimize(score,x0, method = 'L-BFGS-B', bounds = bounds, options = {'maxiter':1000000})

### New parameter determination (rates and ICs) ###
# Store rate parameters results in an empty list that is appended iteratively with the best fit parameters from that simulation
k_result = []
for i,p in enumerate(model.parameters[3:]):
    print p,result.x[i] # For visualization
    k_result.append(result.x[i])
    
# Store IC results in an empty list that is appended with an empty list for each iteration of ICs. The inner empty list is filled with the 
# IC values from the previous simulation (or initial for first - index + IC list (1 list) * 3 (items in list) + 7 (after all rate parameters).
# The previous result is then overwritten with the new parameters (using same equation above) - values are printed to screen
ic_result = []
for n in range(len(ICs())):
    ic_result.append([])
    print ICs()[n]," = [",
    for i,p in enumerate(model.parameters[:3]):
        print str(result.x[i+n*3+7]) + ",",
        ic_result[-1].append(result.x[i+n*3+7])
    print "]"

### Plotting parameter optimized model output ###
# Append a numpy array of 0s with rate parameters (last 7) and ICs (first 3). The rate parameters are the same for all the model outputs, 
# while the ICs change. In each iteration, the parameter values are saved, normalized, and plotted.
param_values = np.zeros(len(model.parameters))
param_values[3:] = k_result
for i, ic in enumerate(ICs()):
    param_values[:3] = ic_result[i] # N lists that contain 3 numbers
    solver.verbose = False
    solver.run(param_values=param_values) # saved in array - species y, obs in yobs, ic is 3 element list
    nl2_sim = np.log2(solver.yobs["Obs_All"]/solver.yobs["Obs_All"][0]) # making it nl2
    plt.figure(i)
    plt.plot(time, nl2_sim, 'go-', ms = 12, mfc = "None", mew = 2, mec = "g")
    plt.xlabel('Time')
    plt.ylabel('Population Doublings (nl2)')
    plt.title('Parameter Estimation')
plt.show()