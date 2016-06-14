"""
Multi-objective Parameter Estimation for Three-State Melanoma model and subclone mixing experiments.
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
mixing_all = pd.read_csv('/Users/Corey/Documents/QuarantaLab/mixing_all_rep', sep = "\t")

# Open an empty list which will be appended with necessary experimental data
mixing = []

mixing.append(mixing_all.loc[mixing_all['CellLine'] == "SKMEL5 Parental"])
#mixing.append(mixing_all.loc[mixing_all['CellLine'] == "SKMEL5 Subclone01"])
mixing.append(mixing_all.loc[mixing_all['CellLine'] == "SKMEL5 Subclone07"])
mixing.append(mixing_all.loc[mixing_all['CellLine'] == "SKMEL5 Subclone10"])
#mixing.append(mixing_all.loc[mixing_all['CellLine'] == "Subclone01+Subclone07"])
mixing.append(mixing_all.loc[mixing_all['CellLine'] == "Subclone01+Subclone10"])
mixing.append(mixing_all.loc[mixing_all['CellLine'] == "Subclone07+Subclone10"])
mixing.append(mixing_all.loc[mixing_all['CellLine'] == "All"])


# Open an empty list of averaged (from 3 replicates) normalized log2 cell counts
nl2_average = []
time_points = 17 # Replace with however many time points in experimental data

# Loop over the mixing list and create sublists corresponding to the experimental replicates, and
# take the average, which will be appended to the empty list above. Loop for all conditions, and 
# plot the experimental replicates and average.
for i,m in enumerate(mixing):
    sublist = [m[x:x+time_points] for x in range(0, len(m), time_points)]
    rep1 = sublist[0] # first replicate data points
    rep2 = sublist[1]
    rep3 = sublist[2]
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
    
    
nl2_average = np.array(nl2_average) # Make nl2_average list an array - easier to work with

### Define new initial conditions ###
# Function that creates an empty list that is appended with initial cell counts and returns
# that list. This list will be used to continuously change simulation ICs
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
   
### Score Function ###
param_values = np.array([p.value for p in model.parameters]) # Do first so can manipulate in function
#
def score(params):
    #13 parameters now
    #params = [kdipA, kdipB, kdipB, kAB, kBA, kBC, kCB, 
    #          A_0[0], B_0[0], C_0[0], A_0[1], B_0[1], C_0[1], ...]
    #eq changed to solver.run
    params = np.array(params)
    param_values = np.zeros(len(model.parameters))
    #param_values = [A_0, B_0, C_0, kdipA, kdipB, kdipB, kAB, kBA, kBC, kCB]
    k_params = params[:7]
    param_values[3:] = k_params
    init_params =  params[7:]
    n_sims = len(init_params) / 3 #integer
    sse = 0
    for n in range(n_sims):
        start = n*3
        end = start+3
        param_values[:3] = init_params[start:end]
        solver.run(param_values = param_values)
        nl2_sim = np.log2(solver.yobs["Obs_All"] / solver.yobs["Obs_All"][0]) 
        sse += ((nl2_average[n] - nl2_sim)**2).sum() # Sum of squared errors - cost function
    return sse

#print model.parameters
#print model.parameters[0].value # or name
#fit_score=score([p.value for p in model.parameters[3:-1]]) # this is for three_state (non-dip)
#fit_score=score([p.value for p in model.parameters[:]]) # No indecies reported so all parameters can float (including initial conditions)
#print fit_score

### Optimize Fit of parameters ###

#x0 = np.array([p.value for p in model.parameters])[3:-1] # this is for three_state (non-dip)
x0 = [p.value for p in model.parameters[3:]] # No indecies reported so all parameters can float (including initial conditions)
for ic in ICs():
    x0 += ic

print len(x0)
print x0


#
#print x0
#print len(x0)
#bounds = [(0,None)]*len(model.parameters[3:-1])
bounds = [(None, 0), (-0.5, 0.5), (0, None), (0, None), (0, None), (0, None), (0, None)]
for i in range(len(bounds), len(x0)):
    bounds.append((x0[i]-0.1*x0[i], x0[i]+0.1*x0[i]))
 # hard coded - should change
#print len(bounds)

# Method L-BFGS-B works by identifying fixed and free variables at every step 
# (using a simple gradient method), and then using the L-BFGS method on the free
# variables only to get higher accuracy, and then repeating the process.
# See https://en.wikipedia.org/wiki/Limited-memory_BFGS#L-BFGS-B for more details
# FOR CONSTRAINED OPTIMIZATION
result = minimize(score,x0, method = 'L-BFGS-B', bounds = bounds, options = {'maxiter':1000000})
#for i,p in enumerate(model.parameters[3:-1]): # this is for three_state (non-dip)
#    print p,result.x[i] # this is for three_state (non-dip)

k_result = []
for i,p in enumerate(model.parameters[3:]): # No indecies reported so all parameters can float (including initial conditions)
    print p,result.x[i]
    k_result.append(result.x[i])
    
# ic_result = []
# for i,p in enumerate(model.parameters[:3]):
#     ic_result.append([])
#     print p,
#     for n in range(len(ICs())):
#         print result.x[i+n*3+7],
#         ic_result[-1].append(result.x[i+n*3+7])
#     print

ic_result = []
for n in range(len(ICs())):
    ic_result.append([])
    print ICs()[n]," = [",
    for i,p in enumerate(model.parameters[:3]):
        print str(result.x[i+n*3+7]) + ",",
        ic_result[-1].append(result.x[i+n*3+7])
    print "]"

    


# xlist = [1,2,3,4,5]
# x_array = np.array(xlist)
# ylist = [6,7,8,9,10,11,12]
# y_array = np.array(ylist)
# 
# x_y_list = [x_array, y_array]
# print x_y_list
# print type(x_y_list)
# x_y_array = np.array(x_y_list)
# print x_y_array
# print type(x_y_array)

 
#print xopt
#print type(xopt)

### New simulation

#param_values[3:-1] = result.x
param_values = np.zeros(len(model.parameters))
param_values[3:] = k_result
for i, ic in enumerate(ICs()):
    param_values[:3] = ic_result[i] # 2 lists that contain 3 numbers
    solver.verbose = False
    solver.run(param_values=param_values) # saved in array - species y, obs in yobs, ic is 3 element list
    nl2_sim = np.log2(solver.yobs["Obs_All"]/solver.yobs["Obs_All"][0]) # making it nl2
    plt.figure(i)
    plt.plot(time, nl2_sim, 'go-', ms = 12, mfc = "None", mew = 2, mec = "g")


    plt.xlabel('Time')
    plt.ylabel('Population Doublings (nl2)')
    plt.title('Parameter Estimation')
plt.show()


### testing ideas ###


# print x_array + y_array
# print (x_array + y_array) / 2
# print 3/2
# print [(0,None)]*10


#def mean(all_lists):
#    return sum(all_lists) / len(all_lists)
#all_lists = [nl2_rep1, nl2_rep2, nl2_rep3]
#average_lists = map(mean, zip(*all_lists))
#plt.plot(time, average_lists, ms = 12, mfc = "0.50")

#for i in range(len(time)):
#    print time[i], solver.yobs["Obs_All"][i]



#sum_lists = [sum(i) for i in zip(*all_lists)]
#print all_lists
#print sum_lists
#array = np.array([[nl2_rep1], [nl2_rep2], [nl2_rep3]])
#rep_average = np.mean(array, axis = 0)
#print rep_average
#plt.plot(time, rep_average, ms = 12, mfc = "0.01")

#def sublist(list, n):
#    for i in range(1, len(list), n):
#        yield list[i: i+n]

#print sublist(time, 10)

#exp_data = np.array([[time[i], nl2[i]] for i in range(1, len(time)) if time[i] < 10])
#print exp_data

#for i in range(1, len(time), 10):
#    time_new = time[i, i+10],
#    print time_new

#def sublist(exp_data,10):
#    sub = []; result = []
#    for i in exp_data:
#        sub+=[i]
#        if len(sub)==n: result+=[sub]; sub =[]
#    if sub:result += [sub]
#    return result    
        
#cell_count = data_parental_treated["Cell.Nucleus"]
#print cell_count

# def eq(par,initial_cond,start_t,end_t,incr):
#     # Create time space
#     t  = np.linspace(start_t, end_t,incr)
#     # Create ODE system
#     def funct(y,t):
#        A=y[0]
#        B=y[1]
#        C=y[2]
#        par = k_div_A, k_div_B, k_div_C, k_death_A, k_death_B, k_death_C, k_AB, k_BA, k_CB, k_BC
#        # the model equations from PySB model
#        f0 = -A*k_AB - A*k_death_A + A*k_div_A + B*k_BA
#        f1 = A*k_AB - B*k_BA - B*k_BC - B*k_death_B + B*k_div_B + C*k_CB
#        f2 = B*k_BC - C*k_CB -C*k_death_C + C*k_div_C
#        return [f0, f1, f2]
     # Integrate over the equations
#     ds = integrate.odeint(funct,initial_cond,t)
#     print ds
#     return (ds[:,0],ds[:,1],ds[:,2],t)

#rates = (k_div_A, k_div_B, k_div_C, k_death_A, k_death_B, k_death_C, k_AB, k_BA, k_CB, k_BC)

#### Define Initial Conditions ###
#A_0 = 1000
#B_0 = 1000
#C_0 = 1000
#y_0 = [A_0, B_0, C_0]

#### Redefine model steps - encompassed in PySB model ###
#start_time=0.0
#end_time=120.0
#intervals=10
#mt=np.linspace(start_time,end_time,intervals)

#### Model index to compare to data ###
#findindex=lambda x:np.where(mt>=x)[0][0]
#mindex=map(findindex,time)
#print mindex

# time_sim = [time[0]]
# 
# for i in range(1, len(time)):
#     if time[i] < time[i-1]:
#         break
#     time_sim.append(time[i])

# solver.run # saved in array - species y, obs in yobs, 