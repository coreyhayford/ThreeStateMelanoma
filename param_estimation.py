"""
Model-data overlay with fit score from Three-State Melanoma model and mixing experiment.
The parameters are estimated using the fmin() function in python.

The current parameter estimation only includes data from parental replicates from 
mixing experiment.

Expectation: Need to add data from additional mixing of subclones.

Based on code from: http://adventuresinpython.blogspot.com/2012/08/fitting-differential-equation-system-to.html

Implemented by: Corey E. Hayford (Vito Quaranta Lab, Vanderbilt University)

"""
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

# from three_state import model
#import three_state
import three_state_dip as three_state

### test

xlist = [1,2,3,4,5]

x_array = np.array(xlist)

ylist = [6,7,8,9,10]

y_array = np.array(ylist)

print xlist + ylist
print x_array + y_array

print (x_array + y_array) / 2

print 3/2

print [(0,None)]*10

#xtest = np.array([11,12,13])
#x_array[2:-1] = xtest
#print x_array


Model()  
    
three_state.declare_monomers()
three_state.declare_parameters()
three_state.declare_initial_conditions()
three_state.declare_observables()
three_state.declare_functions()

for p in model.parameters:
    print p

### Data - need nl2 and time ###
data = pd.read_csv('/Users/Corey/Documents/QuarantaLab/mixing_all', sep = "\t")
print data

data_parental = pd.read_csv('/Users/Corey/Documents/QuarantaLab/mixing_parental', sep = "\t")
print data_parental

data_parental_treated = data_parental.loc[data_parental['conc'] == 8]
print data_parental_treated


## Separate the 3 replicates into 3 datasets (lists) -- calculate the average of the lists and plot
print(range(0, len(data_parental_treated), 10))


sublist = [data_parental_treated[x:x+10] for x in range(0, len(data_parental_treated), 10)]

rep1 = sublist[0]
rep2 = sublist[1]
rep3 = sublist[2]



#time_all = [x for x in data_parental_treated["Time"]]
time = [x for x in rep1["Time"]] # Same time for all
#nl2 = [x for x in data_parental_treated["nl2"]]

nl2_rep1 = np.array([x for x in rep1["nl2"]])
nl2_rep2 = np.array([x for x in rep2["nl2"]])
nl2_rep3 = np.array([x for x in rep3["nl2"]])

nl2_average = (nl2_rep1 + nl2_rep2 + nl2_rep3) / 3.0

plt.plot(time, nl2_average, ms = 0.50, mfc = "0.25")

#def mean(all_lists):
#    return sum(all_lists) / len(all_lists)

#all_lists = [nl2_rep1, nl2_rep2, nl2_rep3]
#average_lists = map(mean, zip(*all_lists))

plt.plot(time, nl2_rep1, '*', ms = 12, mfc = "0.75")
plt.plot(time, nl2_rep2, '*', ms = 12, mfc = "0.75")
plt.plot(time, nl2_rep3, '*', ms = 12, mfc = "0.75")
#plt.plot(time, average_lists, ms = 12, mfc = "0.50")


time_sim = [time[0]]

for i in range(1, len(time)):
    if time[i] < time[i-1]:
        break
    time_sim.append(time[i])

#print time_sim

solver = Solver(model, time_sim, verbose = True)
solver.run() # saved in array - species y, obs in yobs
nl2_sim = np.log2(solver.yobs["Obs_All"]/solver.yobs["Obs_All"][0])
plt.plot(time_sim, nl2_sim, 'ro-', ms = 12, mfc = "None", mew = 2, mec = "r")

for i in range(len(time_sim)):
    print time[i], solver.yobs["Obs_All"][i]

#plt.show()

#exp_data = nl2


#for i in range(1, len(time)):
#    if data_parental_treated["rep"] == 1:
#        print 
#rep_1 = data_parental_treated["rep"] == 1
#print rep_1

#quit()


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

#define solver object
solver = Solver(model, time_sim)
# need solver yobs -- score function



### Score Fit of System ###
param_values = np.array([p.value for p in model.parameters])
def score(k_params): # Need initial conditions?
    #a.Get Solution to system
    #eq changed to solver.run
    k_params = np.array(k_params)
    param_values[3:] = k_params
    #param_values[3:-1] = k_params
    solver.run(param_values = param_values)
    nl2_sim = np.log2(solver.yobs["Obs_All"] / solver.yobs["Obs_All"][0]) 
    sse = ((nl2_average - nl2_sim)**2).sum() #need to be np array
    #b.Pick of Model Points to Compare - total population output of model
    ### model_output = np.log2("Obs_All")
    
    #c.Score Difference between model and data points
    # Sep the 3 reps into 3 datasets (lists) -- calculate the average of the lists and plot
    
    return sse


print model.parameters
#print model.parameters[3:-1]



print model.parameters[0]
print model.parameters[0].value # or name
#fit_score=score([p.value for p in model.parameters[3:-1]])
fit_score=score([p.value for p in model.parameters[3:]])
print fit_score

### Optimize Fit of parameters ###
#x0 = np.array([p.value for p in model.parameters])[3:-1]
x0 = np.array([p.value for p in model.parameters])[3:]
print x0
print len(x0)
#bounds = [(0,None)]*len(model.parameters[3:-1])
bounds = [(None, None), (0, 0), (None, None), (0, None), (0, None), (0, None), (0, None)]
print len(bounds)
result = minimize(score,x0, method = 'L-BFGS-B', bounds = bounds, options = {'maxiter':1000000})
#for i,p in enumerate(model.parameters[3:-1]):
#    print p,result.x[i]
for i,p in enumerate(model.parameters[3:]):
    print p,result.x[i]
#print xopt
#print type(xopt)


### New simulation

#param_values[3:-1] = result.x
param_values[3:] = result.x
solver.run(param_values = param_values)
nl2_sim = np.log2(solver.yobs["Obs_All"]/solver.yobs["Obs_All"][0])
plt.plot(time_sim, nl2_sim, 'go-', ms = 12, mfc = "None", mew = 2, mec = "g")




#plt.legend(('A','B','C','Experimental'),'upper center')
plt.xlabel('Time')
plt.ylabel('Population Doublings')

plt.show()