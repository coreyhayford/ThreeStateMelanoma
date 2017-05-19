from pysb import *
from pysb.util import *
from pysb.macros import *
from pysb.bng import *
from pysb.integrate import odesolve
import pylab as pl
from numpy import linspace
from numpy import arange
from sympy import sympify
import scipy
from scipy.stats import multivariate_normal
from scipy.stats import norm
import numpy as np
import matplotlib.pyplot as plt

#fig, ax = plt.subplots(1, 1)

dips = np.random.normal(-0.033, 0.01, 100)
# print(dips)
count, bins, ignored = plt.hist(dips, 10, normed=True)
# plt.show()
# print(count)
# print(bins)

# print(sum(count))

count_normalized = count/sum(count)
# print(count_normalized)
# print(sum(count_normalized))

bins_avg = []
for i in range(1,len(bins+1)):
#     print(i)
#     print(a[i])
#     print(a[i-1])
#     print((bins[i]+bins[i-1])/2.)
    bins_avg.append((bins[i]+bins[i-1])/2.)
bins_avg_array = np.array(bins_avg)
    # print(bins_avg)
# print(bins_avg_array)


## Number of cells in each state - NOT necessary is explicitly modeled
picks = np.random.multinomial(10, count_normalized)
# print(picks)


import preDrugDynamics
# 
# import decimal
# 
# def drange(x, y, jump):
#   while x < y:
#     yield float(x)
#     x += decimal.Decimal(jump)

Model()

preDrugDynamics.declare_monomers()
# print(model.monomers)

def declare_parameters():
    
    Parameter("k_div_A", 0.35)
    Parameter("k_div_B", 0.35)
    Parameter("k_div_C", 0.35)
    Parameter("k_div_D", 0.35)
    Parameter("k_div_E", 0.35)
    Parameter("k_div_F", 0.35)
    Parameter("k_div_G", 0.35)
    Parameter("k_div_H", 0.35)
    Parameter("k_div_I", 0.35)
    Parameter("k_div_J", 0.35)
    
    Parameter("k_death_A", -bins_avg_array[0]+0.35)
    Parameter("k_death_B", -bins_avg_array[1]+0.35)
    Parameter("k_death_C", -bins_avg_array[2]+0.35)
    Parameter("k_death_D", -bins_avg_array[3]+0.35)
    Parameter("k_death_E", -bins_avg_array[4]+0.35)
    Parameter("k_death_F", -bins_avg_array[5]+0.35)
    Parameter("k_death_G", -bins_avg_array[6]+0.35)
    Parameter("k_death_H", -bins_avg_array[7]+0.35)
    Parameter("k_death_I", -bins_avg_array[8]+0.35)
    Parameter("k_death_J", -bins_avg_array[9]+0.35)

declare_parameters()
preDrugDynamics.declare_newparameters()

# print(model.parameters)


### each monomer is not globally defined
def declare_ICs():
    Initial(A, postDrugIC_A)
    Initial(B, postDrugIC_B)
    Initial(C, postDrugIC_C)
    Initial(D, postDrugIC_D)
    Initial(E, postDrugIC_E)
    Initial(F, postDrugIC_F)
    Initial(G, postDrugIC_G)
    Initial(H, postDrugIC_H)
    Initial(I, postDrugIC_I)
    Initial(J, postDrugIC_J)

declare_ICs()
# print(model.initial_conditions)

preDrugDynamics.declare_observables()
# print(model.observables)

def declare_rules():
    Rule('div_A', A() >> A() + A(), k_div_A)
    Rule('div_B', B() >> B() + B(), k_div_B)
    Rule('div_C', C() >> C() + C(), k_div_C)
    Rule('div_D', D() >> D() + D(), k_div_D)
    Rule('div_E', E() >> E() + E(), k_div_E)
    Rule('div_F', F() >> F() + F(), k_div_F)
    Rule('div_G', G() >> G() + G(), k_div_G)
    Rule('div_H', H() >> H() + H(), k_div_H)
    Rule('div_I', I() >> I() + I(), k_div_I)
    Rule('div_J', J() >> J() + J(), k_div_J)
    
    Rule('death_A', A() >> None, k_death_A)
    Rule('death_B', B() >> None, k_death_B)
    Rule('death_C', C() >> None, k_death_C)
    Rule('death_D', D() >> None, k_death_D)
    Rule('death_E', E() >> None, k_death_E)
    Rule('death_F', F() >> None, k_death_F)
    Rule('death_G', G() >> None, k_death_G)
    Rule('death_H', H() >> None, k_death_H)
    Rule('death_I', I() >> None, k_death_I)
    Rule('death_J', J() >> None, k_death_J)
 
declare_rules()
# print(model.rules)   

    
## DONE: established initial conditions as a function of gillespie sims
## DONE: Pick DIP rates from multinomial distribution
## DONE: Plot before and after drug simulated results

## GOAL: Apply new plot function to stitch together results
## GOAL: Put in line for drug introduction
## GOAL: Play with different DIP rates/distributions to model diversification phenotype

## GOAL: Include state transitions in the model

## GOAL: Add more barcodes and plot like Hata


# plt.figure("Post-drug Stochastic Simulations")
####
for i in range(1,2):
    
#     Model()  
#     declare_parameters()
#     declare_ICs()
#     declare_observables()
#     declare_rules()

    t = linspace(0, 200, 201)
    
#     y = odesolve(model, t, verbose=True)
    y = run_ssa(model,t_end = t[-1], n_steps = len(t)-1, verbose=True)
    preDrugDynamics.plot(t, t_start=t[-1], label=True)
#     plot(t, t_start=t[-1], label=True)
    #plt.subplot(2,2,i)
#     plt.plot(t, y["Obs_A"], 'r-', lw=2, label = "A")
#     plt.plot(t, y["Obs_B"], 'b-', lw=2, label = "B")
#     plt.plot(t, y["Obs_C"], 'g-', lw=2, label = "C")
#     plt.plot(t, y["Obs_D"], 'm-', lw=2, label = "D")
#     plt.plot(t, y["Obs_E"], 'c-', lw=2, label = "E")
#     plt.plot(t, y["Obs_F"], 'r:', lw=2, label = "F")
#     plt.plot(t, y["Obs_G"], 'b:', lw=2, label = "G")
#     plt.plot(t, y["Obs_H"], 'g:', lw=2, label = "H")
#     plt.plot(t, y["Obs_I"], 'm:', lw=2, label = "I")
#     plt.plot(t, y["Obs_J"], 'c:', lw=2, label = "J")
#     plt.plot(t, y["Obs_All"], 'k-', lw=4, label = "All")
#     plt.xlabel("Time")
#     plt.ylabel("Number of Barcodes")
#     plt.legend(loc = 0)
#     plt.title("Replicate " + str(i))
# print(model.parameters)
plt.show()
# x = np.linspace(-0.033, 0.033, 10)
# # x = np.linspace(norm.ppf(0.01), norm.ppf(0.1), 10)
# print(x)
# print(norm.pdf(x))
# ax.plot(x, norm.pdf(x), 'r-', lw=5, alpha=0.6, label='norm pdf')
# # plt.show()
# # quit()
# 
# def declare_parameters():
#     
#     dips = np.random.poisson(0, 10)
#     count, bins, ignored = plt.hist(dips, 20, normed=True)
#     probabilities = np.random.multinomial(dips)
#     print(probabilities)
# #     mean = 0.25
# #     var = 1
# #     d_range = np.array(list(drange(0, 2*mean, '0.05')))
# #     print(d_range)
# #     pmf_dist = scipy.stats.poisson.pmf(x, mu = 0)
# #     print(pmf_dist)
# #     quit()
# #     print(len(pmf_dist))
# #     print(np.sum(pmf_dist))
#     #print(np.array(range(mean+1)))
#     plt.bar(left = np.array(arange(2*mean)), height = pmf_dist, width = 0.01)
# declare_parameters()
# plt.show()
# quit()
#     x_nvals = np.linspace(0, 3*mean+1, 3*mean*100)
#     pdf_dist = scipy.stats.norm.pdf(x_nvals, loc = mean, scale = np.sqrt(var)) #, endpoint, retstep, dtype))
#     #a = np.random.multinomial(int(model.parameters["postDrugIC_All"].value), pdf_dist)
# #     print(model.parameters["postDrugIC_All"].value)
#     print(pdf_dist)
#     b = pdf_dist/pdf_dist.sum()
#     print(b)
#     print(pdf_dist.sum())
#     #print(np.sum(pdf_dist)*(x_nvals[1]-x_nvals[0]))
#     a = np.random.multinomial(10, pdf_dist)
#     a1 = a/10.
#     print(a)
#     print(a1)
#     #print(counts)
#     plt.plot(range(len(pdf_dist)), a1, "g-o")


