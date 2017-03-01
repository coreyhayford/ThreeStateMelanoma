from pysb import *
from pysb.util import *
from pysb.macros import *
from pysb.bng import *
from pysb.integrate import odesolve
import pylab as pl
from numpy import linspace
from sympy import sympify
from scipy.stats import multivariate_normal
import numpy as np
import matplotlib.pyplot as plt

Model()

def declare_monomers():
    Monomer("A")
    Monomer("B")
    Monomer("C")
    Monomer("D")
    Monomer("E")
    Monomer("F")
    Monomer("G")
    Monomer("H")
    Monomer("I")
    Monomer("J")
    

def declare_parameters():
    Parameter("A_0", 1)
    Parameter("B_0", 1)
    Parameter("C_0", 1)
    Parameter("D_0", 1)
    Parameter("E_0", 1)
    Parameter("F_0", 1)
    Parameter("G_0", 1)
    Parameter("H_0", 1)
    Parameter("I_0", 1)
    Parameter("J_0", 1)

    Parameter("k_div_noDrug", 0.05)
    Parameter("k_death_noDrug", 0.025)
    
    alias_model_components()

def declare_ICs():
    Initial(A, A_0)
    Initial(B, B_0)
    Initial(C, C_0)
    Initial(D, D_0)
    Initial(E, E_0)
    Initial(F, F_0)
    Initial(G, G_0)
    Initial(H, H_0)
    Initial(I, I_0)
    Initial(J, J_0)
    
def declare_observables():
    Observable("Obs_A", A)
    Observable("Obs_B", B)
    Observable("Obs_C", C)
    Observable("Obs_D", D)
    Observable("Obs_E", E)
    Observable("Obs_F", F)
    Observable("Obs_G", G)
    Observable("Obs_H", H)
    Observable("Obs_I", I)
    Observable("Obs_J", J)
    Observable("Obs_All", A()+B()+C()+D()+E()+F()+G()+H()+I()+J())

def declare_rules():
    Rule('div_A', A() >> A() + A(), k_div_noDrug)
    Rule('div_B', B() >> B() + B(), k_div_noDrug)
    Rule('div_C', C() >> C() + C(), k_div_noDrug)
    Rule('div_D', D() >> D() + D(), k_div_noDrug)
    Rule('div_E', E() >> E() + E(), k_div_noDrug)
    Rule('div_F', F() >> F() + F(), k_div_noDrug)
    Rule('div_G', G() >> G() + G(), k_div_noDrug)
    Rule('div_H', H() >> H() + H(), k_div_noDrug)
    Rule('div_I', I() >> I() + I(), k_div_noDrug)
    Rule('div_J', J() >> J() + J(), k_div_noDrug)
    
    Rule('death_A', A() >> None, k_death_noDrug)
    Rule('death_B', B() >> None, k_death_noDrug)
    Rule('death_C', C() >> None, k_death_noDrug)
    Rule('death_D', D() >> None, k_death_noDrug)
    Rule('death_E', E() >> None, k_death_noDrug)
    Rule('death_F', F() >> None, k_death_noDrug)
    Rule('death_G', G() >> None, k_death_noDrug)
    Rule('death_H', H() >> None, k_death_noDrug)
    Rule('death_I', I() >> None, k_death_noDrug)
    Rule('death_J', J() >> None, k_death_noDrug)
    

# plt.figure("Pre-drug Stochastic Simulations")
####

def plot(t, t_start=0, label=False):
    plt.plot(t+t_start, y["Obs_A"], 'c-', lw=2, 
             label= ('A DIP = %g' %(k_divide_A.value-k_death_A.value)) if label else None)
    plt.plot(t+t_start, y["Obs_B"], 'b-', lw=2, 
             label= ('B DIP = %g' %(k_divide_B.value-k_death_B.value)) if label else None)
    plt.plot(t+t_start, y["Obs_C"], 'g-', lw=2, 
             label= ('C DIP = %g' %(k_divide_C.value-k_death_C.value)) if label else None)
    plt.plot(t+t_start, y["Obs_D"], 'r-', lw=2, 
             label= ('D DIP = %g' %(k_divide_D.value-k_death_D.value)) if label else None)
    plt.plot(t+t_start, y["Obs_E"], 'm-', lw=2, 
             label= ('E DIP = %g' %(k_divide_E.value-k_death_E.value)) if label else None)
    plt.plot(t+t_start, y["Obs_F"], 'c:', lw=2, 
             label= ('F DIP = %g' %(k_divide_F.value-k_death_F.value)) if label else None)
    plt.plot(t+t_start, y["Obs_G"], 'b:', lw=2, 
             label= ('G DIP = %g' %(k_divide_G.value-k_death_G.value)) if label else None)
    plt.plot(t+t_start, y["Obs_H"], 'g:', lw=2, 
             label= ('H DIP = %g' %(k_divide_H.value-k_death_H.value)) if label else None)
    plt.plot(t+t_start, y["Obs_I"], 'r:', lw=2, 
             label= ('I DIP = %g' %(k_divide_I.value-k_death_I.value)) if label else None)
    plt.plot(t+t_start, y["Obs_J"], 'm:', lw=2, 
             label= ('J DIP = %g' %(k_divide_J.value-k_death_J.value)) if label else None)
    
#     plt.plot(t+t_start, y["Obs_All"], 'k-', lw=4, 
#              label = ('All') if label else None)
    plt.xlabel("Time")
    plt.ylabel("Population Doublings")
    plt.legend(loc = 0) #, fontsize = 6)


# def plot(t, t_start=0, label=False):
#     plt.plot(t+t_start, np.log2(y["Obs_A"]), 'c-', lw=2, 
#              label= ('A DIP = %g' %(k_divide_A.value-k_death_A.value)) if label else None)
#     plt.plot(t+t_start, np.log2(y["Obs_B"]), 'b-', lw=2, 
#              label= ('B DIP = %g' %(k_divide_B.value-k_death_B.value)) if label else None)
#     plt.plot(t+t_start, np.log2(y["Obs_C"]), 'g-', lw=2, 
#              label= ('C DIP = %g' %(k_divide_C.value-k_death_C.value)) if label else None)
#     plt.plot(t+t_start, np.log2(y["Obs_D"]), 'r-', lw=2, 
#              label= ('D DIP = %g' %(k_divide_D.value-k_death_D.value)) if label else None)
#     plt.plot(t+t_start, np.log2(y["Obs_E"]), 'm-', lw=2, 
#              label= ('E DIP = %g' %(k_divide_E.value-k_death_E.value)) if label else None)
#     plt.plot(t+t_start, np.log2(y["Obs_F"]), 'c:', lw=2, 
#              label= ('F DIP = %g' %(k_divide_F.value-k_death_F.value)) if label else None)
#     plt.plot(t+t_start, np.log2(y["Obs_G"]), 'b:', lw=2, 
#              label= ('G DIP = %g' %(k_divide_G.value-k_death_G.value)) if label else None)
#     plt.plot(t+t_start, np.log2(y["Obs_H"]), 'g:', lw=2, 
#              label= ('H DIP = %g' %(k_divide_H.value-k_death_H.value)) if label else None)
#     plt.plot(t+t_start, np.log2(y["Obs_I"]), 'r:', lw=2, 
#              label= ('I DIP = %g' %(k_divide_I.value-k_death_I.value)) if label else None)
#     plt.plot(t+t_start, np.log2(y["Obs_J"]), 'm:', lw=2, 
#              label= ('J DIP = %g' %(k_divide_J.value-k_death_J.value)) if label else None)
#     
#     plt.plot(t+t_start, np.log2(y["Obs_All"]), 'k-', lw=4, 
#              label = ('All') if label else None)
#     plt.xlabel("Time")
#     plt.ylabel("Population Doublings")
#     plt.legend(loc = 0) #, fontsize = 6)

for i in range(1,2):
    
    Model()  
    declare_monomers()
    declare_parameters()
    declare_ICs()
    declare_observables()
    declare_rules()
    
#     for m in model.monomers:
#         print m
#     quit()
    
    t = linspace(0, 200, 201)
    
    #y = odesolve(model, t, verbose=True)
    y = run_ssa(model,t_end = t[-1], n_steps = len(t)-1, verbose=True)
    
    
    plot(t, label = False)
    
    plt.axvline(x=t[-1], color='red', ls='--', lw=4)

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

#plt.tight_layout()
# plt.show()
# quit()
# print(y["Obs_A"][-1])
# print(y["Obs_B"][-1])
# print(y["Obs_C"][-1])
# print(y["Obs_D"][-1])
# print(y["Obs_E"][-1])
# print(y["Obs_F"][-1])
# print(y["Obs_G"][-1])
# print(y["Obs_H"][-1])
# print(y["Obs_I"][-1])
# print(y["Obs_J"][-1])
# 
# print(y["Obs_A"][-1]/y["Obs_All"][-1])
# print(y["Obs_B"][-1]/y["Obs_All"][-1])
# print(y["Obs_C"][-1]/y["Obs_All"][-1])
# print(y["Obs_D"][-1]/y["Obs_All"][-1])
# print(y["Obs_E"][-1]/y["Obs_All"][-1])
# print(y["Obs_F"][-1]/y["Obs_All"][-1])
# print(y["Obs_G"][-1]/y["Obs_All"][-1])
# print(y["Obs_H"][-1]/y["Obs_All"][-1])
# print(y["Obs_I"][-1]/y["Obs_All"][-1])
# print(y["Obs_J"][-1]/y["Obs_All"][-1])

# These will be ICs for Post Drug...
def declare_newparameters():
    Parameter("postDrugIC_A", y["Obs_A"][-1])
    Parameter("postDrugIC_B", y["Obs_B"][-1])
    Parameter("postDrugIC_C", y["Obs_C"][-1])
    Parameter("postDrugIC_D", y["Obs_D"][-1])
    Parameter("postDrugIC_E", y["Obs_E"][-1])
    Parameter("postDrugIC_F", y["Obs_F"][-1])
    Parameter("postDrugIC_G", y["Obs_G"][-1])
    Parameter("postDrugIC_H", y["Obs_H"][-1])
    Parameter("postDrugIC_I", y["Obs_I"][-1])
    Parameter("postDrugIC_J", y["Obs_J"][-1])
 
    Parameter("postDrugIC_All", y["Obs_A"][-1] + y["Obs_B"][-1] + y["Obs_C"][-1] + y["Obs_D"][-1] + y["Obs_E"][-1] + y["Obs_F"][-1] + y["Obs_G"][-1] + y["Obs_H"][-1] + y["Obs_I"][-1] + y["Obs_J"][-1])
    
    alias_model_components()


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


# import preDrugDynamics
# 
# import decimal
# 
# def drange(x, y, jump):
#   while x < y:
#     yield float(x)
#     x += decimal.Decimal(jump)


# print(model.monomers)

def declare_drug_parameters():
    
    Parameter("k_divide_A", 0.35)
    Parameter("k_divide_B", 0.35)
    Parameter("k_divide_C", 0.35)
    Parameter("k_divide_D", 0.35)
    Parameter("k_divide_E", 0.35)
    Parameter("k_divide_F", 0.35)
    Parameter("k_divide_G", 0.35)
    Parameter("k_divide_H", 0.35)
    Parameter("k_divide_I", 0.35)
    Parameter("k_divide_J", 0.35)
    
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

declare_drug_parameters()
declare_newparameters()

# print(model.parameters)


A_0.value = y["Obs_A"][-1]
B_0.value = y["Obs_B"][-1]
C_0.value = y["Obs_C"][-1]
D_0.value = y["Obs_D"][-1]
E_0.value = y["Obs_E"][-1]
F_0.value = y["Obs_F"][-1]
G_0.value = y["Obs_G"][-1]
H_0.value = y["Obs_H"][-1]
I_0.value = y["Obs_I"][-1]
J_0.value = y["Obs_J"][-1]

### each monomer is not globally defined
# def declare_drug_ICs():
#     Initial(A, postDrugIC_A)
#     Initial(B, postDrugIC_B)
#     Initial(C, postDrugIC_C)
#     Initial(D, postDrugIC_D)
#     Initial(E, postDrugIC_E)
#     Initial(F, postDrugIC_F)
#     Initial(G, postDrugIC_G)
#     Initial(H, postDrugIC_H)
#     Initial(I, postDrugIC_I)
#     Initial(J, postDrugIC_J)
# 
# declare_drug_ICs()
# print(model.initial_conditions)

# preDrugDynamics.declare_observables()
# print(model.observables)

def declare_drug_rules():
    Rule('divide_A', A() >> A() + A(), k_divide_A)
    Rule('divide_B', B() >> B() + B(), k_divide_B)
    Rule('divide_C', C() >> C() + C(), k_divide_C)
    Rule('divide_D', D() >> D() + D(), k_divide_D)
    Rule('divide_E', E() >> E() + E(), k_divide_E)
    Rule('divide_F', F() >> F() + F(), k_divide_F)
    Rule('divide_G', G() >> G() + G(), k_divide_G)
    Rule('divide_H', H() >> H() + H(), k_divide_H)
    Rule('divide_I', I() >> I() + I(), k_divide_I)
    Rule('divide_J', J() >> J() + J(), k_divide_J)
    
    Rule('die_A', A() >> None, k_death_A)
    Rule('die_B', B() >> None, k_death_B)
    Rule('die_C', C() >> None, k_death_C)
    Rule('die_D', D() >> None, k_death_D)
    Rule('die_E', E() >> None, k_death_E)
    Rule('die_F', F() >> None, k_death_F)
    Rule('die_G', G() >> None, k_death_G)
    Rule('die_H', H() >> None, k_death_H)
    Rule('die_I', I() >> None, k_death_I)
    Rule('die_J', J() >> None, k_death_J)
 
declare_drug_rules()
# print(model.rules)   

    
## DONE: established initial conditions as a function of gillespie sims
## DONE: Pick DIP rates from multinomial distribution
## DONE: Plot before and after drug simulated results
## DONE: Apply new plot function to stitch together results

## GOAL: Put in line for drug introduction
## GOAL: Play with different DIP rates/distributions to model diversification phenotype
## GOAL: Only pick time points similar to experiment??
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

#     t = linspace(0, 200, 201)
    
#     y = odesolve(model, t, verbose=True)
    y = run_ssa(model,t_end = t[-1], n_steps = len(t)-1, verbose=True)
    plot(t, t_start=t[-1], label=True)
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



# declare_newparameters()    
# print(model.parameters)

# def plot(t, t_start=0, label=False):
#     plt.plot(t+t_start, np.log2(y["Obs_A"]/A_init), 'c-', lw=2, 
#              label= ('A DIP = %g, IC = %d (%d)' %(k_div_A.value-k_death_A.value, A0.value, picks[0])) if label else None)
#     plt.plot(t+t_start, np.log2(y["Obs_B"]/B_init), 'r-', lw=2, 
#              label = ('B DIP = %g, IC = %d (%d)' %(k_div_B.value-k_death_B.value, B0.value, picks[1])) if label else None)
#     plt.plot(t+t_start, np.log2(y["Obs_C"]/C_init), 'g-', lw=2, 
#              label = ('C DIP = %g, IC = %d (%d)' %(k_div_C.value-k_death_C.value, C0.value, picks[2])) if label else None)
#     plt.plot(t+t_start, np.log2(y["Obs_All"]/(A_init+B_init+C_init)), 'k-', lw=4, 
#              label = ('All'  if label else None))
#     plt.xlabel("Time")
#     plt.ylabel("Population Doublings")
#     plt.legend(loc = 0) #, fontsize = 6)


# import numpy as np
# 
# a = np.array([1,2,3,4,5,6,7,8,9,10,11])
# print (a)
# b = []
# for i in range(1,len(a+1)):
#     print(i)
#     print(a[i])
#     print(a[i-1])
#     print((a[i]+a[i-1])/2.)
#     b.append((a[i]+a[i-1])/2.)
# c = np.array(b)
# print(b)
# print(c)
