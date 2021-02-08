from pysb import *
from pysb.util import *
from pysb.macros import *
from pysb.bng import *
from pysb.integrate import odesolve
import pylab as pl
from numpy import linspace
from sympy import sympify
import numpy as np
import random
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats


Model()

# Randomly picking the initial conditions - "Biased Die" Method
  
def roll(Dist, num):
     
    randRolls = []
    results = []
    for i in range(num):
        randRoll = random.random()
        #print(randRoll)
        randRolls.append(randRoll) # in [0,1)
        sum = 0
        result = 0 #bin
        for val in Dist: # bin, val in enumerate(Dist)
            sum += val
            #print("Mass: " + str(val))
             
            if randRoll < sum:
                #print(result)
                results.append(result) #[bin]
                break 
            result += 1 # delete 
            #print("Sum " + str(sum))
            #print(result)
 
    #print(randRolls)
    #print(results)
    counts = []
    for j in range(len(Dist)):
        count = 0
        for res in results:
            if res == j:
                count += 1
        counts.append(count)
    counts = 1.* np.array(counts)/num #number. = float
    #print(counts)
    return counts
  
#     ps = pd.Series(results)
#     global counts
#     counts = ps.value_counts(sort = False)
#     print(counts)
plt.figure(1)
mean = 10
var = mean
x_vals = np.linspace(0,3*mean+1, 3*mean)
x_nvals = np.linspace(0, 3*mean+1, 3*mean*100)
# x_vals = np.linspace(-3*var+mean, 3.001*var+mean, 300)
pdf_dist = scipy.stats.norm.pdf(x_nvals, loc = mean, scale = np.sqrt(var)) #, endpoint, retstep, dtype))
print(np.sum(pdf_dist)*(x_nvals[1]-x_nvals[0]))

pmf_dist = scipy.stats.poisson.pmf(range(2*mean), mu = mean)
print(pmf_dist)
print(len(pmf_dist))
print(np.sum(pmf_dist))
print(np.array(range(mean+1)))
quit()


plt.plot(x_nvals, pdf_dist, "r-", lw = 2)
plt.bar(left = np.array(range(2*mean))-0.5, height = pmf_dist, width = 1)
#plt.xlim(xmin = 50, xmax = 150)

# plt.figure(2)
# sampleDist = (0.02, 0.04, 0.08, 0.16, 0.20, 0.20, 0.16, 0.08, 0.04, 0.02)
# plt.plot(range(len(sampleDist)), sampleDist, "r-*", lw = 2)

#counts = roll(pmf_dist, 10000)

#plt.bar(left = np.array(range(len(counts)))-0.5, height = counts, width = 1) #, bottom, hold, data) 

a = np.random.multinomial(10000, pmf_dist)
a1 = a/10000.
print(a)
print(a1)
#print(counts)
plt.plot(range(len(pmf_dist)), a1, "g-o")
plt.show()
#plt.show()

#print(results)
#print(counts)
quit()
def declare_monomers():
    
    Monomer('Cell', ['barcode', 'dip'], {'barcode':['1', '2', '3', '4', '5', '6', '7', '8', '9', '10'], 
                                         'dip':['negA', 'negB', 'negC', 'negD', 'negE', 'posE', 
                                                'posD','posC', 'posB', 'posA']})

def declare_parameters():
    
### Initial Conditions - number of cells in each DIP rate bin


    for i,j in enumerate(counts):
#     print("I: " + str(i))
#     print("J: " + str(j))
        Parameter("Cell%d_0" %i, j)

#     [Parameter("Cell%d_0" %i, ) for i in range(11)]




#     Parameter("CellA_0", 10)
#     Parameter("CellB_0", 10)
#     Parameter("CellC_0", 10)
#     Parameter("CellD_0", 10)
#     Parameter("CellE_0", 10)
#     Parameter("CellF_0", 10)
#     Parameter("CellG_0", 10)
#     Parameter("CellH_0", 10)
#     Parameter("CellI_0", 10)
#     Parameter("CellJ_0", 10)
    
### Rate parameters - Non overlapping states (small SD) and 10 different DIP rates

    Parameter("k_dip_A", 0.02 * np.random.randn() - 0.25)
    Parameter("k_dip_B", 0.02 * np.random.randn() - 0.2)
    Parameter("k_dip_C", 0.02 * np.random.randn() - 0.15)
    Parameter("k_dip_D", 0.02 * np.random.randn() - 0.1)
    Parameter("k_dip_E", 0.02 * np.random.randn() - 0.05)
    Parameter("k_dip_F", 0.02 * np.random.randn() + 0.05)
    Parameter("k_dip_G", 0.02 * np.random.randn() + 0.1)
    Parameter("k_dip_H", 0.02 * np.random.randn() + 0.15)
    Parameter("k_dip_I", 0.02 * np.random.randn() + 0.2)
    Parameter("k_dip_J", 0.02 * np.random.randn() + 0.25)
    
    alias_model_components()
    
declare_parameters()
print(model.parameters)


def declare_initial_conditions():
    # Have to define barcode and dip for each IC
    ## Should I randomly distribute those into different DIP rates?
    Initial(Cell(barcode = '1', dip = 'negA'), Cell0_0)
    Initial(Cell(barcode = '2'), Cell1_0)
    Initial(Cell(barcode = '3'), Cell2_0)
    Initial(Cell(barcode = '4'), Cell3_0)
    Initial(Cell(barcode = '5'), Cell4_0)
    Initial(Cell(barcode = '6'), Cell5_0)
    Initial(Cell(barcode = '7'), Cell6_0)
    Initial(Cell(barcode = '8'), Cell7_0)
    Initial(Cell(barcode = '9'), Cell8_0)
    Initial(Cell(barcode = '10'), Cell9_0)

def declare_observables():
    
    Observable("Obs_A", Cell(barcode = '1'))
    Observable("Obs_B", Cell(barcode = '2'))
    Observable("Obs_C", Cell(barcode = '3'))
    Observable("Obs_D", Cell(barcode = '4'))
    Observable("Obs_E", Cell(barcode = '5'))
    Observable("Obs_F", Cell(barcode = '6'))
    Observable("Obs_G", Cell(barcode = '7'))
    Observable("Obs_H", Cell(barcode = '8'))
    Observable("Obs_I", Cell(barcode = '9'))
    Observable("Obs_J", Cell(barcode = '10'))
    Observable("Obs_All", Cell())
    
def declare_functions():
    
#     barcodes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10']
#     dips = ['negA', 'negB', 'negC', 'negD', 'negE', 'posE', 'posD','posC', 'posB', 'posA']
#     names = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']
#     for name in names: 
#         for bar in barcodes:
#             for d in dips:
#                 Rule('DIP_%s' %name, Cell(barcode = bar, dip = d) >> Cell(barcode = bar, dip = d) +
#                     Cell(barcode = bar, dip = d), 'k_dip_%s' %name)
    Rule('DIP_A', Cell(dip = 'negA') >> Cell(dip = 'negA') + 
         Cell(dip = 'negA'), k_dip_A)
    Rule('DIP_B', Cell(dip = 'negB') >> Cell(dip = 'negB') + 
         Cell(dip = 'negB'), k_dip_B)
    Rule('DIP_C', Cell(dip = 'negC') >> Cell(dip = 'negC') + 
         Cell(dip = 'negC'), k_dip_C)
    Rule('DIP_D', Cell(dip = 'negD') >> Cell(dip = 'negD') + 
         Cell(dip = 'negD'), k_dip_D)
    Rule('DIP_E', Cell(dip = 'negE') >> Cell(dip = 'negE') + 
         Cell(dip = 'negE'), k_dip_E)
    Rule('DIP_F', Cell(dip = 'posE') >> Cell(dip = 'posE') + 
         Cell(dip = 'posE'), k_dip_F)
    Rule('DIP_G', Cell(dip = 'posD') >> Cell(dip = 'posD') + 
         Cell(dip = 'posD'), k_dip_G)
    Rule('DIP_H', Cell(dip = 'posC') >> Cell(dip = 'posC') + 
         Cell(dip = 'posC'), k_dip_H)
    Rule('DIP_I', Cell(dip = 'posB') >> Cell(dip = 'posB') + 
         Cell(dip = 'posB'), k_dip_I)
    Rule('DIP_J', Cell(dip = 'posA') >> Cell('10', dip = 'posA') + 
         Cell(dip = 'posA'), k_dip_J)
    
#     Rule('DIP_A', Cell(barcode = '1', dip = 'negA') >> Cell(barcode = '1', dip = 'negA') + 
#          Cell(barcode = '1', dip = 'negA'), k_dip_A)
#     Rule('DIP_B', Cell(barcode = '2', dip = 'negB') >> Cell(barcode = '2', dip = 'negB') + 
#          Cell(barcode = '2', dip = 'negB'), k_dip_B)
#     Rule('DIP_C', Cell(barcode = '3', dip = 'negC') >> Cell(barcode = '3', dip = 'negC') + 
#          Cell(barcode = '3', dip = 'negC'), k_dip_C)
#     Rule('DIP_D', Cell(barcode = '4', dip = 'negD') >> Cell(barcode = '4', dip = 'negD') + 
#          Cell(barcode = '4', dip = 'negD'), k_dip_D)
#     Rule('DIP_E', Cell(barcode = '5', dip = 'negE') >> Cell(barcode = '5', dip = 'negE') + 
#          Cell(barcode = '5', dip = 'negE'), k_dip_E)
#     Rule('DIP_F', Cell(barcode = '6', dip = 'posE') >> Cell(barcode = '6', dip = 'posE') + 
#          Cell(barcode = '6', dip = 'posE'), k_dip_F)
#     Rule('DIP_G', Cell(barcode = '7', dip = 'posD') >> Cell(barcode = '7', dip = 'posD') + 
#          Cell(barcode = '7', dip = 'posD'), k_dip_G)
#     Rule('DIP_H', Cell(barcode = '8', dip = 'posC') >> Cell(barcode = '8', dip = 'posC') + 
#          Cell(barcode = '8', dip = 'posC'), k_dip_H)
#     Rule('DIP_I', Cell(barcode = '9', dip = 'posB') >> Cell(barcode = '9', dip = 'posB') + 
#          Cell(barcode = '9', dip = 'posB'), k_dip_I)
#     Rule('DIP_J', Cell(barcode = '10', dip = 'posA') >> Cell(barcode = '10', dip = 'posA') + 
#          Cell(barcode = '10', dip = 'posA'), k_dip_J)
declare_functions()
print(model.rules)