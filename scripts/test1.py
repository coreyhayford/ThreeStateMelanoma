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


import random
import pandas as pd



# q = np.random.normal(size = 100)
# print(q)
# plt.hist(q)
# plt.show()
x = linspace(0,10,10)
a = np.random.multinomial(100, [1/50., 1/25., 2/25., 4/25., 1/5., 1/5., 4/25., 2/25., 1/25., 1/50.])
print(a)
plt.bar(x,a)
plt.show()

quit()
b = np.linspace(0,10,10, endpoint=False)
c = multivariate_normal.pdf(b, mean=2.5, cov=0.5)
print(c)
#plt.bar(x,c)



sampleDist = (0.02, 0.04, 0.08, 0.16, 0.20, 0.20, 0.16, 0.08, 0.04, 0.02)

randRolls = []
results = []

def roll(Dist, num):
    for i in range(num):
        randRoll = random.random()
        #print(randRoll)
        randRolls.append(randRoll) # in [0,1)
        sum = 0
        result = 1
        for val in Dist:
            sum += val
            #print("Mass: " + str(val))
            
            if randRoll < sum:
                #print(result)
                results.append(result)
                break
            result += 1
            #print("Sum " + str(sum))
            #print(result)

    print(randRolls)
    print(results)
    
    ps = pd.Series(results)
    global counts
    counts = ps.value_counts(sort = False)
    print(counts)
    
roll(sampleDist, 100)

counts.plot.bar()
plt.show()
quit()
# print(results)
# results = np.asarray(results)
# print(results)

# declare Monomers function with value

# [Parameter("Cell%d_0" %i) for i in range(11)]
# print(model.parameters)
# quit()
# [Parameter("Cell%d_0" %i) for i in range(10)]
# print(model.parameters)
print(results)
print(counts)
for i,j in enumerate(counts):
#     print("I: " + str(i))
#     print("J: " + str(j))
    Parameter("Cell%d_0" %i, j)
#     for k in range(10):
#         if k in counts.index:
#             Parameter("Cell%d_0" %i, j)
#         else:
#             Parameter("Cell%d_0" %i, 0)

print(model.parameters)
quit()
# for i,j in enumerate(results):
#     Parameter("Cell%d" %i, j)
    #[Parameter("Cell%d_0" %i, res) for i in range(11)]
    
# print(model.parameters)

quit()


[Monomer("M_%d" %i) for i in range(100)]

print model.monomers