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

import preDrugDynamics

Model()

preDrugDynamics.declare_monomers()
preDrugDynamics.declare_observables()
preDrugDynamics.declare_newparameters()

def declare_parameters():
    Parameter("postDrugA_0", postDrugIC_A)
    Parameter("postDrugB_0", postDrugIC_B)
    Parameter("postDrugC_0", postDrugIC_C)
    Parameter("postDrugD_0", postDrugIC_D)
    Parameter("postDrugE_0", postDrugIC_E)
    Parameter("postDrugF_0", postDrugIC_F)
    Parameter("postDrugG_0", postDrugIC_G)
    Parameter("postDrugH_0", postDrugIC_H)
    Parameter("postDrugI_0", postDrugIC_I)
    Parameter("postDrugJ_0", postDrugIC_J)
    
    mean = 0
    var = 0.01
    pdf_dist = scipy.stats.norm.pdf(x_nvals, loc = mean, scale = np.sqrt(var)) #, endpoint, retstep, dtype))
    a = np.random.multinomial(post_drugIC_All, pdf_dist)
    print(a)

declare_parameters()