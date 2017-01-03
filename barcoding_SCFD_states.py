from pysb import *
from pysb.util import *
from pysb.macros import *
from pysb.bng import *
from pysb.integrate import odesolve
import pylab as pl
from numpy import linspace
from sympy import sympify
import numpy as np


Model()

def declare_monomers():
    
    Monomer('Cell', ['barcode', 'dip'], {'barcode':['1', '2', '3', '4', '5', '6', '7', '8', '9', '10'], 
                                         'dip':['negA', 'negB', 'negC', 'negD', 'negE', 'posE', 
                                                'posD','posC', 'posB', 'posA']})

def declare_parameters():
    
### Initial Conditions

    Parameter("CellA_0", 10)
    Parameter("CellB_0", 10)
    Parameter("CellC_0", 10)
    Parameter("CellD_0", 10)
    Parameter("CellE_0", 10)
    Parameter("CellF_0", 10)
    Parameter("CellG_0", 10)
    Parameter("CellH_0", 10)
    Parameter("CellI_0", 10)
    Parameter("CellJ_0", 10)
    
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

def declare_initial_conditions():
    
    Initial(Cell(barcode = '1', dip = 'negA'), CellA_0)
    Initial(Cell(barcode = '2', dip = 'negB'), CellB_0)
    Initial(Cell(barcode = '3', dip = 'negC'), CellC_0)
    Initial(Cell(barcode = '4', dip = 'negD'), CellD_0)
    Initial(Cell(barcode = '5', dip = 'negE'), CellE_0)
    Initial(Cell(barcode = '6', dip = 'posE'), CellF_0)
    Initial(Cell(barcode = '7', dip = 'posD'), CellG_0)
    Initial(Cell(barcode = '8', dip = 'posC'), CellH_0)
    Initial(Cell(barcode = '9', dip = 'posB'), CellI_0)
    Initial(Cell(barcode = '10', dip = 'posA'), CellJ_0)

def declare_observables():
    
    Observable("Obs_A", Cell(barcode = '1', dip = 'negA'))
    Observable("Obs_B", Cell(barcode = '2', dip = 'negB'))
    Observable("Obs_C", Cell(barcode = '3', dip = 'negC'))
    Observable("Obs_D", Cell(barcode = '4', dip = 'negD'))
    Observable("Obs_E", Cell(barcode = '5', dip = 'negE'))
    Observable("Obs_F", Cell(barcode = '6', dip = 'posE'))
    Observable("Obs_G", Cell(barcode = '7', dip = 'posD'))
    Observable("Obs_H", Cell(barcode = '8', dip = 'posC'))
    Observable("Obs_I", Cell(barcode = '9', dip = 'posB'))
    Observable("Obs_J", Cell(barcode = '10', dip = 'posA'))
    Observable("Obs_All", Cell(barcode = '1', dip = 'negA') + Cell(barcode = '2', dip = 'negB') +
               Cell(barcode = '3', dip = 'negC') + Cell(barcode = '4', dip = 'negD') +
               Cell(barcode = '5', dip = 'negE') + Cell(barcode = '6', dip = 'posE') +
               Cell(barcode = '7', dip = 'posD') + Cell(barcode = '8', dip = 'posC') +
               Cell(barcode = '9', dip = 'posB') + Cell(barcode = '10', dip = 'posA'))
    
def declare_functions():
    
    Rule('DIP_A', Cell(barcode = '1', dip = 'negA') >> Cell(barcode = '1', dip = 'negA') + 
         Cell(barcode = '1', dip = 'negA'), k_dip_A)
    Rule('DIP_B', Cell(barcode = '2', dip = 'negB') >> Cell(barcode = '2', dip = 'negB') + 
         Cell(barcode = '2', dip = 'negB'), k_dip_B)
    Rule('DIP_C', Cell(barcode = '3', dip = 'negC') >> Cell(barcode = '3', dip = 'negC') + 
         Cell(barcode = '3', dip = 'negC'), k_dip_C)
    Rule('DIP_D', Cell(barcode = '4', dip = 'negD') >> Cell(barcode = '4', dip = 'negD') + 
         Cell(barcode = '4', dip = 'negD'), k_dip_D)
    Rule('DIP_E', Cell(barcode = '5', dip = 'negE') >> Cell(barcode = '5', dip = 'negE') + 
         Cell(barcode = '5', dip = 'negE'), k_dip_E)
    Rule('DIP_F', Cell(barcode = '6', dip = 'posE') >> Cell(barcode = '6', dip = 'posE') + 
         Cell(barcode = '6', dip = 'posE'), k_dip_F)
    Rule('DIP_G', Cell(barcode = '7', dip = 'posD') >> Cell(barcode = '7', dip = 'posD') + 
         Cell(barcode = '7', dip = 'posD'), k_dip_G)
    Rule('DIP_H', Cell(barcode = '8', dip = 'posC') >> Cell(barcode = '8', dip = 'posC') + 
         Cell(barcode = '8', dip = 'posC'), k_dip_H)
    Rule('DIP_I', Cell(barcode = '9', dip = 'posB') >> Cell(barcode = '9', dip = 'posB') + 
         Cell(barcode = '9', dip = 'posB'), k_dip_I)
    Rule('DIP_J', Cell(barcode = '10', dip = 'posA') >> Cell(barcode = '10', dip = 'posA') + 
         Cell(barcode = '10', dip = 'posA'), k_dip_J)