"""
Model of Barcoded Cell Dynamics in the absence and presence of drug treatment.
Implemented By: Corey Hayford

"""

from pysb import *
from pysb.util import *

def declare_monomers():
    Monomer("Cell")
    alias_model_components()
def declare_parameters():
    Parameter("Cell_0", 1)
    Parameter("k_div", 0.05)
    Parameter("k_death", 0.025)
    alias_model_components()

def declare_ICs():
    Initial(Cell, Cell_0)
    alias_model_components()

def declare_observables():
    Observable("Obs_Cell", Cell)
    alias_model_components()

def declare_rules():
    Rule('div_Cell', Cell() >> Cell() + Cell(), k_div)
    Rule('death_Cell', Cell() >> None, k_death)
    alias_model_components()