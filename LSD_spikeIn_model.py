from pysb import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sp
from pysb.bng import run_ssa
from pysb.simulator.bng import BngSimulator
from pysb.bng import generate_equations
import pandas as pd
import seaborn as sns

choice = np.random.choice(2, np.random.poisson(num_cells), p=[percent, 1 - percent])


def declare_monomers():
    [Monomer("Cell", ['dip'], {'dip': ["%d" % i for i in range(n_cell_types)]})]
    Monomer("Cell_resistant")

def declare_parameters():
    [Parameter('cellInit_%d' %i)for i in range(1,n_cell_types)]
    Parameter('Cell_resistant_0', 0)