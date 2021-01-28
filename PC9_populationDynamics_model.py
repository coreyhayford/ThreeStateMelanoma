from pysb import *
from pysb.util import *
import numpy as np

num_cells = 5000

Model()


def declare_monomers():
    Monomer("PC9_VU")
    Monomer("PC9_BR1")


def declare_parameters():
    ### Initial Conditions ###

    Parameter("PC9_VU_0", num_cells * 1)
    Parameter("PC9_BR1_0", num_cells * 0)

    ## Growth Parameters ###
    Parameter("k_dip_PC9_VU", 0.0025 * np.log(2))
    Parameter("k_dip_PC9_BR1",  0.023 * np.log(2))

    #     Parameter("k_div_PC9_VU", 0)
    #     Parameter("k_div_PC9_BR1", 0)
    #
    #     Parameter("k_death_PC9_VU", 0)
    #     Parameter("k_death_PC9_BR1",  0)

    alias_model_components()


def declare_initial_conditions():
    Initial(PC9_VU, PC9_VU_0)
    Initial(PC9_BR1, PC9_BR1_0)


def declare_observables():
    Observable("Obs_PC9_VU", PC9_VU)
    Observable("Obs_PC9_BR1", PC9_BR1)
    Observable("Obs_All", PC9_VU() + PC9_BR1())


def declare_functions():

    Rule('DIP_PC9_VU', PC9_VU() >> PC9_VU() + PC9_VU(), k_dip_PC9_VU)
    Rule('DIP_PC9_BR1', PC9_BR1() >> PC9_BR1() + PC9_BR1(), k_dip_PC9_BR1)

    # Rule('Divide_PC9_VU', PC9_VU() >> PC9_VU() + PC9_VU(), k_div_PC9_VU)
    # Rule('Death_PC9_VU', PC9_VU() >> None, k_death_PC9_VU)
    # Rule('Divide_PC9_BR1', PC9_BR1() >> PC9_BR1() + PC9_BR1(), k_div_PC9_BR1)
    # Rule('Death_PC9_BR1', PC9_BR1() >> None, k_death_PC9_BR1)

declare_monomers()
declare_parameters()
declare_initial_conditions()
declare_observables()
declare_functions()
#
# print(model.monomers)
# print(model.parameters)
# print(model.initial_conditions)
# print(model.observables)
# print(model.rules)