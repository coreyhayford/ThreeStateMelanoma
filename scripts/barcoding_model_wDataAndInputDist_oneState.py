from pysb import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sp
from pysb.bng import run_ssa
from pysb.simulator.bng import BngSimulator
from pysb.bng import generate_equations
import pandas as pd


def run_barcoding_model(n_exp, n_barcodes, n_cell_types):
    all_post_drug_counts = []

    # Modified variable names for cFP distribution
    # slopes = pd.read_csv("~/Documents/QuarantaLab/SKMEL5_H2BGFP_cFPslopes_rep1.csv")
    slopes = pd.read_csv("/data/lola/hayforc/ParamScan/barcoding/SKMEL5_H2BGFP_cFPslopes_rep1.csv")
    slopes = np.array(slopes["x"])
    low_dips = min(slopes)
    high_dips = max(slopes)
    dips = np.linspace(low_dips, high_dips, n_cell_types)
    dip_mean = np.mean(slopes)
    dip_var = np.std(slopes)
    normal = sp.norm.pdf(dips, dip_mean, dip_var)
    normal_hist = normal * (dips[1] - dips[0])
    normal_hist /= normal_hist.sum()

    Model()
    Monomer("Cell", ["dip"], {"dip": ["x" + str(i) for i in range(5)]})
    # [Monomer("Cell", ['dip'], {'dip': ["%d" %i for i in range(n_cell_types)]})]
    print(model.monomers)

    Parameter('cellInit_0', 1)
    [Parameter('cellInit_%d' % i) for i in range(1, n_cell_types)]
    print(model.parameters)

    Initial(Cell(dip="x0"), cellInit_0)  # could not be a string - parameter only - didn't know how to set up
    [Initial(Cell(dip="x" + str(i)), model.parameters["cellInit_%d" % i]) for i in range(1, n_cell_types)]
    print(model.initial_conditions)

    [Observable("Obs_Cell%d" % i, Cell(dip="x" + str(i))) for i in range(n_cell_types)]
    Observable("Obs_All", Cell())
    print(model.observables)

    print(dips)
    k_div = 0.03
    [Parameter("k_div_%d" % i, k_div) for i in range(len(dips))]
    k_death = -dips + k_div
    [Parameter("k_death_%d" % i, k) for i, k in enumerate(k_death)]
    print(model.parameters)

    [Rule("Cell%d_Div" % i, Cell(dip="x" + str(i)) >> Cell(dip="x" + str(i)) + Cell(dip="x" + str(i)),
          model.parameters["k_div_%d" % i]) for i in range(len(dips))]
    [Rule("Cell%d_Death" % i, Cell(dip="x" + str(i)) >> None,
          model.parameters["k_death_%d" % i]) for i in range(len(k_death))]

    print(model.rules)

    for exp in range(n_exp):
        t1 = np.linspace(0, 336, 337)  # 14 days
        sim = BngSimulator(model, tspan=t1, verbose=False)  # , seed = 1095205711)
        # Why do the n_cells = 0 start off with 1?
        x1 = sim.run(n_runs=n_barcodes, param_values={"k_death_0": 0.005}, verbose=False)

        ICs = pd.read_csv("/data/lola/hayforc/ParamScan/barcoding/barcodeAbundance_forModel.csv")
        ICs = np.array(ICs["x"])
        ICs_top = ICs[:500]
        ICs_top = list(ICs_top.astype(np.float))
        cell_tot = ICs_top

        post_drug_counts = []

        for ind, cell in enumerate(cell_tot):
            print("%d:%d" % (exp, ind))
            # rounding and taking integer is safe way to handle number.0
            b = np.random.choice(n_cell_types, 1, p=normal_hist)
            c = np.zeros((n_cell_types,));
            c[b] = 1
            e = int(round(cell))
            cell_pop = np.random.multinomial(e, c)
            # cell_pop = np.random.multinomial(int(round(cell)), normal_hist)  # *100 for true normal
            t2 = np.linspace(0, 192, 193)  # 8 days in drug
            # sim = BngSimulator(model, tspan=t2, verbose=False)
            x2 = sim.run(tspan=t2, param_values={"k_death_0": model.parameters['k_death_0'].value},
                         n_sim=1, initials={model.species[i]: pop for i, pop in enumerate(cell_pop)}, verbose=False)

            post_drug_counts.append(x2.all["Obs_All"][-1])
        all_post_drug_counts.append(post_drug_counts)

    np.save("barcoding_500barcodes_Simdata_%dexp%dstates.npy" % (e, s), all_post_drug_counts)


states = [2, 3, 4, 5, 6]
experiments = [100]
barcodes = [500]

for s in states:
    for e in experiments:
        for b in barcodes:
            run_barcoding_model(e, b, s)
