from pysb import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sp
from pysb.bng import run_ssa
from pysb.simulator.bng import BngSimulator
from pysb.bng import generate_equations
import pandas as pd
import seaborn as sns

sns.set(font_scale = 1.25)
sns.set_style("whitegrid")

# Setting Key Model Parameters
# n_exp = 5
# n_barcodes = 5
# n_cell_types = 5

def run_barcoding_model(n_exp, n_barcodes, n_cell_types):

    all_post_drug_counts = []
    low_dips = -0.062
    high_dips = -0.002
    dips = np.linspace(low_dips, high_dips, n_cell_types)
    dip_mean = -0.033
    dip_var = 0.01

    # Discretize normal distribution of dip rates - used in post drug simulation
    normal = sp.norm.pdf(dips, dip_mean, dip_var)
    print(normal)
    sum = 0
    for i in range(1, n_cell_types):
        sum += normal[i] * (dips[i] - dips[i - 1])
    print(sum)

    normal_hist = normal * (dips[1] - dips[0])
    print(normal_hist)
    print(dips)


    Model()
    [Monomer("Cell", ['dip'], {'dip': ["%d" %i for i in range(n_cell_types)]})]
    print(model.monomers)

    # Monomer("Cell", ['dip'], {"dip": ["0","1"]})


    Parameter('cellInit_0', 1)
    [Parameter('cellInit_%d' %i)for i in range(1,n_cell_types)]
    # Parameter('Cell_1init')
    print(model.parameters)

    Initial(Cell(dip = "0"), cellInit_0) # could not be a string - parameter only - didn't know how to set up
    [Initial(Cell(dip = "%d" %i), model.parameters["cellInit_%d" %i]) for i in range(1,n_cell_types)]
    # for ic in model.initial_conditions:
    #     print ic
    print(model.initial_conditions)

    [Observable("Obs_Cell%d" %i, Cell(dip = "%d" %i)) for i in range(n_cell_types)]
    Observable("Obs_All", Cell())
    print(model.observables)

    print(dips)
    k_div = 0.03
    [Parameter("k_div_%d" % i, k_div) for i in range(len(dips))]
    k_death = -dips+k_div
    [Parameter("k_death_%d" % i, k) for i,k in enumerate(k_death)]
    print(model.parameters)

    [Rule("Cell%d_Div" %i, Cell(dip = "%d" %i) >> Cell(dip = "%d" %i) + Cell(dip = "%d" %i),
          model.parameters["k_div_%d" %i]) for i in range(len(dips))]
    [Rule("Cell%d_Death" %i, Cell(dip = "%d" %i) >> None,
          model.parameters["k_death_%d" %i]) for i in range(len(k_death))]

    print(model.rules)

    for exp in range(n_exp):

        t1 = np.linspace(0,336,337) # 14 days
        sim = BngSimulator(model, tspan=t1, verbose=False)#, seed = 1095205711)
        # Why do the n_cells = 0 start off with 1?
        x1 = sim.run(n_runs = n_barcodes, param_values={"k_death_0": 0.005}, verbose = False) # returns np.array with species and obs
        # all_cells = np.array(x1.observables)
        # cell_tot = all_cells["Obs_All"].T[-1]
        # print(np.array(x1.all).shape)
        # quit()
        # print(x1.all)
        # print(type(x1.all))
        cell_tot = [x[-1]["Obs_Cell0"] for x in x1.all]
        print(cell_tot)

        print("Sim 1 Finished.")

        post_drug_counts = []

        for ind, cell in enumerate(cell_tot):
            print("%d:%d" % (exp,ind))
            if cell == 0:
                post_drug_counts.append(0.)
            else:
                # rounding and taking integer is safe way to handle number.0
                cell_pop = np.random.multinomial(int(round(cell)), normal_hist) # *100 for true normal
                # print(cell_pop)

                # Normal Histogram Multinomial Selection of DIP Rates
                plt.figure("DIP Rate Distribution Barcode #%d" % ind)
                plt.bar(dips-0.5*(dips[1]-dips[0]), 1.*cell_pop/int(round(cell)), dips[1]-dips[0], label = "n_cells = %d" %int(cell)) # *100 for true normal
                plt.plot(dips, normal_hist, lw = 2, color = "r")
                plt.legend(loc = 0)
                t2 = np.linspace(0,168,169) # 7 days in drug

                x2 = sim.run(tspan = t2, param_values={"k_death_0": model.parameters['k_death_0'].value},
                             n_sim=1, initials={model.species[i]:pop for i,pop in enumerate(cell_pop)}, verbose= False)
                # ,"cellInit_0": 0}, - FIXED THIS PROBLEM

                post_drug_counts.append(x2.all["Obs_All"][-1])
                print(x1.all[ind]["Obs_Cell0"])
                plt.figure()
                plt.plot(x1.tout[ind], np.log2(x1.all[ind]["Obs_Cell0"]), '0.5', lw = 2)
                a = np.array(x1.observables)["Obs_Cell0"]
                print(np.array([tr[:] for tr in a]))
                # quit()
                # plt.plot(x1.tout, np.array([tr[:] for tr in a]), '0.5')
                plt.plot(t1[-1] + x2.tout[0], np.log2(x2.observables["Obs_All"]), '0.5', lw = 2)
                plt.axvline(x = 336, linewidth = 4, color = 'red')
                plt.axvline(x = 504, linewidth = 4, color = 'blue')
                plt.xlabel("Time (hours)")
                plt.ylabel("log2(cell count)")
                plt.xlim([0,550])
                plt.title("Barcoding Model - 100 barcodes, 10 experiments, 15 states")
                # plt.figure("Experiment #%d Barcode #%d" % (exp,ind))
                # plt.plot(t1, x1.all[int(cell)]["Obs_Cell0"],'0.5', lw=4, alpha=0.25)
                # for i,obs in enumerate(model.observables):
                    # print(x1.all[int(cell)])
                    # quit()
                    # plt.plot(t1, x1.all[obs.name], lw = 2, label = obs.name, color = colors[i])
                    # plt.plot(t1, x1.all[obs.name], '0.5', lw=4, alpha=0.25)
                    # plt.plot(t1[-1]+t2, x2.all[obs.name], '0.5', lw = 4, alpha = 0.25)
                    # plt.plot(t1, x1.all[obs.name], lw = 2, label = obs.name, color = colors[i])
                    # plt.plot(t1[-1]+t2, x2.all[obs.name], lw = 2, label = obs.name)
                # plt.annotate(s = "n_cells = %d" % cell, xy = (0.1,0.1), xycoords = "axes fraction", fontsize = 24)
                # plt.text(-0.04, 2800, r'$\mu=%g$' % (np.mean(distr_all)))
                # plt.legend(loc = 0)
        print(post_drug_counts)
        all_post_drug_counts.append(post_drug_counts)

    print(all_post_drug_counts)
    # np.save("barcoding_data_%dbar%dexp%dstates.npy" % (b, e, s), all_post_drug_counts)
    np.save("barcoding_data_100bar10exp15states_forpres.npy", all_post_drug_counts)

    plt.show()

run_barcoding_model(n_exp=1, n_barcodes=10, n_cell_types=10)
quit()
states = [15]
experiments = range(1,11)
barcodes = [5]

for s in states:
    for e in experiments:
        for b in barcodes:
            run_barcoding_model(e,b,s)
