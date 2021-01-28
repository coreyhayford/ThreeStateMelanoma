from pysb import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sp
from pysb.simulator.bng import BngSimulator
# from pysb.bng import generate_equations
# import pandas as pd
# import seaborn as sns
import random

# sns.set(font_scale = 1.25)
# sns.set_style("whitegrid")

# Setting Key Model Parameters
# n_exp = 5
# n_barcodes = 5
# n_cell_types = 5

def run_barcoding_model(n_exp, n_barcodes, n_cell_types, percent):

    all_post_drug_counts = []

    ### SKMEL5 SC01 DIP Rate Distribution
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
    Monomer("Cell_resistant")
    print(model.monomers)

    # Monomer("Cell", ['dip'], {"dip": ["0","1"]})


    Parameter('cellInit_0', 1)
    [Parameter('cellInit_%d' %i)for i in range(1,n_cell_types)]
    # Resistant initial parameter will be overwritten later
    Parameter('Cell_resistant_0', 0)  # list(choice).count(0))
    # Parameter('Cell_1init')
    print(model.parameters)

    Initial(Cell(dip = "0"), cellInit_0) # could not be a string - parameter only - didn't know how to set up
    [Initial(Cell(dip = "%d" %i), model.parameters["cellInit_%d" %i]) for i in range(1,n_cell_types)]
    Initial(Cell_resistant, Cell_resistant_0)
    # for ic in model.initial_conditions:
    #     print ic
    print(model.initial_conditions)

    [Observable("Obs_Cell%d" %i, Cell(dip = "%d" %i)) for i in range(n_cell_types)]
    Observable("Obs_All", Cell())
    Observable("Obs_Cell_resistant", Cell_resistant())
    print(model.observables)

    print(dips)
    k_div = 0.04 * np.log(2)
    [Parameter("k_div_%d" % i, k_div) for i in range(len(dips))]
    k_death = k_div - (dips*np.log(2))
    [Parameter("k_death_%d" % i, k) for i,k in enumerate(k_death)]
    Parameter('k_div_resistant', 0.04 * np.log(2))
    Parameter('k_death_resistant', 0.005 * np.log(2))

    [Rule("Cell%d_Div" %i, Cell(dip = "%d" %i) >> Cell(dip = "%d" %i) + Cell(dip = "%d" %i),
          model.parameters["k_div_%d" %i]) for i in range(len(dips))]
    [Rule("Cell%d_Death" %i, Cell(dip = "%d" %i) >> None,
          model.parameters["k_death_%d" %i]) for i in range(len(k_death))]
    Rule('Cell_resistant_Div', Cell_resistant() >> Cell_resistant() + Cell_resistant(), k_div_resistant)
    Rule('Cell_resistant_Death', Cell_resistant() >> None, k_death_resistant)


    print(model.rules)

    test = np.arange(n_barcodes)
    print(test)
    print(list(test))
    test1 = random.sample(list(test), int(round(percent * len(test))))
    print(test1)

    for exp in range(n_exp):

        # Pre-drug run - all have the same DIP rate +/- stochasticity
        t1 = np.linspace(0,336,337) # 14 days
        sim = BngSimulator(model, tspan=t1, verbose=False)#, seed = 1095205711)
        # Why do the n_cells = 0 start off with 1? SOLVED
        # Use n_runs = n_barcodes
        x1 = sim.run(n_runs = n_barcodes, param_values={"k_death_0": 0.005 * np.log(2)},
                     verbose = False) # returns np.array with species and obs
        # all_cells = np.array(x1.observables)
        # cell_tot = all_cells["Obs_All"].T[-1]
        # print(np.array(x1.all).shape)
        # quit()
        # print(x1.all)
        # print(type(x1.all))
        cell_tot = [x[-1]["Obs_Cell0"] for x in x1.all]
        print(cell_tot)
        quit()
        # print(len(cell_tot))
        # print(int(round(percent*len(cell_tot))))
        # print(list(enumerate(cell_tot)))

        # cell_tot_resistant_list = random.sample(list(enumerate(cell_tot)), int(round(percent*len(cell_tot))))

        # cell_tot_resistant = np.random.choice(list(enumerate(cell_tot)), int(round(percent*len(cell_tot))))


        # resistant_indecies = []
        # for index,elem in cell_tot_resistant_list:
        #     resistant_indecies.append(index)

        # sensitive_indecies = [i for i,x in enumerate(cell_tot) if i not in resistant_indecies]
        # print(resistant_indecies)
        # print(sensitive_indecies)
        #
        # # print(cell_tot_resistant[[0]])
        # cell_tot_resistant = [x for i,x in cell_tot_resistant_list]
        # cell_tot_sensitive = [x for i,x in enumerate(cell_tot) if i not in resistant_indecies]
        #
        # print(cell_tot_resistant)
        # print(cell_tot_sensitive)
        #
        # print(model.species)
        # quit()

        print("Sim 1 Finished.")
        post_drug_counts = []

        for ind, cell in enumerate(cell_tot):
            print("%d:%d" % (exp,ind))
            if cell == 0:
                post_drug_counts.append(0.)
            elif ind in test1:
                t2 = np.linspace(0, 168, 169)  # 7 days in drug
                x2 = sim.run(tspan=t2, param_values={"k_death_0": 0.005 * np.log(2)},
                             n_sim=1, initials={model.species[0]: cell}, verbose=False)

                post_drug_counts.append(x2.all["Obs_All"][-1])

                # plt.figure()
                # plt.plot(x1.tout[ind], np.log2(x1.all[ind]["Obs_Cell0"]), '0.5', lw=2)
                # plt.plot(t1[-1] + x2.tout[0], np.log2(x2.observables["Obs_All"]), '0.5', lw=2)
                # plt.axvline(x=336, linewidth=4, color='red')
                # plt.axvline(x=504, linewidth=4, color='blue')
                # plt.xlabel("Time (hours)")
                # plt.ylabel("log2(cell count)")
                # plt.xlim([0, 550])
                # plt.title("Barcoding Model - 10 barcodes, 1 experiments, 10 states")
            else:
                # rounding and taking integer is safe way to handle number.0
                cell_pop = np.random.multinomial(int(round(cell)), normal_hist) # *100 for true normal
                # print(cell_pop)

                # # Normal Histogram Multinomial Selection of DIP Rates
                # plt.figure("DIP Rate Distribution Barcode #%d" % ind)
                # plt.bar(dips-0.5*(dips[1]-dips[0]), 1.*cell_pop/int(round(cell)), dips[1]-dips[0],
                #         label = "n_cells = %d" %int(cell)) # *100 for true normal
                # plt.plot(dips, normal_hist, lw = 2, color = "r")
                # plt.legend(loc = 0)
                t2 = np.linspace(0,168,169) # 7 days in drug

                x2 = sim.run(tspan = t2, param_values={"k_death_0": model.parameters['k_death_0'].value},
                             n_sim=1, initials={model.species[i]:pop for i,pop in enumerate(cell_pop)}, verbose= False)
                # ,"cellInit_0": 0}, - FIXED THIS PROBLEM

                post_drug_counts.append(x2.all["Obs_All"][-1])
                print(x1.all[ind]["Obs_Cell0"])
                # plt.figure()
                # plt.plot(x1.tout[ind], np.log2(x1.all[ind]["Obs_Cell0"]), '0.5', lw = 2)
                # a = np.array(x1.observables)["Obs_Cell0"]
                # print(np.array([tr[:] for tr in a]))
                # quit()
                # plt.plot(x1.tout, np.log2(np.array([tr[:] for tr in a])), '0.5')
                # plt.plot(t1[-1] + x2.tout[0], np.log2(x2.observables["Obs_All"]), '0.5', lw = 2)
                # plt.axvline(x = 336, linewidth = 4, color = 'red')
                # plt.axvline(x = 504, linewidth = 4, color = 'blue')
                # plt.xlabel("Time (hours)")
                # plt.ylabel("log2(cell count)")
                # plt.xlim([0,550])
                # plt.title("Barcoding Model - 10 barcodes, 1 experiments, 10 states")
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
    np.save("barcoding_data_%dbar%dexp%dstates%.6fspikein.npy" % (b, e, s, si), all_post_drug_counts)
    # np.save("barcoding_data_10bar5exp10states_0.2spikein.npy", all_post_drug_counts)

    # plt.show()

# run_barcoding_model(n_exp=5, n_barcodes=10, n_cell_types=10, percent=0.2)
# quit()

states = [15]
experiments = [1,5,10,15,20] #range(1,11)
barcodes = [100,500,1000]
spike_ins = [0,0.000001,0.00001,0.0001,0.001,0.01,0.05,0.10,0.20,0.25,0.50,1]

for s in states:
    for e in experiments:
        for b in barcodes:
            for si in spike_ins:
                run_barcoding_model(e,b,s,si)