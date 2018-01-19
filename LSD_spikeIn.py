from pysb import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sp
from pysb.bng import run_ssa
from pysb.simulator.bng import BngSimulator
from pysb.bng import generate_equations
import pandas as pd
import seaborn as sns
import math

# a = [2]
# print(len(a))
# print(np.sum(a) != 0)
# print(len(a) != 1)
# quit()

# rs = np.random.RandomState(0)
# num_elements = 1
# x = rs.normal(10, 1, num_elements)
# sns.distplot(np.r_[x, x], kde=False, norm_hist=True)
# plt.show()
# quit()

sns.set(font_scale = 1.25)
sns.set_style("white")


# Setting Key Model Parameters
# n_exp = 5
# n_barcodes = 5
# n_cell_types = 5

def run_LSDspikeIn_model(n_wells, n_cell_types, percent, n_cells):

    low_dips = -0.04
    high_dips = 0.02
    dips = np.linspace(low_dips, high_dips, n_cell_types)
    dip_mean = -0.01
    dip_var = 0.01
    # print("DIP RATES")
    # print(dips)
    # print("Death rates")
    # print(-dips+(0.03*np.log(2)))


    # Discretize normal distribution of dip rates - used in post drug simulation
    normal = sp.norm.pdf(dips, dip_mean, dip_var)
    sum = 0
    for i in range(1, n_cell_types):
        sum += normal[i] * (dips[i] - dips[i - 1])
    normal_hist = normal * (dips[1] - dips[0])

    Model()
    [Monomer("Cell", ['dip'], {'dip': ["%d" %i for i in range(n_cell_types)]})]
    Monomer("Cell_resistant")
    # print(model.monomers)

    Parameter('cellInit_0', 0) # list(choice).count(1))
    [Parameter('cellInit_%d' %i)for i in range(1,n_cell_types)]
    Parameter('Cell_resistant_0', 0) # list(choice).count(0))
    # print(model.parameters)

    Initial(Cell(dip = "0"), cellInit_0) # could not be a string - parameter only - didn't know how to set up
    [Initial(Cell(dip = "%d" %i), model.parameters["cellInit_%d" %i]) for i in range(1,n_cell_types)]
    Initial(Cell_resistant, Cell_resistant_0)
    # print(model.initial_conditions)

    [Observable("Obs_Cell%d" %i, Cell(dip = "%d" %i)) for i in range(n_cell_types)]
    Observable("Obs_All", Cell())
    Observable("Obs_Cell_resistant", Cell_resistant())
    # print(model.observables)

    # print(dips)
    k_div = 0.04 * np.log(2)
    [Parameter("k_div_%d" % i, k_div) for i in range(len(dips))]
    k_death = k_div - (dips*np.log(2))
    [Parameter("k_death_%d" % i, k) for i,k in enumerate(k_death)]
    Parameter('k_div_resistant', 0.04 * np.log(2))
    Parameter('k_death_resistant', 0.005 * np.log(2))

    # print(model.parameters)

    [Rule("Cell%d_Div" %i, Cell(dip = "%d" %i) >> Cell(dip = "%d" %i) + Cell(dip = "%d" %i),
          model.parameters["k_div_%d" %i]) for i in range(len(dips))]
    [Rule("Cell%d_Death" %i, Cell(dip = "%d" %i) >> None,
          model.parameters["k_death_%d" %i]) for i in range(len(k_death))]
    Rule('Cell_resistant_Div', Cell_resistant() >> Cell_resistant() + Cell_resistant(), k_div_resistant)
    Rule('Cell_resistant_Death', Cell_resistant() >> None, k_death_resistant)

    # print(model.rules)
    # quit()

    distr_normal = []
    distr_resistant = []
    count = 0

    # preDrug_cellNormal = []
    # preDrug_cellResistant = []
    for exp in range(n_wells):

        #cFP
        # num_cells = 1
        #LSD
        num_cells = np.random.poisson(n_cells)
        # the 2 provides two choices (0 and 1), which correspond to two populations
        choice = np.random.choice(2, num_cells, p=[percent, 1-percent])
        # print(choice)

        if num_cells != 0:

            normal_init = list(choice).count(1)
            resistant_init = list(choice).count(0)

            cellInit_0.value = normal_init
            Cell_resistant_0.value = resistant_init

            # print(cellInit_0)
            # print(Cell_resistant_0)
            t1 = np.linspace(0,168,169) # 7 days
            sim = BngSimulator(model, tspan=t1, verbose=False)#, seed = 1095205711)
            # Why do the n_cells = 0 start off with 1?
            x1 = sim.run(n_runs = 1, param_values={"k_death_0": 0.005 * np.log(2)},
                         verbose = False) # returns np.array with species and obs
            # all_cells = np.array(x1.observables)
            # cell_tot = all_cells["Obs_All"].T[-1]
            # print(np.array(x1.all).shape)
            # quit()
            # print(x1.all)
            # print(type(x1.all))

            # The total well subpopulation count before drug addition
            cell_tot = x1.observables["Obs_Cell0"][-1]
            cell_tot_resistant = x1.observables["Obs_Cell_resistant"][-1]

            # preDrug_cellNormal.append(cell_tot)
            # preDrug_cellResistant.append(cell_tot_resistant)
            # For when n_wells > 1
            # cell_tot = [x[-1]["Obs_Cell0"] for x in x1.all]
            # cell_tot_resistant = [x[-1]["Obs_Cell_resistant"] for x in x1.all]
            # print(cell_tot)
            # print(cell_tot_resistant)

            # print("Sim 1 Finished.")

            cell_pop = np.random.multinomial(int(round(cell_tot)), normal_hist)  # *100 for true normal
            # print(cell_pop)
            t2 = np.linspace(0, 144, 145)  # 6 days in drug


            x2 = sim.run(tspan=t2, param_values={"k_death_0": model.parameters['k_death_0'].value,
                                                 "Cell_resistant_0": cell_tot_resistant},
                         n_sim=1, initials={model.species[i]: pop for i, pop in enumerate(cell_pop)},
                         verbose=False)
            # print(cell_tot_resistant)

            if np.sum(cell_pop) != 0:
                slope, intercept, r_value, p_value, std_err = sp.stats.linregress(t2, np.log2(
                    x2.observables["Obs_All"] / x2.observables["Obs_All"][0]))
                if math.isnan(slope) == False:
                    distr_normal = distr_normal + [slope]
                # print "Normal", slope

            if cell_tot_resistant != 0:
                slope, intercept, r_value, p_value, std_err = sp.stats.linregress(t2, np.log2(
                    x2.observables["Obs_Cell_resistant"] / x2.observables["Obs_Cell_resistant"][0]))
                if math.isnan(slope) == False:
                    distr_resistant = distr_resistant + [slope]
                # print "Resistant", slope

            # plt.figure()
            # plt.plot(t1, np.log2(x1.observables["Obs_Cell0"]), 'g', lw=2)
            # plt.plot(t1, np.log2(x1.observables["Obs_Cell_resistant"]), 'r', lw=2)
            # plt.plot(t1[-1] + x2.tout[0], np.log2(x2.observables["Obs_All"]), 'g', lw = 2)
            # plt.plot(t1[-1] + x2.tout[0], np.log2(x2.observables["Obs_Cell_resistant"]), 'r', lw = 2)
            # plt.axvline(x=168, linewidth=4, color='k')
            # plt.xlabel("Time (hours)", weight = "bold")
            # plt.ylabel("log2(cell count)", weight = "bold")
            # plt.title("LSD %d Resistant Spike-In: Well %d" % (percent*100, exp), weight = "bold")
        count = count + 1
        print("Well")
        print(count)
        print(distr_normal)
        print(distr_resistant)

    plt.figure()
    bins = np.linspace(-0.04, 0.06, 101)
    # bins = np.histogram(np.hstack((distr_resistant, distr_normal)), bins=100)[1]
    if len(distr_normal) >= 2:
        sns.distplot(distr_normal, bins, kde=False, color="g", kde_kws={"shade": True}, label = "Parental")
    # elif len(distr_normal) == 1:
    #     rs = np.random.RandomState(0)
    #     num_elements = len(distr_normal)
    #     x = rs.normal(10, 1, num_elements)
    #     sns.distplot(np.r_[x,x], bins, kde=False, color="g", kde_kws={"shade": True}, label="Parental")
    if len(distr_resistant) >= 2:
        sns.distplot(distr_resistant, bins, kde=False, color="r", kde_kws={"shade": True}, label = "Resistant")
    else:
        pass

    # elif len(distr_resistant) == 1:
    #     rs1 = np.random.RandomState(0)
    #     num_elements1 = len(distr_resistant)
    #     x1 = rs1.normal(10, 1, num_elements1)
    #     sns.distplot(np.r_[x1,x1], bins, kde=False, color="r", kde_kws={"shade": True},label="Resistant")
    plt.axvline(x=(k_div-((k_death_4.value+k_death_5.value)/2))/np.log(2), color = 'g', ls = '--', lw = 3,
                label = "Parental DIP Mean")
    plt.axvline(x=(k_div_resistant.value-k_death_resistant.value)/np.log(2), color = 'r', ls = '--', lw = 3,
                label = "Resistant DIP Mean")
    plt.axvline(x=0, color = 'k', ls = '--', lw = 3)

    # bins = np.histogram(np.hstack((distr_resistant, distr_normal)), bins=100)[1]
    plt.xlabel("DIP Rate", weight="bold")
    plt.ylabel("Frequency", weight = "bold")
    plt.xlim(-0.04, 0.06)
    plt.ylim(0,400)
    # plt.ylim(0,100*(n_wells/384))
    plt.title("LSD DIP Distribution %.4f%% Resistant Spike-In (%d Cells, %d Wells) " % (percent*100, n_cells, n_wells),
              weight = "bold")
    plt.legend(loc="upper left")
    np.save("LSDspikein_normalDIPs_%.6fresistant_%dcells_%dwells.npy" % (percent, n_cells, n_wells), distr_normal)
    np.save("LSDspikein_resistantDIPs_%.6fresistant_%dcells_%dwells.npy" % (percent, n_cells, n_wells), distr_resistant)
    plt.savefig("LSDspikein_DIPdistr_%.6fresistant_%dcells_%dwells.pdf" % (percent, n_cells, n_wells))

# per = 0.10
# wells = 384
# cells = 1
# run_LSDspikeIn_model(n_wells=wells, n_cell_types=10, percent=per, n_cells=cells)
# plt.show()
# quit()

# 80,065,152 pure simulations
# 145,152 simulations if only 96,384,1536 plates
# 100 replicates for each unique simulation to get CIs
percentages = [0,0.000001,0.00001,0.0001,0.001,0.01,0.05,0.10,0.20,0.25,0.50,1] #[0,0.000001,0.00001,0.0001,
cell_num_list = [5,10,25,50,100]
n_wells_list = [1536] #[96,384,1536] #,10000,100000,1000000]

exp_count = 1
for per in percentages:
    for cells in cell_num_list:
        for wells in n_wells_list:
            print("Experiment")
            print(exp_count)
            exp_count = exp_count + 1
            run_LSDspikeIn_model(n_wells=wells, n_cell_types=10, percent=per, n_cells=cells)

#
# plt.show()
quit()
            # plt.figure()
            # for i in range(n_wells):
            #     plt.plot(t1, np.log2(x1.observables["Obs_Cell0"]), 'g', lw=2)
            #     plt.plot(t1, np.log2(x1.observables["Obs_Cell_resistant"]), 'r', lw=2)
            #     # plt.axvline(x = 336, linewidth = 4, color = 'red')
            # # plt.axvline(x = 504, linewidth = 4, color = 'blue')
            # plt.xlabel("Time (hours)")
            # plt.ylabel("log2(cell count)")
            # post_drug_counts = []

    # print(preDrug_cellNormal)
    # print(preDrug_cellResistant)

    # preDrug_cellNormal_noZeroes = filter(lambda a: a != 0.0, preDrug_cellNormal)
    # preDrug_cellResistant_noZeroes = filter(lambda a: a != 0.0, preDrug_cellResistant)
    #
    # print(preDrug_cellNormal_noZeroes)
    # print(preDrug_cellResistant_noZeroes)

### Do a different for loop, or combine loops with resistant

    # for ind, cell in enumerate(preDrug_cellNormal):
    #     print("%d:%d" % (exp,ind))
    #     if cell != 0:
    #
    #         # rounding and taking integer is safe way to handle number
    #         cell_pop = np.random.multinomial(int(round(cell)), normal_hist)  # *100 for true normal
    #         # print(cell_pop)
    #
    #         # Normal Histogram Multinomial Selection of DIP Rates
    #         plt.figure("DIP Rate Distribution Normal #%d" % ind)
    #         plt.bar(dips-0.5*(dips[1]-dips[0]), 1.*cell_pop/int(round(cell)), dips[1]-dips[0],
    #                 label = "n_cells = %d" %int(cell)) # *100 for true normal
    #         plt.plot(dips, normal_hist, lw = 2, color = "g")
    #         plt.legend(loc = 0)
    #
    #
    #         t2 = np.linspace(0, 168, 169)  # 7 days in drug
    #
    #         x2 = sim.run(tspan=t2, param_values={"k_death_0": model.parameters['k_death_0'].value},
    #                      n_sim=1, initials={model.species[i]: pop for i, pop in enumerate(cell_pop)},
    #                      verbose=False)
    #
    #         plt.figure()
    #         plt.plot(t1, np.log2(x1.observables["Obs_Cell0"]), 'g', lw=2)
    #         plt.plot(t1, np.log2(x1.observables["Obs_Cell_resistant"]), 'r', lw=2)
    #
    #
    #         # plt.plot(x1.tout[ind], np.log2(x1.all[ind]["Obs_Cell0"]), '0.5', lw = 2)
    #         # a = np.array(x1.observables)["Obs_Cell0"]
    #         # print(np.array([tr[:] for tr in a]))
    #         # # plt.plot(x1.tout, np.array([tr[:] for tr in a]), '0.5')
    #         #
    #         # plt.plot(t1[-1] + x2.tout[0], np.log2(x2.observables["Obs_All"]), '0.5', lw = 2)
    #         # plt.axvline(x = 168, linewidth = 4, color = 'red')
    #         # plt.axvline(x = 336, linewidth = 4, color = 'blue')
    #         plt.xlabel("Time (hours)")
    #         plt.ylabel("log2(cell count)")
    #         plt.xlim([0,350])
    #         plt.title("LSD Spike-In Normal Population")
            # post_drug_counts.append(x2.all["Obs_All"][-1])
            # print(x1.all[ind]["Obs_Cell0"])

    # for ind, cell in enumerate(preDrug_cellResistant):
    #     print("%d:%d" % (exp, ind))
    #     if cell != 0:
    #         print((k_div_resistant.value-k_death_resistant.value)*np.ones(int(cell)))
    #         # Histogram (line) of Resistant DIP Rates
    #         plt.figure("DIP Rate Distribution Resistant #%d" % ind)
    #         sns.distplot((k_div_resistant.value - k_death_resistant.value) * np.ones(int(cell)), color = '0.5', kde = False,
    #                      label="n_cells = %d" % int(cell))
    #         # plt.bar(dips - 0.5 * (dips[1] - dips[0]), 1. * cell_pop / int(round(cell)), dips[1] - dips[0],
    #         #         label="n_cells = %d" % int(cell))  # *100 for true normal
    #         # plt.plot(dips, normal_hist, lw=2, color="r")
    #         plt.legend(loc=0)
    #
    #         t3 = np.linspace(0, 168, 169)  # 7 days in drug
    #
    #         x3 = sim.run(tspan=t3, param_values={"k_death_0": model.parameters['k_death_0'].value},
    #                      n_sim=1, initials={model.species[i]: pop for i, pop in enumerate(cell_pop)},
    #                      verbose=False)
    #
    #         # post_drug_counts.append(x2.all["Obs_All"][-1])
    #         print(x1.all[ind]["Obs_Cell0"])
                    # plt.figure()
                    # plt.plot(x1.tout[ind], np.log2(x1.all[ind]["Obs_Cell0"]), '0.5', lw = 2)
                    # a = np.array(x1.observables)["Obs_Cell0"]
                    # print(np.array([tr[:] for tr in a]))
                    # # plt.plot(x1.tout, np.array([tr[:] for tr in a]), '0.5')
                    #
                    # plt.plot(t1[-1] + x2.tout[0], np.log2(x2.observables["Obs_All"]), '0.5', lw = 2)
                    # plt.axvline(x = 336, linewidth = 4, color = 'red')
                    # plt.axvline(x = 504, linewidth = 4, color = 'blue')
                    # plt.xlabel("Time (hours)")
                    # plt.ylabel("log2(cell count)")
                    # plt.xlim([0,550])
                    # plt.title("Barcoding Model - 100 barcodes, 10 experiments, 15 states")



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
        #     print(post_drug_counts)
        #     all_post_drug_counts.append(post_drug_counts)
        #
        # print(all_post_drug_counts)
        #             post_drug_counts.append(0.)
        #         else:





                    # ,"cellInit_0": 0}, - FIXED THIS PROBLEM


#          # Same loop for resistant population
#
#         for ind, cell in enumerate(cell_tot_resistant):
#             print("%d:%d" % (exp, ind))
#             if cell == 0:
#                 post_drug_counts.append(0.)
#             else:
#                 # print(cell_pop)
#
#                 # Normal Histogram Multinomial Selection of DIP Rates
#                 # plt.figure("DIP Rate Distribution Barcode #%d" % ind)
#                 # plt.bar(dips-0.5*(dips[1]-dips[0]), 1.*cell_pop/int(round(cell)), dips[1]-dips[0],
#                 #         label = "n_cells = %d" %int(cell)) # *100 for true normal
#                 # plt.plot(dips, normal_hist, lw = 2, color = "r")
#                 # plt.legend(loc = 0)
#
#                 t3 = np.linspace(0, 168, 169)  # 7 days in drug
#
#                 x3 = sim.run(tspan=t3, param_values={"k_death_0": model.parameters['k_death_0'].value},
#                              n_sim=1, initials={model.species[i]: pop for i, pop in enumerate(cell_pop)},
#                              verbose=False)
#                 # ,"cellInit_0": 0}, - FIXED THIS PROBLEM
#
#                 post_drug_counts.append(x2.all["Obs_All"][-1])
#                 print(x1.all[ind]["Obs_Cell0"])
#                 plt.figure()
#                 plt.plot(x1.tout[ind], np.log2(x1.all[ind]["Obs_Cell0"]), '0.5', lw=2)
#                 a = np.array(x1.observables)["Obs_Cell0"]
#                 print(np.array([tr[:] for tr in a]))
#                 # quit()
#                 # plt.plot(x1.tout, np.array([tr[:] for tr in a]), '0.5')
#                 plt.plot(t1[-1] + x2.tout[0], np.log2(x2.observables["Obs_All"]), '0.5', lw = 2)
#                 plt.axvline(x = 336, linewidth = 4, color = 'red')
#                 plt.axvline(x = 504, linewidth = 4, color = 'blue')
#                 plt.xlabel("Time (hours)")
#                 plt.ylabel("log2(cell count)")
#                 plt.xlim([0,550])
#                 plt.title("Barcoding Model - 100 barcodes, 10 experiments, 15 states")
#                 # plt.figure("Experiment #%d Barcode #%d" % (exp,ind))
#                 # plt.plot(t1, x1.all[int(cell)]["Obs_Cell0"],'0.5', lw=4, alpha=0.25)
#                 # for i,obs in enumerate(model.observables):
#                     # print(x1.all[int(cell)])
#                     # quit()
#                     # plt.plot(t1, x1.all[obs.name], lw = 2, label = obs.name, color = colors[i])
#                     # plt.plot(t1, x1.all[obs.name], '0.5', lw=4, alpha=0.25)
#                     # plt.plot(t1[-1]+t2, x2.all[obs.name], '0.5', lw = 4, alpha = 0.25)
#                     # plt.plot(t1, x1.all[obs.name], lw = 2, label = obs.name, color = colors[i])
#                     # plt.plot(t1[-1]+t2, x2.all[obs.name], lw = 2, label = obs.name)
#                 # plt.annotate(s = "n_cells = %d" % cell, xy = (0.1,0.1), xycoords = "axes fraction", fontsize = 24)
#                 # plt.text(-0.04, 2800, r'$\mu=%g$' % (np.mean(distr_all)))
#                 # plt.legend(loc = 0)
#         print(post_drug_counts)
#         all_post_drug_counts.append(post_drug_counts)
#
#     print(all_post_drug_counts)
#     # np.save("barcoding_data_%dbar%dexp%dstates.npy" % (b, e, s), all_post_drug_counts)
#     # np.save("barcoding_data_100bar10exp15states_forpres.npy", all_post_drug_counts)
#
#     plt.show()
#

# states = [15]
# experiments = range(1,11)
# barcodes = [5]
#
# for s in states:
#     for e in experiments:
#         for b in barcodes:
#             run_barcoding_model(e,b,s)
