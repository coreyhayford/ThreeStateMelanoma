from pysb import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sp
from pysb.simulator.bng import BngSimulator
from pysb.bng import generate_equations
import pandas as pd
import seaborn as sns
import math
import matplotlib


cFP_rates = pd.read_csv("cFP_rates_VUlines.csv")

div = [0.032, 0.033]
dth = [0.0311, 0.0265]

def dist_compare(div, dth):
    dat = cFP_rates[cFP_rates['Cell_Line'] == 'PC9-DS8']
    num_cells = 1
    # kdiv = [0.026, 0.026]
    # kdth = [0.025, 0.020]
    kdiv = div
    kdth = dth



    Model()

    Monomer('Cell', ['type'], {'type': ["x"+str(i) for i in range(len(kdiv))]})

    [Initial(Cell(type="x"+str(i)), Parameter('Cell%d_Init' % i, num_cells)) for i in range(len(kdiv))]

    [Rule('divide_Cell%d' % i, Cell(type="x"+str(i)) >> Cell(type="x"+str(i)) + Cell(type="x"+str(i)),
          Parameter('kdiv_%d' % i, kdiv[i])) for i in range(len(kdiv))]

    [Rule('die_Cell%d' % i, Cell(type="x"+str(i)) >> None,
          Parameter('kdth_%d' % i, kdth[i])) for i in range(len(kdiv))]

    # [degrade(Cell(type=str(i)), Parameter('kdeg_%d' % i, kdeg[i])) for i in range(len(kdiv))]
    Observable("Cell_total", Cell())
    [Observable("Cell_t_%s" % i, Cell(type="x"+str(i))) for i in range(len(kdiv))]


    sim_dist = []
    plt.figure(figsize=[4,3])

    t1 = np.linspace(0, 200, 201) #hrs
    sim1 = BngSimulator(model, tspan=t1, verbose=False)  # , seed = 1095205711)
    x1 = sim1.run(tspan=t1, param_values={'kdiv_%d' % i: 0.04 * np.log(2),
                                          'kdth_%d' % i: 0.005 * np.log(2)},
                 n_runs=len(dat), verbose=False)

    trajs = np.array(np.array([tr[:]["Cell_total"] for tr in np.array(x1.observables)]).T)

    cell_tot = trajs[-1]
    cell_tot_num = len(cell_tot)
    cell_tot_s = np.random.multinomial(cell_tot_num, [3/4.,1/4.])


    prob_list = ([0]*int(cell_tot_s[0])) + ([1]*int(cell_tot_s[1]))

    cell_pop_n = [np.random.multinomial(int(round(cell)), [1,0]) for i,cell in enumerate(cell_tot)]
    cell_pop1 = [np.random.multinomial(int(round(cell)), [1,0]) for i,cell in enumerate(cell_tot)][:int(cell_tot_s[0])]
    cell_pop2 = [np.random.multinomial(int(round(cell)), [0,1]) for i,cell in enumerate(cell_tot)][int(cell_tot_s[0]):]

    trajects = []
    for ind, cell in enumerate(cell_tot[:int(cell_tot_s[0])]):
        print("%d" % ind)

        # cell_pop1 = np.random.multinomial(int(round(cell)), [1,0])
        # cell_pop2 = np.random.multinomial(int(round(cell)), [0,1])
        # print(cell_pop1)
        cell_pop_1n = cell_pop1[ind]
        t2 = np.linspace(0, 225, 226)  # in drug
        x2 = sim1.run(tspan=t2,
                      initials={model.species[i]: pop for i, pop in enumerate(cell_pop_1n)},
                      verbose=False)
        if cell != 0:
            slope, intercept, r_value, p_value, std_err = sp.stats.linregress(t2, np.log2(
                x2.observables["Cell_total"] / x2.observables["Cell_total"][0]))
            if math.isnan(slope) == False:
                sim_dist = sim_dist + [slope]
        if cell > 50:
            trajects.append(np.log2(x2.observables["Cell_total"] / x2.observables["Cell_total"][0]))
        # plt.plot(t1[-1] + x2.tout[0], np.log2(x2.observables["Cell_total"] / x2.observables["Cell_total"][0]),
        #          color = 'grey', lw=0.5)

    for ind, cell in enumerate(cell_tot[int(cell_tot_s[0]):]):
        print("%d" % ind)

        # cell_pop1 = np.random.multinomial(int(round(cell)), [1,0])
        # cell_pop2 = np.random.multinomial(int(round(cell)), [0,1])
        # print(cell_pop1)
        cell_pop_2n = cell_pop2[ind]
        t3 = np.linspace(0, 225, 226)  # in drug
        x3 = sim1.run(tspan=t3,
                      initials={model.species[i]: pop for i, pop in enumerate(cell_pop_2n)},
                      verbose=False)
        if cell != 0:
            slope, intercept, r_value, p_value, std_err = sp.stats.linregress(t3, np.log2(
                x3.observables["Cell_total"] / x3.observables["Cell_total"][0]))
            if math.isnan(slope) == False:
                sim_dist = sim_dist + [slope]
        if cell > 50:
            trajects.append(np.log2(x3.observables["Cell_total"] / x3.observables["Cell_total"][0]))
        # plt.plot(t1[-1] + x3.tout[0], np.log2(x3.observables["Cell_total"] / x3.observables["Cell_total"][0]),
        #          color = 'grey', lw=0.5)


    # plt.xlabel("Time (hours)")
    # plt.ylabel("Norm Log2 Cell Count")
    # plt.ylim(-5.25,4)
    # plt.title("Post-drug Dynamics - DS8", fontsize = 10)
    # plt.text(210, -3.5, r"$k_{division1}=%.3f$" % kdiv[0], fontsize = 9)
    # plt.text(210, -4, r"$k_{division2}=%.3f$" % kdiv[1], fontsize = 9)
    # plt.text(210, -4.5, r"$k_{death1}=%.4f$" % kdth[0], fontsize = 9)
    # plt.text(210, -5, r"$k_{death2}=%.4f$" % kdth[1], fontsize = 9)
    # plt.savefig("PostDrug_withMetrics_DS8_normalized.svg")


    D_stat, p_val = sp.ks_2samp(dat['DIP_Rate'], sim_dist)
    # mean_exp = np.mean(dat['DIP_Rate'])
    # sd_exp = np.std(dat['DIP_Rate'])
    # mean_sim = np.mean(sim_dist)
    # sd_sim = np.std(sim_dist)

    return trajects, sim_dist, p_val

data_DS8 = dist_compare([0.032, 0.033],[0.0311, 0.0265])
trajectories_DS8 = pd.DataFrame(data_DS8[0])
trajectories_DS8 = trajectories_DS8.transpose()
trajectories_DS8.to_csv('trajectories_DS8_G50.csv')

distributions_DS8 = pd.DataFrame({'DS8': data_DS8[1]})
distributions_DS8.to_csv('distributions_DS8_G50.csv')

# data_DS8[2].to_csv('pval_DS8.csv')
print(data_DS8[2])
quit()


plt.figure(figsize=[4,3])
sns.distplot(sim_dist, kde = True, hist= False, color = "grey", kde_kws={"shade": True})
sns.distplot(dat['DIP_Rate'], kde = True, hist= False, color = "seagreen", kde_kws={"shade": True})
plt.xlabel("DIP Rate")
plt.ylabel("Density")
plt.xlim(-0.025, 0.025)
plt.ylim(0,200)
plt.legend(labels=['Simulated','Experimental'], loc = 'upper left')
# plt.text(0.007, 200, "Mean (Exp)=%.3f" % round(mean_exp, 3), fontsize = 8)
# plt.text(0.007, 180, "SD (Exp)=%.4f" % round(sd_exp, 4), fontsize = 8)
# plt.text(0.007, 160, "Mean (Sim)=%.3f" % round(mean_sim, 3), fontsize = 8)
# plt.text(0.007, 140, "SD (Sim)=%.4f" % round(sd_sim, 4), fontsize = 8)
# plt.text(0.007, 180, r"$k_{division1}=%.3f$" % kdiv[0], fontsize = 9)
# plt.text(0.007, 160, r"$k_{division2}=%.3f$" % kdiv[1], fontsize = 9)
# plt.text(0.007, 140, r"$k_{death1}=%.4f$" % kdth[0], fontsize = 9)
# plt.text(0.007, 120, r"$k_{death2}=%.4f$" % kdth[1], fontsize = 9)
plt.text(0.01, 150, "p=%.3f" % round(p_val, 3), fontsize = 11)
# plt.title("DIP Rate Distribution")
plt.savefig("DS8_model_example_FINAL.svg")

# dist_compare(div = [0.032, 0.033], dth = [0.0311, 0.0265])

# di1 = np.linspace(0.024,0.027, 4)
# di2 = np.linspace(0.024,0.027, 4)
# dip1 = np.linspace(-0.001, 0.001, 3)
# dip2 = np.linspace(0.005, 0.007, 3)
#
# ls_div = []
# ls_dip = []
# ls_dth = []
# p_vals_DS8 = []
# z = 0
# for d1,d1_val in enumerate(di1):
#     for d2, d2_val in enumerate(di2):
#         for p1, p1_val in enumerate(dip1):
#             for p2, p2_val in enumerate(dip2):
#                 ls_div.append([d1_val, d2_val])
#                 ls_dip.append([p1_val,p2_val])
#                 ls_dth.append([d1_val-p1_val, d2_val-p2_val])
#                 p = dist_compare(div = [d1_val, d2_val], dth = [d1_val-p1_val, d2_val-p2_val])
#                 p_vals_DS8.append(p)
#                 z = z + 1
#                 print(z)
#                 # print("%d;%d;%d;%d" % (d1,d2,p1,p2))
#
# dict = {'DIP rate': ls_dip, 'division rate': ls_div,
#         'death rate': ls_dth, 'p-value': p_vals_DS8}
#
# df = pd.DataFrame(data=dict)
#
# df.to_pickle('PC9-DS8_param-scan_twoState.pkl')
