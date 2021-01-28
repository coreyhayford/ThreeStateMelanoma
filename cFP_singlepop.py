from pysb import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sp
from pysb.simulator.bng import BngSimulator
import seaborn as sns
import math

sns.set(font_scale = 1.25)
sns.set_style("whitegrid")

# def get_cmap(n, name='hsv'):
#     '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
#     RGB color; the keyword argument name must be a standard mpl colormap name.'''
#     return plt.cm.get_cmap(name, n)
#
# cmap = get_cmap(10)

colors = ['brown', 'magenta', 'yellow', 'black']

def seg_model(num_cells, k_death, k_division, k_death_drug, n_sims, col):

    # dip_dist = []
    # dip_dist1 = []
    colors = []

    for i in range(n_sims):

        Model()

        Monomer("Cell")
        Parameter('CellInit', num_cells)
        Initial(Cell(), CellInit)
        Observable("Obs_Cell", Cell())
        Parameter("k_div", k_division)
        Parameter("k_dth", k_death)
        Parameter("k_dth_drug", k_death_drug)
        Rule("Cell_Div", Cell() >> Cell() + Cell(), k_div)
        Rule("Cell_Death", Cell() >> None, k_dth)

        t1 = np.linspace(0, 168, 169)  # 7 days
        sim = BngSimulator(model, tspan=t1, verbose=False)  # , seed = 1095205711)
        x1 = sim.run(n_runs=1, verbose=False)  # returns np.array with species and obs

        plt.plot(x1.tout[0], np.log2(x1.all["Obs_Cell"] / x1.all["Obs_Cell"][0]),
                 color = col, lw=2)

        CellInit.value = x1.all["Obs_Cell"][-1]
        print(CellInit.value)
        k_div.value = k_division
        k_dth.value = k_death_drug

        t2 = np.linspace(0, 75, 76)
        sim1 = BngSimulator(model, tspan=t1, verbose=False)
        x2 = sim1.run(nsims=1, verbose=False)

        plt.plot(x1.tout[0][-1] + x2.tout[0],
                 np.log2(x2.all["Obs_Cell"] / x1.all["Obs_Cell"][0]),
                 color = col, lw=2)

    # plt.figure()
    # for sim in range(n_sims):
    #     plt.plot(x1.tout[sim], (np.log2(x1.all[sim]["Obs_Cell"]/x1.all[sim]["Obs_Cell"][0])), '0.5', lw=2,
    #              label = "Simulated Data" if sim == 0 else "")
    #     plt.plot(x1.tout[sim][80:], (np.log2(x1.all[sim]["Obs_Cell"][80:] / x1.all[sim]["Obs_Cell"][0])), 'r', lw=2,
    #              label="Simulated Data" if sim == 0 else "")
    #     print(x1.tout[1][80:])
    #
    #     slope, intercept, r_value, p_value, std_err = sp.stats.linregress(x1.tout[sim][80:], np.log2(
    #         x1.all[sim]["Obs_Cell"][80:]/x1.all[sim]["Obs_Cell"][0]))
    #     if math.isnan(slope) == False:
    #         dip_dist = dip_dist + [slope]
    #     plt.plot(x1.tout[sim], intercept + slope * x1.tout[sim], 'k', label = "Fitted Lines" if sim == 0 else "", lw = 2)
    #
    #     slope, intercept, r_value, p_value, std_err = sp.stats.linregress(x1.tout[sim], np.log2(
    #         x1.all[sim]["Obs_Cell"] / x1.all[sim]["Obs_Cell"][0]))
    #     if math.isnan(slope) == False:
    #         dip_dist1 = dip_dist1 + [slope]
    #     plt.plot(x1.tout[sim], intercept + slope * x1.tout[sim], 'k', label="Fitted Lines" if sim == 0 else "", lw=2)

    plt.legend(loc = "upper left")
    plt.xlabel("Time (hours)", weight = "bold")
    plt.ylabel("Population Doublings", weight = "bold")
    plt.axvline(x = 168, color = 'k', ls = '--', lw = 3)
    plt.title("Model Trajectories", weight = "bold")
    plt.savefig("stochsims_Fig1_10cells_%dcells.pdf" % num_cells)

    # plt.figure()
    # sns.distplot(dip_dist, hist = False, color = "k", kde_kws={"shade": True})
    # sns.distplot(dip_dist1, hist = False, color = "r", kde_kws={"shade": True})
    # plt.axvline(x=(k_division-k_death)/np.log(2), color = 'k', ls = '--', lw = 3)
    # plt.xlabel("DIP Rate", weight = "bold")
    # plt.ylabel("Density", weight = "bold")
    # plt.title("DIP Rate Distribution - 100 cells", weight = "bold")

    # print(dip_dist)
    # np.save("stoch_expgrowth_%dcells_%ddivrate_%ddthrate_%dsims.npy" % (cn,di_r,de_r,ns), dip_dist)

cells = [1,5,10,25]
for i,c in enumerate(cells):
    plt.figure()
    seg_model(num_cells=c, k_division = 0.035, k_death=0.005,
              k_death_drug=0.045, n_sims=10, col=colors[i])

# cell_numbers = [1,5,10,50,100]
# div_rates = [0.03]
# dth_rates = [0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.055]
# n_exps = [100, 1000, 10000]

# for cn in cell_numbers:
#     for di_r in div_rates:
#         for de_r in dth_rates:
#             for ns in n_exps:
#                 seg_model(num_cells=cn, k_death=de_r, k_division=di_r, n_sims=ns)

plt.show()
