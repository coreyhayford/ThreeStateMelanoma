from pysb import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sp
from pysb.simulator.bng import BngSimulator
import seaborn as sns
import math
from scipy.stats import sem

sns.set(font_scale = 1.25)
sns.set_style("whitegrid")


def seg_model(num_cells, k_death, k_division, n_sims, color):

    def plot_mean_min_max(name, title=None, color="red"):
        trajectories = y1
        x = np.array([tr[:][name] for tr in trajectories]).T
        tout = x1.tout.T
        plt.title(title)
        plt.subplot(111)
        mean = np.log2(x.mean(1) / x.mean(1)[0])
        print(np.log2(x / x[0]))
        ci_95 = 1.96 * sem(np.log2(x / x[0]).T)
        # print(ci_95)
        ci_min = mean - ci_95
        ci_max = mean + ci_95
        # print(ci_min)
        # print(ci_max)
        # plt.fill_between(t1, ci_min, ci_max, facecolor=color, alpha=0.25)
        plt.plot(tout, np.log2(x / x[0]), '0.5', color=color, lw=1, alpha=0.25)  # individual trajectories
        plt.ylim(0,20)
        plt.legend(loc="upper left")
        plt.xlabel("Time (hours)", weight="bold")
        plt.ylabel("log2(Cell Count)", weight="bold")
        plt.title("Model Trajectories", weight="bold")
        # plt.savefig('stochasticGrowth_%dcells_%.3fdivrate_%.3fdivrate_%dsims.pdf' % (num_cells, k_death, k_division, n_sims))

    Model()

    Monomer("Cell")
    Parameter('CellInit', num_cells)
    Initial(Cell(), CellInit)
    Observable("Obs_Cell", Cell())
    Parameter("k_div", k_division)
    Parameter("k_dth", k_death)
    Rule("Cell_Div", Cell() >> Cell() + Cell(), k_div)
    Rule("Cell_Death", Cell() >> None, k_dth)

    t1 = np.linspace(0, 336, 337)  # 14 days
    sim = BngSimulator(model, tspan=t1, verbose=False)  # , seed = 1095205711)
    x1 = sim.run(n_runs=n_sims, verbose=False)  # returns np.array with species and obs
    y1 = np.array(x1.observables)

    plt.figure()
    plot_mean_min_max('Obs_Cell', color=color)
    # for sim in range(n_sims):
        # plt.plot(x1.tout[sim], x1.all[sim]["Obs_Cell"], lw=2,
        #          label="Simulated Data" if sim == 0 else "")
        # plt.plot(x1.tout[sim], (np.log2(x1.all[sim]["Obs_Cell"]/x1.all[sim]["Obs_Cell"][0])), lw=2,
        #          label = "Simulated Data" if sim == 0 else "")

        # slope, intercept, r_value, p_value, std_err = sp.stats.linregress(x1.tout[sim][80:], np.log2(
        #     x1.all[sim]["Obs_Cell"][80:]/x1.all[sim]["Obs_Cell"][0]))
        # if math.isnan(slope) == False:
        #     dip_dist = dip_dist + [slope]
        # plt.plot(x1.tout[sim], intercept + slope * x1.tout[sim], 'k', label = "Fitted Lines" if sim == 0 else "", lw = 2)
        #
        # slope, intercept, r_value, p_value, std_err = sp.stats.linregress(x1.tout[sim], np.log2(
        #     x1.all[sim]["Obs_Cell"] / x1.all[sim]["Obs_Cell"][0]))
        # if math.isnan(slope) == False:
        #     dip_dist1 = dip_dist1 + [slope]
        # plt.plot(x1.tout[sim], intercept + slope * x1.tout[sim], 'k', label="Fitted Lines" if sim == 0 else "", lw=2)



seg_model(num_cells=32, k_division = 0.035, k_death=0.005, n_sims=100, color="red")
# seg_model(num_cells=16, k_division = 0.035, k_death=0.005, n_sims=100, color="green")
# seg_model(num_cells=8, k_division = 0.035, k_death=0.005, n_sims=100, color="blue")
# seg_model(num_cells=4, k_division = 0.035, k_death=0.005, n_sims=100, color="magenta")
# seg_model(num_cells=2, k_division = 0.035, k_death=0.005, n_sims=100, color="yellow")
# seg_model(num_cells=1, k_division = 0.035, k_death=0.005, n_sims=100, color="black")



# cell_numbers = [32,16,8,4,2,1]
# div_rates = [0.03]
# dth_rates = [0.005]#, 0.01, 0.02, 0.03, 0.04, 0.05, 0.055]
# n_exps = [100]#, 1000, 10000]
#
# for cn in cell_numbers:
#     for di_r in div_rates:
#         for de_r in dth_rates:
#             for ns in n_exps:
#                 seg_model(num_cells=cn, k_death=de_r, k_division=di_r, n_sims=ns)

plt.show()