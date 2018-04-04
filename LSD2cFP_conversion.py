from pysb import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sp
from pysb.simulator import ScipyOdeSimulator
from pysb.simulator.bng import BngSimulator
import pandas as pd
import seaborn as sns

sns.set(font_scale = 1.25)
sns.set_style("whitegrid")

# Setting Key Model Parameters
n_exp = 384
# n_barcodes = 5
n_cell_types = 3
n_cells = 5
cols = ['r','b','g','y','k']


def adjusted_r_squared(r, n, p):
    """
    Calculate adjusted r-squared value from r value
    Parameters
    ----------
    r: float
        r value (between 0 and 1)
    n: int
        number of sample data points
    p: int
        number of free parameters used in fit
    Returns
    -------
    float
        Adjusted r-squared value
    """
    if n <= p:
        return np.nan
    return 1 - (1 - r ** 2) * ((n - 1) / (n - p - 1))

def tyson1(adj_r_sq, rmse, n):
    """
    Tyson1 algorithm for selecting optimal DIP rate fit
    Parameters
    ----------
    adj_r_sq: float
        Adjusted r-squared value
    rmse: float
        Root mean squared error of fit
    n: int
        Number of data points used in fit
    Returns
    -------
    float
        Fit value (higher is better)
    """
    return adj_r_sq * ((1 - rmse) ** 2) * ((n - 3) ** 0.25)


def expt_dip(t_hours, assay_vals, selector_fn=tyson1):
    # t_hours = np.array(df_timecourses.index.get_level_values(
    #     level='timepoint').total_seconds()) / SECONDS_IN_HOUR
    #
    # assay_vals = np.log2(np.array(df_timecourses))
    n_total = len(t_hours)

    dip = None
    final_std_err = None
    first_timepoint = None
    final_intercept = None
    dip_selector = -np.inf
    if n_total < 3:
        return None
    for i in range(n_total - 2):
        x = t_hours[i:]
        y = assay_vals[i:]
        slope, intercept, r_value, p_value, std_err = sp.linregress(x, y)

        n = len(x)
        adj_r_sq = adjusted_r_squared(r_value, n, 1)
        predictions = np.add(np.multiply(x, slope), intercept)
        rmse = np.linalg.norm(predictions - y) / np.sqrt(n)
        new_dip_selector = selector_fn(adj_r_sq, rmse, n)
        if new_dip_selector > dip_selector:
            dip_selector = new_dip_selector
            dip = slope
            final_std_err = std_err
            first_timepoint = x[0]
            final_intercept = intercept

    return dip, final_std_err, first_timepoint, final_intercept



def LSD(prob, n_exp, n_cell_types, n_cells):

    dip_dist = []

    low_dips = -0.06 * np.log(2)
    high_dips =  0.01 * np.log(2)
    dips = np.linspace(low_dips, high_dips, n_cell_types)
    dip_mean = -0.02 * np.log(2)
    dip_var = 0.01 * np.log(2)

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
    [Monomer("Cell", ['dip'], {'dip': ["%d" % i for i in range(n_cell_types)]})]
    print(model.monomers)

    Parameter('cellInit_0', 1)
    [Parameter('cellInit_%d' % i) for i in range(1, n_cell_types)]
    # Parameter('Cell_1init')
    print(model.parameters)

    Initial(Cell(dip="0"), cellInit_0)  # could not be a string - parameter only - didn't know how to set up
    [Initial(Cell(dip="%d" % i), model.parameters["cellInit_%d" % i]) for i in range(1, n_cell_types)]
    print(model.initial_conditions)

    [Observable("Obs_Cell%d" % i, Cell(dip="%d" % i)) for i in range(n_cell_types)]
    Observable("Obs_All", Cell())
    print(model.observables)

    k_div = 0.025 * np.log(2)
    [Parameter("k_div_%d" % i, k_div) for i in range(len(dips))]
    k_death = -dips
    [Parameter("k_death_%d" % i, k) for i, k in enumerate(k_death)]
    print(model.parameters)

    [Rule("Cell%d_Div" % i, Cell(dip="%d" % i) >> Cell(dip="%d" % i) + Cell(dip="%d" % i),
          model.parameters["k_div_%d" % i]) for i in range(len(dips))]
    [Rule("Cell%d_Death" % i, Cell(dip="%d" % i) >> None,
          model.parameters["k_death_%d" % i]) for i in range(len(k_death))]
    print(model.rules)

    for exp in range(n_exp):

        num_cells = np.random.poisson(n_cells)
        picks = np.random.multinomial(num_cells, prob)
        picks_total = np.sum(picks)
        if picks_total != 0:
            # ['cellInit_%d' % i for i in range(n_cell_types)] =
            cellInit_0.value = picks[0]
            cellInit_1.value = picks[1]
            cellInit_2.value = picks[2]

            # A0.value = A_init
            # B0.value = B_init
            # C0.value = C_init

            k_div_0.value = .025 * np.log(2)
            k_div_1.value = .025 * np.log(2)
            k_div_2.value = .025 * np.log(2)

            k_death_0.value = .002 * np.log(2)
            k_death_1.value = .002 * np.log(2)
            k_death_2.value = .002 * np.log(2)

            # print cellInit_0, cellInit_1, cellInit_2

        t1 = np.linspace(0, 169, 168)  # 7 days
        sim1 = ScipyOdeSimulator(model, tspan=t1, verbose=False)
        x1 = sim1.run()  # returns np.array with species and obs
        # plt.plot(x1.tout[0], np.log2(x1.all["Obs_All"]), color = cols[exp], lw=2)
        t2 = np.linspace(0, 336, 337)
        sim2 = ScipyOdeSimulator(model, tspan=t2, verbose=False)
        x2 = sim2.run(param_values={"cellInit_0": x1.all["Obs_Cell0"][-1],
                                    "cellInit_1": x1.all["Obs_Cell1"][-1],
                                    "cellInit_2": x1.all["Obs_Cell2"][-1],
                                    "k_death_0": k_death[0],
                                    "k_death_1": k_death[1],
                                    "k_death_2": k_death[2]})

        # print x1.all["Obs_Cell0"][-1]
        # print x1.all["Obs_Cell1"][-1]
        # print x1.all["Obs_Cell2"][-1]

        print(expt_dip(t_hours=x2.tout[0], assay_vals=np.log2(x2.all["Obs_All"])))

        dip,std_err,first_tp,intercept = expt_dip(t_hours=x2.tout[0], assay_vals=np.log2(x2.all["Obs_All"]))

        dip_dist.append(dip)

        plt.plot(x1.tout[0], np.log2(x1.all["Obs_All"]), '0.5', lw=2)
        plt.plot(x1.tout[0][-1] + x2.tout[0], np.log2(x2.all["Obs_All"]), '0.5', lw=2)
        plt.xlabel("Time (hours)")
        plt.ylabel("Population Doublings")
        plt.title("Deterministic LSD trajectories", weight = "bold")

    plt.figure("DIP Rate Distribution")
    sns.distplot(dip_dist, kde=False, color='k', hist_kws={"alpha":0.5})
    plt.xlabel("DIP Rate")
    plt.ylabel("Frequency")
    plt.title("LSD Biased DIP Rate Distribution", weight = "bold")

LSD(prob=[0.2,0.5,0.3], n_cell_types=3, n_exp=100, n_cells=50)
plt.show()


""""" NEXT STEPS
1. Feed LSD distribution into numerical solver for three states
2. See if can reproduce input cFP distribution
3. Repeat and extend to more states
"""""

