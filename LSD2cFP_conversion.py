from pysb import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sp
from pysb.simulator import ScipyOdeSimulator
from pysb.simulator.bng import BngSimulator
import pandas as pd
import seaborn as sns
from sympy.solvers.solvers import nsolve, Symbol
from sympy.solvers.solveset import nonlinsolve
import sympy
from collections import Counter, OrderedDict
from sympy import log, Eq, pprint
from sympy import symbols, Function, Idx
import mpmath
import pprint


sns.set(font_scale = 1.25)
sns.set_style("whitegrid")
# cols = ['r','b','g','k','y']
# bins = np.linspace(-0.05,0.05,10)
# labels = ['1 cell', '1 cell(ish)', '5 cells', '10 cells', '25 cells']

a = np.load("dip_dist_list_ODE.npy")
b,c,d,e = a[0],a[1],a[2],a[3]
f = list(np.load("dip_dist_1cell_ODE.npy"))
# print(f)
# print(b)
# print(c)
# print(d)
# print(e)

lists = []
lists.append(f)
lists.append(b)
lists.append(c)
lists.append(d)
lists.append(e)

# counts = get_biased_data(input_data = np.array(lists[2]))
# cFP_counts = get_biased_data(input_data = np.load('cFP_MGH_Erl.npy'), bins=np.array([-0.05, -0.015, -0.005, 0.02]))


### Discretizing data into x groups -- BETTER VERSION ###
def get_biased_data(input_data, num_bins, low_lim, high_lim):
    # input_data = np.load('LSD_MGH_Erl.npy')
    # num_bins = 40
    bins=np.linspace(low_lim, high_lim, num_bins+1)
    inds = pd.Series(np.digitize(input_data, bins))
    keys = inds.value_counts().keys()

    all_keys = range(num_bins)
    all_vals = []
    for key in all_keys:
        if key in keys:
            all_vals.append(inds.value_counts()[key])
        else:
            all_vals.append(0)
    all_vals_prop = np.array(all_vals)/float(np.sum(all_vals))
    counts_df = pd.DataFrame(
        {'bin_num': all_keys,
         'num_events': all_vals,
         'prop_events': all_vals_prop}
    )
    return counts_df


# n = 5 # Number of cells per well
LSD_counts = get_biased_data(np.load('LSD_MGH_Erl.npy'), 40, -0.03, 0.02)
num_states = 40 #len(bins)
# p_all = list(symbols('p0:%d'%num_states))

p = sympy.IndexedBase('p')
i = Symbol('i', integer = True)
j = Symbol('j')
n = Symbol('n')
# m = Symbol('m')
m = 14
# e = Eq(n!=3, 1)

p_all_measured = [LSD_counts['prop_events'][z] for z in range(num_states)]

test = sympy.Sum(p[i]**n, (i, 0, m-1))
test1 = sympy.Sum(p[i]**n, (i, m+1, num_states-1))

# test2 = sympy.Product(sympy.Sum())
print(test.doit())
print(test1.doit())
print(1-(test.doit() + test1.doit()))
quit()
ls = list(map(chr, range(num_states)))
print(ls)
quit()
f = lambda N: sum(p_all[N], ())

eq_p1 = sympy.Sum(p_all[p], (p, min(range(num_states)), max(range(num_states))))

# eq_p1 = mpmath.nsum(lambda p: p_all[p], [min(range(num_states)), max(range(num_states))])
print(eq_p1)
quit()

n = 5 # Number of cells per well
p1 = Symbol('p1')
p2 = Symbol('p2')
p3 = Symbol('p3')



P1_measured = LSD_counts[0]/float(np.sum(LSD_counts))
P2_measured = LSD_counts[1]/float(np.sum(LSD_counts))
P3_measured = LSD_counts[2]/float(np.sum(LSD_counts))

# P1_measured = (1-0.8**5) * 0.7**5 * 0.5**5
# P2_measured = (1-0.5**5) * 0.7**5
# P3_measured = 1-0.7**5

P1 = log(1 - (p2 + p3)**n) + n*log(p1+p2) + n*log(p1+p3) - log(P1_measured)
P2 = log(1 - (p1 + p3)**n) + n*log(p1+p2) - log(P2_measured)
P3 = log(1 - (p1 + p2)**n) - log(P3_measured)

# P1 = (1 - (p2 + p3)**n) * ((p1 + p2)**n * (p1 + p3)**n) - P1_measured
# P2 = (1 - (p1 + p3)**n) * ((p1 + p2)**n) - P2_measured
# P3 = (1 - (p1 + p2)**n) - P3_measured
# P1 = 1 - P2 - P3 - counts[0]/float(np.sum(counts))
P4 = p1 + p2 + p3 - 1

# test = sympy.solve([P1-counts_prob[0], P2-counts_prob[1], P3-counts_prob[2]], (p1,p2,p3))
# print(test)
# quit()

from sympy.core import S

test = nsolve((P1,P2,P3,P4), (p1,p2,p3), (0.2,0.7,0.1), verify = False)
print(test)
print(sum(test))
quit()
# print(nsolve((P1,P2,P3), (p1,p2,p3), (0.2,0.5,0.3)))


# for i,vals in enumerate(lists):
#     sns.distplot(vals, kde=False, color=cols[i], bins=bins, hist_kws={"alpha":0.25}, label=labels[i])
#     plt.xlabel("DIP Rate")
#     plt.ylabel("Frequency")
#     plt.title("LSD Biased DIP Rate Distribution", weight = "bold")
#     plt.legend()
#
# plt.show()

bins = [-0.04,0,0.04]
init_prob = [0.2,0.5,0.3]
plt.bar(bins,init_prob,width=0.02, label = "Truth", alpha = 0.5)
plt.bar(bins,LSD_counts/float(np.sum(LSD_counts)),width=0.02, label = "LSD", alpha = 0.5)
plt.bar(bins,test,width=0.02, label = "PDF converted", alpha = 0.5)
plt.legend(loc = 0)
plt.xlabel("Simulated DIP Rate")
plt.ylabel("Probability")
plt.title("Conversion of Biased Probability Distribution", weight = "bold")
plt.savefig('LSD_to_cFP_3state_ODE.pdf', dpi = 600)
plt.show()

quit()



# Setting Key Model Parameters
# n_exp = 384
# # n_barcodes = 5
# n_cell_types = 3
# n_cells = 5

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

    low_dips = -0.04 * np.log(2)
    high_dips =  0.04 * np.log(2)
    dips = np.linspace(low_dips, high_dips, n_cell_types)
    # dip_mean = -0.01 * np.log(2)
    # dip_var = 0.01 * np.log(2)

    # Discretize normal distribution of dip rates - used in post drug simulation
    # normal = sp.norm.pdf(dips, dip_mean, dip_var)
    # print(normal)
    # sum = 0
    # for i in range(1, n_cell_types):
    #     sum += normal[i] * (dips[i] - dips[i - 1])
    # print(sum)

    # normal_hist = normal * (dips[1] - dips[0])
    # print(normal_hist)
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

    k_div = 0.040 * np.log(2)
    [Parameter("k_div_%d" % i, k_div) for i in range(len(dips))]
    k_death = k_div-dips
    [Parameter("k_death_%d" % i, k) for i, k in enumerate(k_death)]
    print(model.parameters)

    [Rule("Cell%d_Div" % i, Cell(dip="%d" % i) >> Cell(dip="%d" % i) + Cell(dip="%d" % i),
          model.parameters["k_div_%d" % i]) for i in range(len(dips))]
    [Rule("Cell%d_Death" % i, Cell(dip="%d" % i) >> None,
          model.parameters["k_death_%d" % i]) for i in range(len(k_death))]
    # print(model.rules)
    # for i in range(n_cell_types):
    #     print(model.parameters['cellInit_%d' %i])
    # quit()

    np.random.seed(5)

    for exp in range(n_exp):
        # num_cells = np.random.poisson(n_cells)
        num_cells = 1
        picks = np.random.multinomial(num_cells, prob)
        picks_total = np.sum(picks)
        if picks_total != 0:
            div_death = {}
            for i in range(n_cell_types):
                model.parameters['cellInit_%d' %i].value = picks[i]
                div_death["k_div_%d" %i] = 0.04 * np.log(2)
                div_death["k_death_%d" %i] = 0.005 * np.log(2)
        else:
            continue
        # print(model.parameters)
        t1 = np.linspace(0, 169, 168)  # 7 days
        sim = ScipyOdeSimulator(model, verbose=False)
        # sim1 = BngSimulator(model, tspan=t1, verbose=False)
        x1 = sim.run(tspan=t1,
                     param_values=div_death)  # returns np.array with species and obs
        # plt.plot(x1.tout[0], np.log2(x1.all["Obs_All"]), color = 'k', lw=2)
        # plt.plot(x1.tout[0], np.log2(x1.all["Obs_Cell0"]), color = 'b', lw=2)
        # plt.plot(x1.tout[0], np.log2(x1.all["Obs_Cell1"]), color = 'g', lw=2)
        # plt.plot(x1.tout[0], np.log2(x1.all["Obs_Cell2"]), color = 'r', lw=2)
        # print(model.parameters)
        # quit()
        t2 = np.linspace(0, 336, 337)
        # sim2 = BngSimulator(model, tspan=t2, verbose=False)
        # x2 = sim2.run(param_values={"cellInit_{}".format(i): x1.all["Obs_Cell%{}".format(i)][-1] for i in range(n_cell_types),
        #                             "k_death_{}".format(i): k_death[i] for i in range(n_cell_types)})
        # print(model.species)
        # print(x1.species[-1])
        # for r in model.reactions:
        #     print(r)
        x2 = sim.run(tspan=t2,
                     initials=x1.species[-1],
                     param_values=[p.value for p in model.parameters])

        # plt.plot(x2.tout[0]+x1.tout[0][-1], np.log2(x2.all["Obs_All"]), color = 'k', lw=2, label = "All")
        # plt.plot(x2.tout[0]+x1.tout[0][-1], np.log2(x2.all["Obs_Cell0"]), color = 'b', lw=2, label = "Cell0")
        # plt.plot(x2.tout[0]+x1.tout[0][-1], np.log2(x2.all["Obs_Cell1"]), color = 'g', lw=2, label = "Cell1")
        # plt.plot(x2.tout[0]+x1.tout[0][-1], np.log2(x2.all["Obs_Cell2"]), color = 'r', lw=2, label = "Cell2")
        # plt.legend()
        # plt.show()
        # print x1.all["Obs_Cell0"][-1]
        # print x1.all["Obs_Cell1"][-1]
        # print x1.all["Obs_Cell2"][-1]

        print(expt_dip(t_hours=x2.tout[0], assay_vals=np.log2(x2.all["Obs_All"])))
        print(dips)

        dip,std_err,first_tp,intercept = expt_dip(t_hours=x2.tout[0], assay_vals=np.log2(x2.all["Obs_All"]))

        if dip is not None:
            dip_dist.append(dip)

        plt.plot(x1.tout[0], np.log2(x1.all["Obs_All"]), '0.5', lw=2)
        plt.plot(x1.tout[0][-1] + x2.tout[0], np.log2(x2.all["Obs_All"]), '0.5', lw=2)
        plt.xlabel("Time (hours)")
        plt.ylabel("Population Doublings")
        plt.title("Deterministic LSD trajectories", weight = "bold")

    plt.figure("DIP Rate Distribution")
    sns.distplot(dip_dist, kde=False, color='k', hist_kws={"alpha":0.5}, bins = 5)
    plt.xlabel("DIP Rate")
    plt.ylabel("Frequency")
    plt.title("LSD Biased DIP Rate Distribution", weight = "bold")

    return dip_dist

# dd = LSD(prob=[0.2,0.5,0.3], n_cell_types=3, n_exp=384, n_cells=1)
# np.save("dip_dist_1cell_ODE.npy", dd)
# dip_list = []
# num_cells = [1,5,10,25]
# for nc in num_cells:
#     dl = LSD(prob=[0.2,0.5,0.3], n_cell_types=3, n_exp=384, n_cells=nc)
#     dip_list.append(dl)
    # plt.show()

# np.save("dip_dist_list_ODE.npy", dip_list)

""""" NEXT STEPS
1. Feed LSD distribution into numerical solver for three states - X
2. See if can reproduce input cFP distribution - X
3. Repeat and extend to more states
"""""

### Discretizing data into x groups ###
# def get_biased_data(input_data, bins):
#     data_cell = input_data
#     # print(data_cell)
#     bins = bins
#     # bins = np.array([-0.06, -0.02, 0.02, 0.06])
#     inds = np.digitize(data_cell, bins)
#     val_counts = Counter(inds)
#     counts = np.array([val_counts[i+1] for i in range(len(val_counts.items()))])
#     # counts_prob = counts/float(np.sum(counts))
#     # dat = []
#     # for z in val_counts.keys():
#     #     dat.append(data_cell[np.equal(inds, z)])
#     return counts