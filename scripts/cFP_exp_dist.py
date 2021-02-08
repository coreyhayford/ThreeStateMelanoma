# Do numerical integral, delta x small enough, run loop and get distance between curves and abs value it
import numpy as np
# import scipy.integrate as integ
from sympy import DiracDelta
from sympy import *
# from sympy.mpmath import *
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as sp
import pandas as pd
from numpy import median
from scipy.stats import gaussian_kde
from scipy.integrate import quad
import math

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

sns.set(font_scale = 1.25)
sns.set_style("white")
# print(test.index)

# test = pd.read_csv('PC9-VU-99_PC9-BR1-BFP-1_Erl-3uM_25cells_nl2counts.csv', index_col = 0)
#
# # exp1.loc[103.2:] = exp1.loc[103.2:] - (exp1.loc[103.2] - exp1.loc[89.4])
# # test.plot()
# # plt.show()
# # print(list(test['X144'][test.index[6]:]))
# # print(list(test.index[6:]))
# # print(test[5][])
# # quit()
# dip_dist = []
# for col in test:
#     # print(exp1_nl2.index[start_index:])
#     # print(exp1_nl2[col][start_time:])
#     dip, final_std_err, first_timepoint, final_intercept = expt_dip(list(test.index[6:]), list(test[col][test.index[6]:]))
#     if math.isnan(dip) == False:
#         dip_dist = dip_dist + [dip]
# print(dip_dist)
#
# sns.distplot(dip_dist, kde = True, color="k", bins = 50)
# plt.xlabel("DIP Rate", weight = "bold")
# plt.ylabel("Density", weight = "bold")
# plt.title("PC9.VU 99% | PC9.BR1 1% LSD", weight = "bold")
# plt.legend(title = "Experiments", loc = "upper left")
# plt.savefig("ExperimentalLSD_PC9-VU-99_PC9-BR1-01_25cells.pdf")
# plt.show()
# quit()

def get_LSD_data(file_name, start_index1, start_time1, start_index2, start_time2):

    # Getting number of columns to loop over later
    with open(file_name) as f:
        ncols = len(f.readline().split(','))

    # Reading data for generating time points and getting cell counts
    data_all = np.genfromtxt(file_name, delimiter=",")
    data = np.genfromtxt(file_name, delimiter=",", usecols=range(1,ncols))

    # Identifying time points
    times = []
    for tp in range(len(data_all)):
        time_pt = data_all[tp][0]
        times.append(time_pt)

    # Identifying data for each time point and removing NaNs
    ## NaNs exist because csv exported data with multiple replicates are shifted in file
    data_noNaN = []
    for pt in range(len(data)):
        data_sub = [value for value in data[pt] if not math.isnan(value)]
        data_noNaN.append(data_sub)

    # Removing NaNs from time list
    times = [value for value in times if not math.isnan(value)]

    # Identifying drug concentration and cell counts
    drug_conc = data_noNaN[0]
    cell_counts = data_noNaN[1:]

    # print(times)
    # print(drug_conc)
    # print(cell_counts)

    # Breaking down data by replicate
    times_exp1 = times[0::2]
    times_exp2 = times[1::2]
    drug_conc_each = drug_conc[0::2]
    cell_counts_exp1 = cell_counts[0::2]
    cell_counts_exp2 = cell_counts[1::2]

    # print(times_exp1)
    # print(drug_conc_each)
    # print(cell_counts_exp1)
    #
    # print(times_exp2)
    # print(drug_conc_each)
    # print(cell_counts_exp2)

    # Establishing data as pandas dataframe
    exp1 = pd.DataFrame(cell_counts_exp1)
    exp1.index = times_exp1
    # print(exp1)

    exp2 = pd.DataFrame(cell_counts_exp2)
    exp2.index = times_exp2
    # print(exp2)
    #
    # exp1.plot()
    # # plt.show()
    # # quit()
    # exp2.plot()
    # plt.show()
    #
    # # Normalizing each post-drug data to pre-drug time point
    exp1.loc[103.2:] = exp1.loc[103.2:] - (exp1.loc[103.2] - exp1.loc[89.4])
    exp2.loc[103.6:] = exp2.loc[103.6:] - (exp2.loc[103.6] - exp2.loc[89.7])
    #
    # # Normalizing each count trajectory to time point 1 and plotting
    exp1_nl2 = np.log2(exp1.loc[:,:].div(exp1.iloc[0][:]))
    # exp1_nl2.plot(legend = False)
    # plt.xlabel("Time (hours)", weight="bold")
    # plt.ylabel("Population Doublings", weight="bold")
    # plt.title("LSD - PC9-MGHd in 3uM Erlotnib - Rep 1", weight="bold")
    #
    exp2_nl2 = np.log2(exp2.loc[:,:].div(exp2.iloc[0][:]))
    # exp2_nl2.plot(legend = False)
    # plt.xlabel("Time (hours)", weight="bold")
    # plt.ylabel("Population Doublings", weight="bold")
    # plt.title("LSD - PC9-MGHd in 3uM Erlotnib - Rep 2", weight="bold")
    #
    # plt.show()
    # print("XXX")
    # print(exp1_nl2[56.6:])
    # print(exp2_nl2[57.3:])
    # print(exp1_nl2.index[start_index:])
    # print(exp1_nl2[1][start_time:])
    #
    # Generating DIP rate distributions for each replicate
    dip_dist1 = []
    for col in exp1_nl2:
        # print(exp1_nl2.index[start_index:])
        # print(exp1_nl2[col][start_time:])
        slope, intercept, r_value, p_value, std_err = sp.stats.linregress(exp1_nl2.index[start_index1:],
                                                                          exp1_nl2[col][start_time1:])
        if math.isnan(slope) == False:
            dip_dist1 = dip_dist1 + [slope]
    # print(dip_dist1)
    #
    dip_dist2 = []
    for col in exp2_nl2:
        # print(exp1_nl2.index[start_index:])
        # print(exp1_nl2[col][start_time:])
        slope, intercept, r_value, p_value, std_err = sp.stats.linregress(exp2_nl2.index[start_index2:],
                                                                          exp2_nl2[col][start_time2:])
        if math.isnan(slope) == False:
            dip_dist2 = dip_dist2 + [slope]
            # print(dip_dist2)
    return dip_dist1 + dip_dist2


cFP_rates = pd.read_csv("Frick_cFP_rates.csv")

cFP_MGH_DMSO = cFP_rates[(cFP_rates['cell.line'] == 'PC9.MGH') & (cFP_rates['drug1'] == 'control')]
cFP_MGH_Erl = cFP_rates[(cFP_rates['cell.line'] == 'PC9.MGH') & (cFP_rates['drug1'] == 'erlotinib')]

LSD_MGH_Erl = get_LSD_data(file_name='LSD_HTS016-MGHd.csv', start_time1=48.3, start_index1=3,
                           start_time2=48.7, start_index2=3)

# np.save('LSD_MGH_Erl.npy', LSD_MGH_Erl)
# np.save('cFP_MGH_Erl.npy', list(cFP_MGH_Erl['rate']))
print(type(LSD_MGH_Erl))

# print(cFP_rates)

print(type(list(cFP_MGH_DMSO['rate'])))
print(type(list(cFP_MGH_Erl['rate'])))
# quit()

all_data = LSD_MGH_Erl+list(cFP_MGH_DMSO['rate'])+list(cFP_MGH_Erl['rate'])
min_lim = -0.05
max_lim = 0.05

labels = ['cFP DMSO','cFP Erl', 'LSD Erl']
bins = np.arange(min(all_data), max(all_data) + 0.0005, 0.0005)
plt.figure()
sns.distplot(cFP_MGH_DMSO['rate'], kde = True, color="r", bins = bins)
sns.distplot(cFP_MGH_Erl['rate'], kde = True, color="b", bins = bins)
sns.distplot(LSD_MGH_Erl, kde = True, color = "g", bins = bins)
# plt.xlim(-0.2,0.1)
plt.xlabel("DIP Rate", weight = "bold")
plt.ylabel("Density", weight = "bold")
plt.title("PC9.MGH", weight = "bold")
plt.legend(title = "Experiments", labels = labels, loc = "upper left")
plt.savefig("MGH_cFP-LSD-comparison.pdf")
plt.show()
quit()
# expData_p07 = np.genfromtxt("SKMEL5_parentalP22_slopesCI.txt", delimiter='\t')
# expData_p07_df = pd.DataFrame(expData_p07)
# expData_p07_df.columns = ['Replicate', 'DIP Rates', 'CI_lower', 'CI_upper']
# p07rates = expData_p07_df[['DIP Rates']]
# p07rates_noNull = p07rates.apply(lambda x: [y for y in x if pd.notnull(y)])
# p07rate_values = p07rates_noNull.values[0]
# p07rate_vals = np.asarray(p07rate_values)

# SLOPE DATA with CIs

# expData_cFP1 = np.genfromtxt("SKMEL5par_H2BGFP_cFP_slopes.txt", delimiter='\t')
# expData_cFP2 = np.genfromtxt("SKMEL5par_H2BGFP_cFP2_slopes.txt", delimiter='\t')
#
# expData_cFPall = np.concatenate((expData_cFP1, expData_cFP2))
#
# labels = ["cFP1", "cFP2"]
#
# plt.figure()
# sns.distplot(expData_cFP1, kde = True, color="r", bins = 50)
# sns.distplot(expData_cFP2, kde = True, color="g", bins = 50)
# plt.xlim(-0.2,0.1)
# plt.xlabel("DIP Rate", weight = "bold")
# plt.ylabel("Density", weight = "bold")
# plt.title("SKMEL5 H2B-GFP cFP Comparison - 100 cell number", weight = "bold")
# plt.legend(title = "Experiments", labels = labels, loc = "upper left")
#
# plt.figure()
# sns.distplot(expData_cFPall, kde = False, color="b", bins = 50)
# plt.xlim(-0.2,0.1)
# plt.xlabel("DIP Rate", weight = "bold")
# plt.ylabel("Density", weight = "bold")
# plt.title("SKMEL5 H2B-GFP cFPs Merged - 100 cell number", weight = "bold")

#
# plt.figure()
# sns.distplot(expData_cFP1, hist = False, color="r", kde_kws={"shade": True})
# sns.distplot(expData_cFP2, hist = False, color="g", kde_kws={"shade": True})
# plt.xlim(-0.1,0.05)
# plt.xlabel("DIP Rate", weight = "bold")
# plt.ylabel("Density", weight = "bold")
# plt.title("SKMEL5 H2B-GFP cFP Comparison - 150 cell number", weight = "bold")
# plt.legend(title = "Experiments", labels = labels, loc = "upper left")

# labels1 = ["cFP", "LSD"]
# plt.figure()
# sns.distplot(expData_cFPall, hist = False, color="b", kde_kws={"shade": True})
# sns.distplot(p07rate_vals, hist = False, color = "m", kde_kws={"shade": True})
# # plt.xlim(-0.1,0.05)
# plt.xlabel("DIP Rate", weight = "bold")
# plt.ylabel("Density", weight = "bold")
# plt.title("SKMEL5 H2B-GFP cFP Merged vs LSD", weight = "bold")
# plt.legend(title = "Assay", labels = labels1, loc = "upper left")
# plt.show()