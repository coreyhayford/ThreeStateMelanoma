# Do numerical integral, delta x small enough, run loop and get distance between curves and abs value it
import numpy as np
# import scipy.integrate as integ
from sympy import DiracDelta
from sympy import *
import math
from sympy.mpmath import *
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp
import pandas as pd
from numpy import median
from scipy.stats import gaussian_kde
from scipy.integrate import quad

sns.set(font_scale = 2)
sns.set_style("whitegrid")

# test1 = np.random.normal(0, 0.1, 10000)
# test2 = np.random.normal(0, 0.1, 10000)
# test3 = np.random.normal(0.1, 0.1, 10000)
#
# sns.distplot(test1, hist = False,color="r", kde_kws={"shade": True})
# sns.distplot(test2, hist = False,color="g", kde_kws={"shade": True})
# sns.distplot(test3, hist = False,color="b", kde_kws={"shade": True})
#
# sns.plt.title("Histogram Distance Comparison", weight = "bold")
# labels = ["Random Normal 1A", "Random Normal 1B", "Random Normal 2"]
# plt.ylabel("Density", weight = "bold")
# plt.xlabel("DIP Rate", weight = "bold")
# # plt.xlim(-0.02,0.02)
# # plt.ylim(0,200)
# plt.legend(labels = labels)
#
#
# plt.show()
#



min = -0.15
max = 0.15
steps = 1000
# x_grid = np.linspace(min, max, steps)
# dx = (max-min) / steps

def bootstrap_resample(X, n=None):
    """ Bootstrap resample an array_like
    Parameters
    ----------
    X : array_like
      data to resample
    n : int, optional
      length of resampled array, equal to len(X) if n==None
    Results
    -------
    returns X_resamples

    Used for experimental data to get idea of "self-distance" (Cao and Petzold, 2006)
    instead of a second ssa simulation from the same PDF
    """
    if n == None:
        n = len(X)
    np.random.seed(123)
    rand_num = np.random.rand(n)
    resample_i = np.floor(rand_num * len(X)).astype(int)
    # print(rand_num)
    # print(resample_i)
    # quit()
    X_resample = X[resample_i]
    return X_resample

def histogram_distance(control_hist, test_hist, name, min = min, max = max, steps = steps):
    # Multiplication factor to ensure values are from 0-1 --> completes "integral"
    dx = (max - min) / steps
    ## Setup grid over which to "integrate" distance difference
    x_grid = np.linspace(min, max, steps)
    dist = 0.5 * dx * sum(abs(test_hist(x_grid) - control_hist(x_grid)))
    # print(str(name), dist)
    return dist

# Random normal different distribution to make sure distance is ~1
# test = np.random.normal(0.07, 0.001, 10000)
# test1 = np.random.normal(0.01, 0.02, 384)
#
# sns.distplot(test1, hist = True, color = "k", kde_kws = {"shade": True})
# plt.xticks([])
# plt.xlabel("DIP Rate", weight = "bold")
# plt.ylabel("Density", weight = "bold")
# plt.title("DIP Rate Distribution - All Wells", weight = "bold")
# plt.show()
# quit()
## Load LCS model produced data - transition model with different initial conditions##
# LCSdata0_trans = np.load("1.0_0.0_0.0_10kData_LCS_forselfdistance.npy")
# LCSdata1_trans = np.load("1.0_0.0_0.0_10kData.npy")
# LCSdata2_trans = np.load("0.95_0.05_0.0_10kData.npy")
# LCSdata3_trans = np.load("0.85_0.1_0.05_10kData.npy")
# LCSdata4_trans = np.load("0.7_0.2_0.1_10kData.npy")
# LCSdata5_trans = np.load("0.5_0.3_0.2_10kData.npy")
# LCSdata6_trans = np.load("0.25_0.5_0.25_10kData.npy")

## Load experimental data - Passage 11, 15, 19, 28, parental SKMEL5 LSD data ## INCORRECT
# expData_p11 = np.loadtxt("experimental_psg11.2.txt", delimiter='\t')
# expData_p15 = np.loadtxt("experimental_psg15.2.txt", delimiter='\t')
# expData_p19 = np.loadtxt("experimental_psg19.2.txt", delimiter='\t')
# expData_p28 = np.loadtxt("experimental_psg28.2.txt", delimiter='\t')
# expData_par = np.loadtxt("experimental_parental2.txt", delimiter='\t')

# SLOPE DATA with CIs

expData_p07_par = np.genfromtxt("SKMEL5-Parental-Passage22_DIPs.txt", delimiter='\t')
expData_Fig1_DIPdist = [value for value in expData_p07_par if not math.isnan(value)]
print(expData_p07_par)
print(type(expData_p07_par))
print(expData_Fig1_DIPdist)
print(type(expData_Fig1_DIPdist))
sns.distplot(expData_Fig1_DIPdist, hist = True, color="k", kde_kws={"shade": True})
plt.xlim(-0.1, 0.1)
plt.xlabel("DIP Rate", weight = "bold")
plt.ylabel("Density", weight = "bold")
plt.title("DIP Rate Distribution - All Wells", weight = "bold")
plt.show()
expData_p07 = np.genfromtxt("SKMEL5_SC01P07_slopesCI.txt", delimiter='\t')
expData_p11 = np.genfromtxt("SKMEL5_SC01P11_slopesCI.txt", delimiter='\t')
expData_p15 = np.genfromtxt("SKMEL5_SC01P15_slopesCI.txt", delimiter='\t')
expData_p19 = np.genfromtxt("SKMEL5_SC01P19_slopesCI.txt", delimiter='\t')
expData_p28 = np.genfromtxt("SKMEL5_SC01P28_slopesCI.txt", delimiter='\t')
expData_par = np.genfromtxt("SKMEL5_parentalP22_slopesCI.txt", delimiter='\t')

expData_p07_df = pd.DataFrame(expData_p07)
expData_p07_df.columns = ['Replicate', 'DIP Rates', 'CI_lower', 'CI_upper']
expData_p11_df = pd.DataFrame(expData_p11)
expData_p11_df.columns = ['Replicate', 'DIP Rates', 'CI_lower', 'CI_upper']
expData_p15_df = pd.DataFrame(expData_p15)
expData_p15_df.columns = ['Replicate', 'DIP Rates', 'CI_lower', 'CI_upper']
expData_p19_df = pd.DataFrame(expData_p19)
expData_p19_df.columns = ['Replicate', 'DIP Rates', 'CI_lower', 'CI_upper']
expData_p28_df = pd.DataFrame(expData_p28)
expData_p28_df.columns = ['Replicate', 'DIP Rates', 'CI_lower', 'CI_upper']
expData_par_df = pd.DataFrame(expData_par)
expData_par_df.columns = ['Replicate', 'DIP Rates', 'CI_lower', 'CI_upper']

rate_data = [expData_p07_df[['DIP Rates']],expData_p11_df[['DIP Rates']],expData_p15_df[['DIP Rates']],
             expData_p19_df[['DIP Rates']],expData_p28_df[['DIP Rates']],expData_par_df[['DIP Rates']]]

rate_data_list = []
for rates in rate_data:
    rates_noNull = rates.apply(lambda x: [y for y in x if pd.notnull(y)])
    print(rates_noNull)
    rate_values = rates_noNull.values[0]
    rate_vals = np.asarray(rate_values)
    rate_data_list.append(rate_vals)

print(rate_data_list)

# data = [expData_p07_df, expData_p11_df, expData_p15_df,
#         expData_p19_df, expData_p28_df, expData_par_df]
#
# DIPs = []
# for dat in data:
#     dat_DIP = dat[['DIP Rates']]
#     DIPs.append(dat_DIP)
# print(DIPs)
# DIPs_cat = pd.concat(DIPs, axis = 1)
# DIPs_cat.columns = ["Passage 07", "Passage 11", "Passage 15", "Passage 19",
#            "Passage 28", "Parental"]
# print(DIPs_cat)

# sns.distplot(DIPs_cat['Passage 07'], hist = False, color="y", kde_kws={"shade": True})
# sns.distplot(DIPs_cat['Passage 11'], hist = False, color="r", kde_kws={"shade": True})
# sns.distplot(DIPs_cat['Passage 15'], hist = False, color="g", kde_kws={"shade": True})
# sns.distplot(DIPs_cat['Passage 19'], hist = False, color="b", kde_kws={"shade": True})
# sns.distplot(DIPs_cat['Passage 28'], hist = False, color="m", kde_kws={"shade": True})
# sns.distplot(DIPs_cat['Parental'], hist = False, color="k", kde_kws={"shade": True})
#
# labels1 = ["Passage 07", "Passage 11", "Passage 15", "Passage 19",
#            "Passage 28", "Parental"]
# sns.plt.title("Comparison Across Passages", weight="bold")
# plt.ylabel("DIP Rate", weight="bold")
# plt.xlabel("Population", weight="bold")
# plt.show()
# quit()

### DIPs range from -0.15 tp +0.15

## Bootstrap resampling of experimental data - used to get 'self-distance' for experimental data ##
## Probably not long term solution...
# expData_p11_resample = bootstrap_resample(expData_p11)
# expData_p15_resample = bootstrap_resample(expData_p15)
# expData_p19_resample = bootstrap_resample(expData_p19)
# expData_p28_resample = bootstrap_resample(expData_p28)
# expData_par_resample = bootstrap_resample(expData_par)

## Load cFP model produced data - transition model with different initial conditions ##
# cFPdata0_trans = np.load("1.0_0.0_0.0_10kData_cFP_forselfdistance.npy")
# cFPdata1_trans = np.load("1.0_0.0_0.0_10kData_cFP.npy")
# cFPdata2_trans = np.load("0.95_0.05_0.0_10kData_cFP.npy")
# cFPdata3_trans = np.load("0.85_0.1_0.05_10kData_cFP.npy")
# cFPdata4_trans = np.load("0.7_0.2_0.1_10kData_cFP.npy")
# cFPdata5_trans = np.load("0.5_0.3_0.2_10kData_cFP.npy")
# cFPdata6_trans = np.load("0.25_0.5_0.25_10kData_cFP.npy")

## Plotting LCS experimental data
# sns.distplot(test, hist = False,color="c", kde_kws={"shade": True}, label = "Different")
# sns.distplot(expData_p11, hist = False,color="k", kde_kws={"shade": True}, label = "Passage 11")
# sns.distplot(expData_p11_resample, hist = False,color="y", kde_kws={"shade": True}, label = "Passage 11 Resample")
# sns.distplot(expData_p15, hist = False,color="r", kde_kws={"shade": True}, label = "Passage 15")
# sns.distplot(expData_p19, hist = False,color="g", kde_kws={"shade": True}, label = "Passage 19")
# sns.distplot(expData_p28, hist = False,color="b", kde_kws={"shade": True}, label = "Passage 28")
# sns.distplot(expData_par, hist = False,color="m", kde_kws={"shade": True}, label = "Parental")

## Plotting LCS model data
# sns.distplot(LCSdata0_trans, hist = False,color="k", kde_kws={"shade": True})
# sns.distplot(LCSdata1_trans, hist = False,color="r", kde_kws={"shade": True})
# sns.distplot(LCSdata2_trans, hist = False,color="g", kde_kws={"shade": True})
# sns.distplot(LCSdata3_trans, hist = False,color="b", kde_kws={"shade": True})
# sns.distplot(LCSdata4_trans, hist = False,color="m", kde_kws={"shade": True})
# sns.distplot(LCSdata5_trans, hist = False,color="c", kde_kws={"shade": True})
# sns.distplot(LCSdata6_trans, hist = False,color="y", kde_kws={"shade": True})

## Get estimated probability density function (PDF) using Kernel Density Estimate (KDE) ##
## When split up into the number of bins in x_grid, can get pairwise differences among the same
## regions, of which the absolute value of those differences can be summed.

# Experimental densities
# density_test = gaussian_kde(test)
min,max,steps = -0.015, 0.015, 100
x_grid = np.linspace(min, max, steps)
dx = (max-min) / steps

density1 = gaussian_kde(rate_data_list[0])
print(density1(x_grid))
integrand = density1(x_grid)
self_dist_mean = 0.5 * np.sqrt(2/(len(rate_data_list[0])*np.pi)) * np.sqrt(sum(integrand*dx))
# self_dist_mean1 = 0.5 * np.sqrt(2/(steps*np.pi)) * sum(integrand*dx)
#quad(np.sqrt(integrand * dx), min, max)
print(self_dist_mean)
# print(self_dist_mean1)
# quit()
density2 = gaussian_kde(rate_data_list[1])
density3 = gaussian_kde(rate_data_list[2])
density4 = gaussian_kde(rate_data_list[3])
density5 = gaussian_kde(rate_data_list[4])
density6 = gaussian_kde(rate_data_list[5])
#
# d_test1 = gaussian_kde(test1)
# d_test2 = gaussian_kde(test2)
# d_test3 = gaussian_kde(test3)
#

## Model densities
# density0 = gaussian_kde(LCSdata0_trans)
# density1 = gaussian_kde(LCSdata1_trans)
# density2 = gaussian_kde(LCSdata2_trans)
# density3 = gaussian_kde(LCSdata3_trans)
# density4 = gaussian_kde(LCSdata4_trans)
# density5 = gaussian_kde(LCSdata5_trans)
# density6 = gaussian_kde(LCSdata6_trans)

## Histogram distances for each distribution ##
# histogram_distance(control_hist = d_test1, test_hist = d_test2, name = "Self Distance")
# histogram_distance(control_hist = d_test1, test_hist = d_test3, name = "Histogram Distance")

hds_allsteps = []
steps = [5, 10,100,1000,10000]
for step in steps:
    hd_p07_p11 = histogram_distance(control_hist = density1, test_hist = density2, name = "Self-Distance",
                                    min = -0.015, max = 0.015, steps = step)
    hd_p07_p15 = histogram_distance(control_hist = density1, test_hist = density3, name = "Distance: P11 to P15",
                                    min = -0.015, max = 0.015, steps = step)
    hd_p07_p19 = histogram_distance(control_hist = density1, test_hist = density4, name = "Distance: P11 to P19",
                                    min = -0.015, max = 0.015, steps = step)
    hd_p07_p28 = histogram_distance(control_hist = density1, test_hist = density5, name = "Distance: P11 to P28",
                                    min = -0.015, max = 0.015, steps = step)
    hd_p07_par = histogram_distance(control_hist = density1, test_hist = density6, name = "Distance: P11 to Parental",
                                    min = -0.015, max = 0.015, steps = step)

    hds = [hd_p07_p11, hd_p07_p15, hd_p07_p19, hd_p07_p28, hd_p07_par]
    hds_allsteps.append(hds)

comparisons = [1,2,3]
comp_ticks = ['p07_p11', 'p07_p19', 'p07_par']
labels = ['5', '10', '100', '1000', '10000']
plt.xticks(comparisons, comp_ticks)
# plt.plot(comparisons, hds_allsteps[0], 'r-o', lw = 3, label = '5')
# plt.plot(comparisons, hds_allsteps[1], 'g-o', lw = 3, label = '10')
getVar = lambda searchList, ind: [searchList[i] for i in ind]
print(hds_allsteps[2])
hd_p07_p11p19par = getVar(hds_allsteps[2], [0,2,4])
plt.plot(comparisons, hd_p07_p11p19par, 'b-o', lw = 5, label = '100')
# plt.plot(comparisons, hds_allsteps[3], 'm-o', lw = 3, label = '1000')
# plt.plot(comparisons, hds_allsteps[4], 'k-o', lw = 3, label = '10000')
plt.plot([0.75,3.25], [self_dist_mean, self_dist_mean], "k--", lw = 5)
plt.text(4,self_dist_mean+0.005, 'Self-Distance')
plt.xlabel('Comparisons', weight = "bold")
plt.ylabel('Histogram Distance', weight = "bold")
plt.title("Histogram Distance Comparisons", weight = "bold")
plt.xlim(0.75,3.25)
# plt.legend(title = 'Delta', loc = 'upper left')
plt.show()

# histogram_distance(control_hist = density1, test_hist = density6, name = "Distance: P11 to test")

# plt.legend()
# plt.show()
quit()

# print(gaussian_kde.integrate_box_1d(density1, -0.1, 0.1))
# print(gaussian_kde.integrate_box_1d(density2, -0.1, 0.1))
# print(gaussian_kde.integrate_box_1d(density3, -0.1, 0.1))

# pdf1 = kde_scipy(cFPdata1_trans, x_grid = x_grid)
# pdf2 = kde_scipy(cFPdata4_trans, x_grid = x_grid)
# pdf3 = kde_scipy(cFPdata6_trans, x_grid = x_grid)

# plt.plot(x_grid, pdf1, color = "red", alpha = 0.4, lw = 3)
# plt.plot(x_grid, pdf2, color = "red", alpha = 0.4, lw = 3)
# plt.plot(x_grid, pdf3, color = "red", alpha = 0.4, lw = 3)

# print('Bootstrap Self Distance', 0.5 * dx * sum(abs(density0(x_grid) - density1(x_grid))))
# print('Distance 1-2', 0.5 * dx * sum(abs(density2(x_grid) - density1(x_grid))))
# print('Distance 1-3',0.5 * dx * sum(abs(density3(x_grid) - density1(x_grid))))
# print('Distance 1-4',0.5 * dx * sum(abs(density4(x_grid) - density1(x_grid))))
# print('Distance 1-5',0.5 * dx * sum(abs(density5(x_grid) - density1(x_grid))))
# print(0.5 * dx * sum(abs(density6(x_grid) - density1(x_grid))))
# print('Distance 1-test',0.5 * dx * sum(abs(density_test(x_grid) - density1(x_grid))))

### Idea: Bootstrap resample data (instead of simulating new dataset)
### Or need to estimate mean and variance of self distance, but that seems tough because not necessarily normal

# def kde_scipy(x, x_grid, bandwidth=0.2, **kwargs):
#     """Kernel Density Estimation with Scipy"""
#     # Note that scipy weights its bandwidth by the covariance of the
#     # input data.  To make the results comparable to the other methods,
#     # we divide the bandwidth by the sample standard deviation here.
#     kde = gaussian_kde(x, bw_method=bandwidth / x.std(ddof=1), **kwargs)
#     return kde.evaluate(x_grid)
#

density1.covariance_factor = lambda : .25
density1._compute_covariance()
density2.covariance_factor = lambda : .25
density2._compute_covariance()
density3.covariance_factor = lambda : .25
density3._compute_covariance()
plt.plot(x_grid,density1(x_grid), color = "red")
plt.plot(x_grid,density2(x_grid), color = "blue")
plt.plot(x_grid,density3(x_grid), color = "green")


nbins = 100
n, bins, _ = plt.hist(cFPdata1_trans, nbins, color = "red")
print(n)
bin_width = bins[1]-bins[0]
integral = bin_width * sum(n)
print(integral)

p, bins, _ = plt.hist(cFPdata4_trans, nbins, color = "green")
print(p)
bin_width = bins[1]-bins[0]
integral1 = bin_width * sum(p)
print(integral1)

m, bins1, _ = plt.hist(cFPdata6_trans, nbins, color = "blue")
print(m)
bin_width1 = bins1[1]-bins1[0]
integral2 = bin_width1 * sum(m)
print(integral2)

print(p-n)
print(m-n)

print(abs(p-n))
print(abs(m-n))

print(bin_width * abs(p-n))
print(bin_width * abs(m-n))

integral_n = bin_width * sum(n)
integral_p = bin_width * sum(p)
integral_m = bin_width * sum(m)

print(integral_n, integral_p, integral_m)

# hd = 0.5
plt.show()
quit()
print(abs(p-n))
print(abs(m-n))

quit()
hist_dist_np = bin_width * 0.5 * sum(abs(p[0:99]-n[0:99]))
hist_dist_nm = bin_width * 0.5 * sum(abs(m[0:99]-n[0:99]))

print(hist_dist_np, hist_dist_nm)
quit()
integrand1 = bin_width * abs(sum(p)-sum(n))

hist_dist = 0.5 * sp.integrate(bin_width * abs(p-n))
print(hist_dist)
plt.show()
quit()
hist_dist = 0.5 * Integral(abs())
a = integral1-integral
b = integral2-integral1

hist_dist1 = 0.5 * a
hist_dist2 = 0.5 * b
print(hist_dist1, hist_dist2)
plt.show()
quit()

# m,n = sns.distplot(cFPdata1_trans).get_lines()[0].get_data()
# print(m)
# print(n)

# support = np.linspace(min(m), max(m), len(m))
# # d1_h = [h.get_height() for h in sns.distplot(cFPdata1_trans).patches]
x = Symbol('x')
i = Symbol('i')
delta = 0.00004
N = 1000

print(cFPdata1_trans)
print(min(cFPdata1_trans), max(cFPdata1_trans))
a = Integral(f(x) * DiracDelta(x-x0), (x,a,b))
print(a)
quit()
print(DiracDelta(x))
print integrate(DiracDelta(x), (x, 0, 5.0))
print integrate(DiracDelta(x), (x, -1, 1))
print integrate(DiracDelta(x), (x, 0, 0 + delta))

exp = Sum(integrate(DiracDelta(x), (x, 0, 0 + delta)), (i, 1, N))
print(exp)
quit()

# print((1/(N*delta)) * )

quit()
N = 1000
delta = 0.00004
x = cFPdata1_trans
print(x)
print(DiracDelta(x))
quit()
integrand = DiracDelta(x) * delta
h = 1./(N*delta) * (np.sum(sp.integrate(integrand, x, x, x+delta)))
print(h)


