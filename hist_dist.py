# Do numerical integral, delta x small enough, run loop and get distance between curves and abs value it
import numpy as np
# import scipy.integrate as integ
from sympy import DiracDelta
from sympy import *
from sympy.mpmath import *
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp
import pandas as pd
from numpy import median
from scipy.stats import gaussian_kde

# Load model produced data
# LCSdata0_trans = np.load("1.0_0.0_0.0_10kData_LCS_forselfdistance.npy")
# LCSdata1_trans = np.load("1.0_0.0_0.0_10kData.npy")
# LCSdata2_trans = np.load("0.95_0.05_0.0_10kData.npy")
# LCSdata3_trans = np.load("0.85_0.1_0.05_10kData.npy")
# LCSdata4_trans = np.load("0.7_0.2_0.1_10kData.npy")
# LCSdata5_trans = np.load("0.5_0.3_0.2_10kData.npy")
# LCSdata6_trans = np.load("0.25_0.5_0.25_10kData.npy")
#
expData_p11 = np.loadtxt("experimental_psg11.2.txt", delimiter='\t')
expData_p15 = np.loadtxt("experimental_psg15.2.txt", delimiter='\t')
expData_p19 = np.loadtxt("experimental_psg19.2.txt", delimiter='\t')
expData_p28 = np.loadtxt("experimental_psg28.2.txt", delimiter='\t')
expData_par = np.loadtxt("experimental_parental2.txt", delimiter='\t')


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
    """
    if n == None:
        n = len(X)

    resample_i = np.floor(np.random.rand(n) * len(X)).astype(int)
    X_resample = X[resample_i]
    return X_resample

# Random different distribution to make sure distance is ~1
test = np.random.normal(0.07, 0.001, 10000)

expData_p11_resample = bootstrap_resample(expData_p11)
# expData_p15_resample = bootstrap_resample(expData_p15)
# expData_p19_resample = bootstrap_resample(expData_p19)
# expData_p28_resample = bootstrap_resample(expData_p28)
# expData_par_resample = bootstrap_resample(expData_par)

# cFPdata0_trans = np.load("1.0_0.0_0.0_10kData_cFP_forselfdistance.npy")
# cFPdata1_trans = np.load("1.0_0.0_0.0_10kData_cFP.npy")
# cFPdata2_trans = np.load("0.95_0.05_0.0_10kData_cFP.npy")
# cFPdata3_trans = np.load("0.85_0.1_0.05_10kData_cFP.npy")
# cFPdata4_trans = np.load("0.7_0.2_0.1_10kData_cFP.npy")
# cFPdata5_trans = np.load("0.5_0.3_0.2_10kData_cFP.npy")
# cFPdata6_trans = np.load("0.25_0.5_0.25_10kData_cFP.npy")

# def kde_scipy(x, x_grid, bandwidth=0.2, **kwargs):
#     """Kernel Density Estimation with Scipy"""
#     # Note that scipy weights its bandwidth by the covariance of the
#     # input data.  To make the results comparable to the other methods,
#     # we divide the bandwidth by the sample standard deviation here.
#     kde = gaussian_kde(x, bw_method=bandwidth / x.std(ddof=1), **kwargs)
#     return kde.evaluate(x_grid)
#
min = -0.1
max = 0.1
steps = 1000

x_grid = np.linspace(min, max, steps)
# pdf1 = kde_scipy(cFPdata1_trans, x_grid = x_grid)
# pdf2 = kde_scipy(cFPdata4_trans, x_grid = x_grid)
# pdf3 = kde_scipy(cFPdata6_trans, x_grid = x_grid)
#
# plt.plot(x_grid, pdf1, color = "red", alpha = 0.4, lw = 3)
# plt.plot(x_grid, pdf2, color = "red", alpha = 0.4, lw = 3)
# plt.plot(x_grid, pdf3, color = "red", alpha = 0.4, lw = 3)
#
sns.distplot(test, hist = False,color="c", kde_kws={"shade": True}, label = "Different")
sns.distplot(expData_p11, hist = False,color="k", kde_kws={"shade": True}, label = "Passage 11")
sns.distplot(expData_p11_resample, hist = False,color="y", kde_kws={"shade": True}, label = "Passage 11 Resample")
sns.distplot(expData_p15, hist = False,color="r", kde_kws={"shade": True}, label = "Passage 15")
sns.distplot(expData_p19, hist = False,color="g", kde_kws={"shade": True}, label = "Passage 19")
sns.distplot(expData_p28, hist = False,color="b", kde_kws={"shade": True}, label = "Passage 28")
sns.distplot(expData_par, hist = False,color="m", kde_kws={"shade": True}, label = "Parental")


# sns.distplot(LCSdata0_trans, hist = False,color="k", kde_kws={"shade": True})
# sns.distplot(LCSdata1_trans, hist = False,color="r", kde_kws={"shade": True})
# sns.distplot(LCSdata2_trans, hist = False,color="g", kde_kws={"shade": True})
# sns.distplot(LCSdata3_trans, hist = False,color="b", kde_kws={"shade": True})
# sns.distplot(LCSdata4_trans, hist = False,color="m", kde_kws={"shade": True})
# sns.distplot(LCSdata5_trans, hist = False,color="c", kde_kws={"shade": True})
# sns.distplot(LCSdata6_trans, hist = False,color="y", kde_kws={"shade": True})

plt.legend()

density_test = gaussian_kde(test)
density0 = gaussian_kde(expData_p11_resample)
density1 = gaussian_kde(expData_p11)
density2 = gaussian_kde(expData_p15)
density3 = gaussian_kde(expData_p19)
density4 = gaussian_kde(expData_p28)
density5 = gaussian_kde(expData_par)

# density0 = gaussian_kde(LCSdata0_trans)
# density1 = gaussian_kde(LCSdata1_trans)
# density2 = gaussian_kde(LCSdata2_trans)
# density3 = gaussian_kde(LCSdata3_trans)
# density4 = gaussian_kde(LCSdata4_trans)
# density5 = gaussian_kde(LCSdata5_trans)
# density6 = gaussian_kde(LCSdata6_trans)

dx = (max-min) / steps
print('Bootstrap Self Distance', 0.5 * dx * sum(abs(density0(x_grid) - density1(x_grid))))
print('Distance 1-2', 0.5 * dx * sum(abs(density2(x_grid) - density1(x_grid))))
print('Distance 1-3',0.5 * dx * sum(abs(density3(x_grid) - density1(x_grid))))
print('Distance 1-4',0.5 * dx * sum(abs(density4(x_grid) - density1(x_grid))))
print('Distance 1-5',0.5 * dx * sum(abs(density5(x_grid) - density1(x_grid))))
# print(0.5 * dx * sum(abs(density6(x_grid) - density1(x_grid))))
print('Distance 1-test',0.5 * dx * sum(abs(density_test(x_grid) - density1(x_grid))))

plt.show()
# print(gaussian_kde.integrate_box_1d(density1, -0.1, 0.1))
# print(gaussian_kde.integrate_box_1d(density2, -0.1, 0.1))
# print(gaussian_kde.integrate_box_1d(density3, -0.1, 0.1))


### Idea: Bootstrap resample data (instead of simulating new dataset)
### Or need to estimate mean and variance of self distance, but that seems tough because not necessarily normal
quit()

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
#
# LCSdata1_notrans = np.load("1.0_0.0_0.0_10kData_LCS_noTrans.npy")
# LCSdata2_notrans = np.load("0.95_0.05_0.0_10kData_LCS_noTrans.npy")
# LCSdata3_notrans = np.load("0.85_0.1_0.05_10kData_LCS_noTrans.npy")
# LCSdata4_notrans = np.load("0.7_0.2_0.1_10kData_LCS_noTrans.npy")
# LCSdata5_notrans = np.load("0.5_0.3_0.2_10kData_LCS_noTrans.npy")
# LCSdata6_notrans = np.load("0.25_0.5_0.25_10kData_LCS_noTrans.npy")
#
# cFPdata1_notrans = np.load("1.0_0.0_0.0_10kData_cFP_noTrans.npy")
# cFPdata2_notrans = np.load("0.95_0.05_0.0_10kData_cFP_noTrans.npy")
# cFPdata3_notrans = np.load("0.85_0.1_0.05_10kData_cFP_noTrans.npy")
# cFPdata4_notrans = np.load("0.7_0.2_0.1_10kData_cFP_noTrans.npy")
# cFPdata5_notrans = np.load("0.5_0.3_0.2_10kData_cFP_noTrans.npy")
# cFPdata6_notrans = np.load("0.25_0.5_0.25_10kData_cFP_noTrans.npy")
#
# expData_p11 = np.loadtxt("experimental_psg11.2.txt", delimiter='\t')
# expData_p15 = np.loadtxt("experimental_psg15.2.txt", delimiter='\t')
# expData_p19 = np.loadtxt("experimental_psg19.2.txt", delimiter='\t')
# expData_p28 = np.loadtxt("experimental_psg28.2.txt", delimiter='\t')
# expData_par = np.loadtxt("experimental_parental2.txt", delimiter='\t')

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


