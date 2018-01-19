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
from scipy.integrate import quad

sns.set(font_scale = 3)
sns.set_style("white")

expData_p07 = np.genfromtxt("SKMEL5_parentalP22_slopesCI.txt", delimiter='\t')
expData_p07_df = pd.DataFrame(expData_p07)
expData_p07_df.columns = ['Replicate', 'DIP Rates', 'CI_lower', 'CI_upper']
p07rates = expData_p07_df[['DIP Rates']]
p07rates_noNull = p07rates.apply(lambda x: [y for y in x if pd.notnull(y)])
p07rate_values = p07rates_noNull.values[0]
p07rate_vals = np.asarray(p07rate_values)

# SLOPE DATA with CIs

expData_cFP1 = np.genfromtxt("SKMEL5par_H2BGFP_cFP_slopes.txt", delimiter='\t')
expData_cFP2 = np.genfromtxt("SKMEL5par_H2BGFP_cFP2_slopes.txt", delimiter='\t')

expData_cFPall = np.concatenate((expData_cFP1, expData_cFP2))

labels = ["cFP1", "cFP2"]

# plt.figure()
# sns.distplot(expData_cFP1, kde = False, color="r", bins = 50)
# sns.distplot(expData_cFP2, kde = False, color="g", bins = 50)
# plt.xlim(-0.1,0.05)
# plt.xlabel("DIP Rate", weight = "bold")
# plt.ylabel("Density", weight = "bold")
# plt.title("SKMEL5 H2B-GFP cFP Comparison - 100 cell number", weight = "bold")
# plt.legend(title = "Experiments", labels = labels, loc = "upper left")
#
# plt.figure()
# sns.distplot(expData_cFPall, kde = False, color="b", bins = 50)
# plt.xlim(-0.1,0.05)
# plt.xlabel("DIP Rate", weight = "bold")
# plt.ylabel("Density", weight = "bold")
# plt.title("SKMEL5 H2B-GFP cFPs Merged - 100 cell number", weight = "bold")
#
#
# plt.figure()
# sns.distplot(expData_cFP1, hist = False, color="r", kde_kws={"shade": True})
# sns.distplot(expData_cFP2, hist = False, color="g", kde_kws={"shade": True})
# plt.xlim(-0.1,0.05)
# plt.xlabel("DIP Rate", weight = "bold")
# plt.ylabel("Density", weight = "bold")
# plt.title("SKMEL5 H2B-GFP cFP Comparison - 150 cell number", weight = "bold")
# plt.legend(title = "Experiments", labels = labels, loc = "upper left")

labels1 = ["cFP", "LSD"]
plt.figure()
sns.distplot(expData_cFPall, hist = False, color="b", kde_kws={"shade": True})
sns.distplot(p07rate_vals, hist = False, color = "m", kde_kws={"shade": True})
# plt.xlim(-0.1,0.05)
plt.xlabel("DIP Rate", weight = "bold")
plt.ylabel("Density", weight = "bold")
plt.title("SKMEL5 H2B-GFP cFP Merged vs LSD", weight = "bold")
plt.legend(title = "Assay", labels = labels1, loc = "upper left")
plt.show()