import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp
import pandas as pd
import re
from numpy import median
from scipy.stats import gaussian_kde

# # Load model produced data
# LCSdata0_trans = np.load("1.0_0.0_0.0_10kData_LCS_forselfdistance.npy")
# LCSdata1_trans = np.load("1.0_0.0_0.0_10kData.npy")
# LCSdata2_trans = np.load("0.95_0.05_0.0_10kData.npy")
# LCSdata3_trans = np.load("0.85_0.1_0.05_10kData.npy")
# LCSdata4_trans = np.load("0.7_0.2_0.1_10kData.npy")
# LCSdata5_trans = np.load("0.5_0.3_0.2_10kData.npy")
# LCSdata6_trans = np.load("0.25_0.5_0.25_10kData.npy")
#
# cFPdata0_trans = np.load("1.0_0.0_0.0_10kData_cFP_forselfdistance.npy")
# cFPdata1_trans = np.load("1.0_0.0_0.0_10kData_cFP.npy")
# cFPdata2_trans = np.load("0.95_0.05_0.0_10kData_cFP.npy")
# cFPdata3_trans = np.load("0.85_0.1_0.05_10kData_cFP.npy")
# cFPdata4_trans = np.load("0.7_0.2_0.1_10kData_cFP.npy")
# cFPdata5_trans = np.load("0.5_0.3_0.2_10kData_cFP.npy")
# cFPdata6_trans = np.load("0.25_0.5_0.25_10kData_cFP.npy")
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

# SLOPE DATA with CIs

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

data = [expData_p07_df, expData_p11_df, expData_p15_df,
        expData_p19_df, expData_p28_df, expData_par_df]

DIPs = []
for dat in data:
    dat_DIP = dat[['DIP Rates']]
    DIPs.append(dat_DIP)
print(DIPs)
DIPs_cat = pd.concat(DIPs, axis = 1)
DIPs_cat.columns = ["Passage 07", "Passage 11", "Passage 15", "Passage 19",
           "Passage 28", "Parental"]
print(DIPs_cat)
melt_DIPs = pd.melt(DIPs_cat)
print(melt_DIPs)


sns.set(font_scale = 2)
sns.set_style("whitegrid")
#
# labels1 = ["Passage 07", "Passage 11", "Passage 15", "Passage 19",
#            "Passage 28", "Parental"]
# sns.pointplot(x = "variable", y = "value", data = melt_DIPs, estimator=median)
# sns.plt.title("Comparison Across Passages", weight = "bold")
# plt.ylabel("DIP Rate", weight = "bold")
# plt.xlabel("Population", weight = "bold")
# # plt.xlim(-0.005,0.005)
# # plt.ylim(0,500)
# # plt.legend(ncol = 3, labels = labels1, loc = 0)
# plt.show()
# quit()

# print(expData_p07_df[['DIP Rates']])



sns.distplot(expData_p07_df[['DIP Rates']], hist = False, color="y", kde_kws={"shade": True})
sns.distplot(expData_p11_df[['DIP Rates']], hist = False, color="r", kde_kws={"shade": True})
sns.distplot(expData_p15_df[['DIP Rates']], hist = False, color="g", kde_kws={"shade": True})
sns.distplot(expData_p19_df[['DIP Rates']], hist = False, color="b", kde_kws={"shade": True})
sns.distplot(expData_p28_df[['DIP Rates']], hist = False, color="m", kde_kws={"shade": True})
sns.distplot(expData_par_df[['DIP Rates']], hist = False, color="k", kde_kws={"shade": True})

sns.plt.title("DIP Rates from LSD", weight = "bold")
labels1 = ["Passage 07", "Passage 11", "Passage 15", "Passage 19",
           "Passage 28", "Parental"]
plt.ylabel("Density", weight = "bold")
plt.xlabel("DIP Rate", weight = "bold")
plt.xlim(-0.02,0.02)
plt.ylim(0,100)
plt.legend(labels = labels1)# ax.legend(labels = labels, fontsize = 16)


plt.show()
quit()
# expData_p11 = np.loadtxt("experimental_psg11.2.txt", delimiter='\t')
# expData_p15 = np.loadtxt("experimental_psg15.2.txt", delimiter='\t')
# expData_p19 = np.loadtxt("experimental_psg19.2.txt", delimiter='\t')
# expData_p28 = np.loadtxt("experimental_psg28.2.txt", delimiter='\t')
# expData_par = np.loadtxt("experimental_parental2.txt", delimiter='\t')
#

frames_LCStrans = [LCSdata1_trans,LCSdata2_trans,LCSdata3_trans,LCSdata4_trans,LCSdata5_trans,LCSdata6_trans]
frames_LCStrans_mod = [LCSdata1_trans,LCSdata2_trans,LCSdata3_trans,LCSdata4_trans, LCSdata6_trans]
frames_cFPtrans = [cFPdata1_trans,cFPdata2_trans,cFPdata3_trans,cFPdata4_trans,cFPdata5_trans,cFPdata6_trans]
frames_cFPtrans_mod = [cFPdata1_trans,cFPdata2_trans,cFPdata3_trans,cFPdata4_trans,cFPdata6_trans]
frames_LCSnotrans = [LCSdata1_notrans,LCSdata2_notrans,LCSdata3_notrans,LCSdata4_notrans,LCSdata5_notrans,LCSdata6_notrans]
frames_cFPnotrans = [cFPdata1_notrans,cFPdata2_notrans,cFPdata3_notrans,cFPdata4_notrans,cFPdata5_notrans,cFPdata6_notrans]
expData = [expData_p11, expData_p15, expData_p19, expData_p28, expData_par]

# labels = ['Passage 11', 'Passage 28']



## Setup grid over which to "integrate" distance difference
min = -0.1
max = 0.1
steps = 100
x_grid = np.linspace(min, max, steps)
# Multiplication factor to ensure values are from 0-1 --> completes "integral"
dx = (max-min) / steps

def histogram_distance(control_hist, test_hist, name, dx = dx, x_grid = x_grid):
    dist = 0.5 * dx * sum(abs(test_hist(x_grid) - control_hist(x_grid)))
    return [name, dist]

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

# # fig, ax = plt.subplots()
# # fig,axs = plt.subplots(nrows=2)
# # plt.figure("cFP Plot")
# # m,n = sns.distplot(cFPdata1_trans).get_lines()[0].get_data()
# # d1_h = [h.get_height() for h in sns.distplot(cFPdata1_trans).patches]
# # print(len(d1_h))
# # print(n)
#
# # ax = None
# a = sns.distplot(cFPdata5_trans).get_lines()[0].get_data()
# # d4_h = [h.get_height() for h in sns.distplot(cFPdata4_trans).patches]
# # print(len(d4_h))
# np.save("cFPdata5.npy", a)
# # np.save("cFPdata1_y.npy")
# quit()
#
# print(p)
# # print(d4_h)
# quit()
# plt.plot(m,n)
# plt.plot(o,p)
# plt.show()
# print(p-n)
# quit()
#
# p = sns.distplot(cFPdata1_trans, hist = False,color="b", kde_kws={"shade": True})
# x,y = p.get_lines()[0].get_data()
# p_cdf = sp.integrate.cumtrapz(y, x, initial=0)
# print(p_cdf)
# quit()

# plt.plot(x,y)
# print(x)
# print(len(x))
# print(y)
# print(len(y))

# q = sns.distplot(cFPdata4_trans, hist = False,color="r", kde_kws={"shade": True})
# a,b = q.get_lines()[0].get_data()
# q_cdf = sp.integrate.cumtrapz(b, a, initial=0)
# print(q_cdf)
# print(q_cdf-p_cdf)
# quit()
# plt.plot(a,b)
# plt.show()
#
# # print(a)
# # print(len(a))
# print(b)
# # print(len(b))
#
# print(y-b)
# quit()
#

expData_p11_resample = bootstrap_resample(expData_p11)

cFPdata1_subsample = np.random.choice(cFPdata1_notrans, 100, replace=False)
cFPdata2_subsample = np.random.choice(cFPdata2_notrans, 100, replace=False)
cFPdata3_subsample = np.random.choice(cFPdata3_notrans, 100, replace=False)
cFPdata4_subsample = np.random.choice(cFPdata4_notrans, 100, replace=False)
cFPdata5_subsample = np.random.choice(cFPdata5_notrans, 100, replace=False)
cFPdata6_subsample = np.random.choice(cFPdata6_notrans, 100, replace=False)

LCSdata1_subsample = np.random.choice(LCSdata1_notrans, 350, replace=False)
LCSdata2_subsample = np.random.choice(LCSdata2_notrans, 350, replace=False)
LCSdata3_subsample = np.random.choice(LCSdata3_notrans, 350, replace=False)
LCSdata4_subsample = np.random.choice(LCSdata4_notrans, 350, replace=False)
LCSdata5_subsample = np.random.choice(LCSdata5_notrans, 350, replace=False)
LCSdata6_subsample = np.random.choice(LCSdata6_notrans, 350, replace=False)

data_nodf = [cFPdata1_notrans, cFPdata2_notrans, cFPdata3_notrans,
        cFPdata4_notrans, cFPdata5_notrans, cFPdata6_notrans,
        cFPdata1_subsample, cFPdata2_subsample, cFPdata3_subsample,
        cFPdata4_subsample, cFPdata5_subsample, cFPdata6_subsample,
        LCSdata1_subsample, LCSdata2_subsample, LCSdata3_subsample,
        LCSdata4_subsample, LCSdata5_subsample, LCSdata6_subsample]

data = [pd.DataFrame(cFPdata1_notrans), pd.DataFrame(cFPdata2_notrans), pd.DataFrame(cFPdata3_notrans),
        pd.DataFrame(cFPdata4_notrans), pd.DataFrame(cFPdata5_notrans), pd.DataFrame(cFPdata6_notrans),
        pd.DataFrame(cFPdata1_subsample), pd.DataFrame(cFPdata2_subsample), pd.DataFrame(cFPdata3_subsample),
        pd.DataFrame(cFPdata4_subsample), pd.DataFrame(cFPdata5_subsample), pd.DataFrame(cFPdata6_subsample),
        pd.DataFrame(LCSdata1_subsample), pd.DataFrame(LCSdata2_subsample), pd.DataFrame(LCSdata3_subsample),
        pd.DataFrame(LCSdata4_subsample), pd.DataFrame(LCSdata5_subsample), pd.DataFrame(LCSdata6_subsample)]

df = pd.concat(data, ignore_index=True, axis = 1)

kdes = []
for kde in data_nodf:
    kdes.append(gaussian_kde(kde))
# print(len(kdes))

hds = []
hds.append([histogram_distance(control_hist=kdes[0], test_hist=kdes[6], name='PDF_cFP_HD1'),
            histogram_distance(control_hist=kdes[0], test_hist=kdes[12], name='PDF_LSD_HD1'),
            histogram_distance(control_hist=kdes[6], test_hist=kdes[12], name='cFP_LSD_HD1')])

hds.append([histogram_distance(control_hist=kdes[1], test_hist=kdes[7], name='PDF_cFP_HD2'),
            histogram_distance(control_hist=kdes[1], test_hist=kdes[13], name='PDF_LSD_HD2'),
            histogram_distance(control_hist=kdes[7], test_hist=kdes[13], name='cFP_LSD_HD2')])

hds.append([histogram_distance(control_hist=kdes[2], test_hist=kdes[8], name='PDF_cFP_HD3'),
            histogram_distance(control_hist=kdes[2], test_hist=kdes[14], name='PDF_LSD_HD3'),
            histogram_distance(control_hist=kdes[8], test_hist=kdes[14], name='cFP_LSD_HD3')])

hds.append([histogram_distance(control_hist=kdes[3], test_hist=kdes[9], name='PDF_cFP_HD4'),
            histogram_distance(control_hist=kdes[3], test_hist=kdes[15], name='PDF_LSD_HD4'),
            histogram_distance(control_hist=kdes[9], test_hist=kdes[15], name='cFP_LSD_HD4')])

hds.append([histogram_distance(control_hist=kdes[4], test_hist=kdes[10], name='PDF_cFP_HD5'),
           histogram_distance(control_hist=kdes[4], test_hist=kdes[16], name='PDF_LSD_HD5'),
           histogram_distance(control_hist=kdes[10], test_hist=kdes[16], name='cFP_LSD_HD5')])

hds.append([histogram_distance(control_hist=kdes[5], test_hist=kdes[11], name='PDF_cFP_HD6'),
           histogram_distance(control_hist=kdes[5], test_hist=kdes[17], name='PDF_LSD_HD6'),
           histogram_distance(control_hist=kdes[11], test_hist=kdes[17], name='cFP_LSD_HD6')])

print pd.DataFrame(hds)
print(hds[0])
print(hds[0][0])
print(hds[0][0][0])
print(hds[0][0][1])

#
df.columns = ['PDF_data1', 'PDF_data2', 'PDF_data3', 'PDF_data4', 'PDF_data5', 'PDF_data6',
             'cFP_data1', 'cFP_data2', 'cFP_data3',
             'cFP_data4', 'cFP_data5', 'cFP_data6',
             'LCS_data1', 'LCS_data2', 'LCS_data3',
             'LCS_data4', 'LCS_data5', 'LCS_data6']
print(df)

lf_df = pd.melt(df)
print(lf_df)
lf_df.columns = ['variable', 'DIP Rate']
lf_df['plot_type'] = lf_df['variable'].str.split('_').str[0]
lf_df['data_type'] = lf_df['variable'].str.split('_').str[1]
del lf_df['variable']
print(lf_df)


labels2 = ["Approximate PDF", "cFP", "LSD"]
fg = sns.FacetGrid(lf_df, col = "data_type", hue = "plot_type", col_wrap=2)
fg = (fg.map(sns.distplot, "DIP Rate", hist = False, kde_kws={"shade": True}))

fg.set_axis_labels("DIP Rate", "Density")
fg.set_titles("{col_name}", weight = "bold")
#

axes_labs = [0,1,2,3,4,5]
col_order = ['Low', 'Medium-Low', 'Medium', 'Medium-High', 'High', 'Equilibrium']

for ax, title in zip(fg.axes.flat, col_order):
    for nums in axes_labs:
        ax.set_title(title, weight = 'bold')
        ax.text(0.85, 0.85,
                "Hello", fontsize = 20)
                # '%s = %1.2f' % (hds[nums][0][0], hds[nums][0][1]),
                # '%s = %1.2f' % (hds[nums][1][0], hds[nums][1][1]),
                # '%s = %1.2f' % (hds[nums][2][0], hds[nums][2][1])) #add text
# fg.add_legend()
# sns.plt.title("Distribution Comparison", weight = "bold")
# plt.ylabel("Density", weight = "bold")
# plt.xlabel("DIP Rate", weight = "bold")
# # plt.xlim(-0.005,0.005)
# # plt.ylim(0,500)
plt.legend(ncol = 3, labels = labels2, loc = 'lower center', bbox_to_anchor=(-0.35,-0.55), borderaxespad=0.)

# sns.plt.legend(loc = 'bottom', bbox_to_anchor=  labels = labels2)
plt.show()
quit()



plt.figure("Comparing Distributions - No Diversity")
sns.distplot(cFPdata1_notrans, hist = False,color="k", kde_kws={"shade": True})
sns.distplot(cFPdata1_subsample, hist = False,color="r", kde_kws={"shade": True})
sns.distplot(LCSdata1_subsample, hist = False,color="b", kde_kws={"shade": True})
sns.plt.title("Distribution Comparison", weight = "bold")
plt.ylabel("Density", weight = "bold")
plt.xlabel("DIP Rate", weight = "bold")
# plt.xlim(-0.005,0.005)
# plt.ylim(0,500)
plt.legend(labels = labels2)

plt.figure("Comparing Distributions - Low Diversity")
sns.distplot(cFPdata2_notrans, hist = False,color="k", kde_kws={"shade": True})
sns.distplot(cFPdata2_subsample, hist = False,color="r", kde_kws={"shade": True})
sns.distplot(LCSdata2_subsample, hist = False,color="b", kde_kws={"shade": True})
sns.plt.title("Distribution Comparison", weight = "bold")
plt.ylabel("Density", weight = "bold")
plt.xlabel("DIP Rate", weight = "bold")
# plt.xlim(-0.005,0.005)
# plt.ylim(0,500)
plt.legend(labels = labels2)

plt.figure("Comparing Distributions - Low-Medium Diversity")
sns.distplot(cFPdata2_notrans, hist = False,color="k", kde_kws={"shade": True})
sns.distplot(cFPdata2_subsample, hist = False,color="r", kde_kws={"shade": True})
sns.distplot(LCSdata2_subsample, hist = False,color="b", kde_kws={"shade": True})
sns.plt.title("Distribution Comparison", weight = "bold")
plt.ylabel("Density", weight = "bold")
plt.xlabel("DIP Rate", weight = "bold")
# plt.xlim(-0.005,0.005)
# plt.ylim(0,500)
plt.legend(labels = labels2)

plt.figure("Comparing Distributions - Medium Diversity")
sns.distplot(cFPdata3_notrans, hist = False,color="k", kde_kws={"shade": True})
sns.distplot(cFPdata3_subsample, hist = False,color="r", kde_kws={"shade": True})
sns.distplot(LCSdata3_subsample, hist = False,color="b", kde_kws={"shade": True})
sns.plt.title("Distribution Comparison", weight = "bold")
plt.ylabel("Density", weight = "bold")
plt.xlabel("DIP Rate", weight = "bold")
# plt.xlim(-0.005,0.005)
# plt.ylim(0,500)
plt.legend(labels = labels2)

plt.figure("Comparing Distributions - Medium-High Diversity")
sns.distplot(cFPdata4_notrans, hist = False,color="k", kde_kws={"shade": True})
sns.distplot(cFPdata4_subsample, hist = False,color="r", kde_kws={"shade": True})
sns.distplot(LCSdata4_subsample, hist = False,color="b", kde_kws={"shade": True})
sns.plt.title("Distribution Comparison", weight = "bold")
plt.ylabel("Density", weight = "bold")
plt.xlabel("DIP Rate", weight = "bold")
# plt.xlim(-0.005,0.005)
# plt.ylim(0,500)
plt.legend(labels = labels2)

plt.figure("Comparing Distributions - High Diversity")
sns.distplot(cFPdata5_notrans, hist = False,color="k", kde_kws={"shade": True})
sns.distplot(cFPdata5_subsample, hist = False,color="r", kde_kws={"shade": True})
sns.distplot(LCSdata5_subsample, hist = False,color="b", kde_kws={"shade": True})
sns.plt.title("Distribution Comparison", weight = "bold")
plt.ylabel("Density", weight = "bold")
plt.xlabel("DIP Rate", weight = "bold")
# plt.xlim(-0.005,0.005)
# plt.ylim(0,500)
plt.legend(labels = labels2)

plt.figure("Comparing Distributions - Equilibrated Diversity")
sns.distplot(cFPdata6_notrans, hist = False,color="r", kde_kws={"shade": True})
sns.distplot(cFPdata6_subsample, hist = False,color="g", kde_kws={"shade": True})
sns.distplot(LCSdata6_subsample, hist = False,color="b", kde_kws={"shade": True})

sns.plt.title("Distribution Comparison", weight = "bold")
labels2 = ["Approximate PDF", "cFP", "LSD"]
plt.ylabel("Density", weight = "bold")
plt.xlabel("DIP Rate", weight = "bold")
# plt.xlim(-0.005,0.005)
# plt.ylim(0,500)
plt.legend(labels = labels2)# ax.legend(labels = labels, fontsize = 16)
plt.show()
quit()


plt.figure("LSD Plot - Experimental")
sns.distplot(expData_p11_resample, hist = False,color="k", kde_kws={"shade": True})
sns.distplot(expData_p11, hist = False,color="r", kde_kws={"shade": True})
sns.distplot(expData_p15, hist = False,color="g", kde_kws={"shade": True})
sns.distplot(expData_p19, hist = False,color="b", kde_kws={"shade": True})
sns.distplot(expData_p28, hist = False,color="m", kde_kws={"shade": True})
sns.distplot(expData_par, hist = False,color="y", kde_kws={"shade": True})


sns.plt.title("DIP Rates from LSD", weight = "bold")
labels1 = ["Passage 11 Bootstrap", "Passage 11", "Passage 15", "Passage 19",
           "Passage 28", "Parental"]
plt.ylabel("Density", weight = "bold")
plt.xlabel("DIP Rate", weight = "bold")
plt.xlim(-0.005,0.005)
plt.ylim(0,500)
plt.legend(labels = labels1)# ax.legend(labels = labels, fontsize = 16)
plt.show()
quit()

plt.figure("cFP Plot")
sns.distplot(cFPdata0_trans, hist = False,color="k", kde_kws={"shade": True})
sns.distplot(cFPdata1_trans, hist = False,color="r", kde_kws={"shade": True})
sns.distplot(cFPdata3_trans, hist = False,color="g", kde_kws={"shade": True})
sns.distplot(cFPdata4_trans, hist = False,color="b", kde_kws={"shade": True})
sns.distplot(cFPdata6_trans, hist = False,color="m", kde_kws={"shade": True})


sns.plt.title("DIP Rates from cFP", weight = "bold")
labels = ["Low Resample", "Low", "Medium-Low", "Medium-High", "High"]
plt.ylabel("Density", weight = "bold")
plt.xlabel("DIP Rate", weight = "bold")
plt.xlim(-0.02,0.02)
plt.ylim(0,225)
plt.legend(labels = labels)# ax.legend(labels = labels, fontsize = 16)

plt.figure("LSD Plot")
sns.distplot(LCSdata0_trans, hist = False,color="k", kde_kws={"shade": True})
sns.distplot(LCSdata1_trans, hist = False,color="r", kde_kws={"shade": True})
sns.distplot(LCSdata3_trans, hist = False,color="g", kde_kws={"shade": True})
sns.distplot(LCSdata4_trans, hist = False,color="b", kde_kws={"shade": True})
sns.distplot(LCSdata6_trans, hist = False,color="m", kde_kws={"shade": True})

sns.plt.title("DIP Rates from LSD", weight = "bold")
plt.ylabel("Density", weight = "bold")
plt.xlabel("DIP Rate", weight = "bold")
plt.xlim(-0.02,0.02)
plt.ylim(0,225)
plt.legend(labels = labels)

plt.show()
quit()

print(np.mean(cFPdata1_trans))
print(np.mean(cFPdata5_trans))

print(np.mean(LCSdata1_trans))
print(np.mean(LCSdata5_trans))

s,p = sp.stats.ttest_ind(cFPdata1_trans, cFPdata2_trans)
print(p)


for exp in expData[1:]:
    s,p = sp.stats.ttest_ind(expData_p11, exp, equal_var=False)
    print(p)
# quit()

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), sp.stats.sem(a)
    h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
    # return m, m-h, m+h
    return h
    # print(m, m-h, m+h)

fig,ax = plt.subplots()
# # plt.plot([np.mean(cFPdata1_trans),np.mean(cFPdata2_trans),np.mean(cFPdata3_trans),np.mean(cFPdata4_trans),
# #           np.mean(cFPdata6_trans)], color = "green")
ax.errorbar(x = [11,15,19,28,35], y = [np.mean(cFPdata1_trans),np.mean(cFPdata2_trans),np.mean(cFPdata3_trans),
                                       np.mean(cFPdata4_trans), np.mean(cFPdata6_trans)],
             yerr = [mean_confidence_interval(cFPdata1_trans), mean_confidence_interval(cFPdata2_trans),
                     mean_confidence_interval(cFPdata3_trans),mean_confidence_interval(cFPdata4_trans),
                     mean_confidence_interval(cFPdata6_trans)], color = "green", lw = 5)

ax.errorbar(x=[11, 15, 19, 28, 35], y=[np.mean(LCSdata1_trans), np.mean(LCSdata2_trans), np.mean(LCSdata3_trans),
                                       np.mean(LCSdata4_trans), np.mean(LCSdata6_trans)],
            yerr=[mean_confidence_interval(LCSdata1_trans), mean_confidence_interval(LCSdata2_trans),
                  mean_confidence_interval(LCSdata3_trans), mean_confidence_interval(LCSdata4_trans),
                  mean_confidence_interval(LCSdata6_trans)], color="orange", lw = 5)
#
# fig,ax = plt.subplots()
# # plt.plot([np.mean(expData_p11),np.mean(expData_p15),np.mean(expData_p19),np.mean(expData_p28)])
# ax.errorbar(x = [11,15,19,28,35], y = [np.mean(expData_p11),np.mean(expData_p15),np.mean(expData_p19),
#                                        np.mean(expData_p28), np.mean(expData_par)],
#              yerr = [mean_confidence_interval(expData_p11), mean_confidence_interval(expData_p15),
#                      mean_confidence_interval(expData_p19),mean_confidence_interval(expData_p28),
#                      mean_confidence_interval(expData_par)], color = 'purple', lw = 5)

expDataLabels = ['Passage 11', 'Passage 15', 'Passage 19', 'Passage 28', 'Parental']
ax.set_xticklabels(expDataLabels)
ax.set_xticks([11,15,19,28,35])

sns.plt.title("DIP Rates from cFP and LSD", weight = "bold")
plt.ylabel("DIP Rate", weight = "bold")
# plt.xlabel("DIP Rate", weight = "bold")
# plt.xlim(-0.02,0.02)
# plt.ylim(0,200)
labels = ['cFP', 'LSD']
plt.legend(labels = labels)

plt.show()
quit()

mean_confidence_interval(expData_p11)

# fig,ax = plt.subplots()
# plt.plot([np.mean(expData_p11),np.mean(expData_p15),np.mean(expData_p19),np.mean(expData_p28)])
# ax.errorbar(x = [11,15,19,28,35], y = [np.mean(expData_p11),np.mean(expData_p15),np.mean(expData_p19),np.mean(expData_p28), np.mean(expData_par)],
#              yerr = [mean_confidence_interval(expData_p11), mean_confidence_interval(expData_p15),
#                      mean_confidence_interval(expData_p19),mean_confidence_interval(expData_p28), mean_confidence_interval(expData_par)])

means = (np.mean(expData_p11), np.mean(expData_p15), np.mean(expData_p19), np.mean(expData_p28),
         np.mean(expData_par))


expDataLabels = ['Passage 11', 'Passage 15', 'Passage 19', 'Passage 28', 'Parental']
ax.set_xticklabels(expDataLabels)
ax.set_xticks([11,15,19,28,35])

ind  = np.arange(5)    # the x locations for the groups
width= 0.7
labels = expDataLabels

# Pull the formatting out here
bar_kwargs = {'width':width,'color':'y','linewidth':2,'zorder':5}
err_kwargs = {'zorder':0,'fmt':None,'lw':2,'ecolor':'k'}

X = ind+width/2

fig, ax = plt.subplots()
ax.errorbar(x=[11, 15, 19, 28, 35],
            y=[np.mean(expData_p11), np.mean(expData_p15), np.mean(expData_p19), np.mean(expData_p28),
               np.mean(expData_par)],
            yerr=[mean_confidence_interval(expData_p11), mean_confidence_interval(expData_p15),
                  mean_confidence_interval(expData_p19), mean_confidence_interval(expData_p28),
                  mean_confidence_interval(expData_par)])
# ax.p1 = plt.bar(ind, menMeans, **bar_kwargs)
# ax.errs = plt.errorbar(X, menMeans, yerr=menStd, **err_kwargs)


# Custom function to draw the diff bars

def label_diff(i,j,text,X,Y):
    x = (X[i]+X[j])/2
    y = 1.1*max(Y[i], Y[j])
    dx = abs(X[i]-X[j])

    props = {'connectionstyle':'bar','arrowstyle':'-',\
                 'shrinkA':20,'shrinkB':20,'lw':2}
    ax.annotate(text, xy=(X[i],y+7), zorder=10)
    ax.annotate('', xy=(X[i],y), xytext=(X[j],y), arrowprops=props)

# Call the function
label_diff(0,1,'p=0.0370',X,means)
label_diff(0,2,'p<0.0001',X,means)
label_diff(0,3,'p=0.0025',X,means)
label_diff(0,4,'p=0.0000',X,means)



sns.set(font_scale = 2)
sns.set_style("whitegrid")
plt.show()
quit()
col_list = ["red", "green", "blue", "purple", "coral"]
col_list_palette = sns.xkcd_palette(col_list)
sns.set_palette(col_list_palette)
sns.despine(offset=10, trim=True)
# labels = ['1_0_0', '0.95_0.05_0', '0.85_0.10_0.05',
#           '0.70_0.20_0.10','0.50_0.30_0.20','0.25_0.50_0.25']
labels = ['1_0_0', '0.95_0.05_0', '0.85_0.10_0.05', '0.7_0.2_0.1' '0.25_0.50_0.25']
expDataLabels = ['Passage 11', 'Passage 15', 'Passage 19', 'Passage 28', 'Parental']

# ax = sns.violinplot(data = frames_LCStrans_mod, cut = 0, inner = 'box')
# ax = sns.violinplot(data = expData, cut = 0, inner = 'box')



# ax = sns.boxplot(data = expData, showfliers = False)
# ax = sns.pointplot(data = frames_LCStrans_mod, estimator= median)
# fig, ax = plt.subplots()
# for a in expData:
#     sns.distplot(a, ax=ax, kde=True, bins=50)
# ax.set_xlim([-0.022, 0.022])

ax.set_xticklabels(expDataLabels)
ax.set_title("DIP Rates from LSD")
ax.set_ylabel("DIP Rate", weight = "bold")
# ax.set_xlabel("Subpopulations", weight = "bold")
plt.ylim(-0.006,0.0065)
# ax.legend(labels = labels, fontsize = 16)

plt.show()
quit()



# mydict = {'1_0_0':data1,
#           '0.95_0.05_0':data2,
#           '0.90_0.10_0':data3,
#           '0.90_0.08_0.02':data4,
#           '0.85_0.10_0.05':data5,
#           '0.80_0.13_0.07':data6,
#           '0.75_0.15_0.10':data7,
#           '0.70_0.20_0.10':data8,
#           '0.55_0.30_0.15':data9,
#           '0.50_0.30_0.20':data10,
#           '0.25_0.50_0.25':data11}



# col_list = ["red", "green", "blue"]
# col_list_palette = sns.xkcd_palette(col_list)
# sns.set_palette(col_list_palette)
# sns.set_palette("husl", 8)

# ax = sns.violinplot(data = frames1, cut = 0, inner = 'box')
# ax = sns.boxplot(data = frames)
# fig, ax = plt.subplots()
# for a in frames:
#     sns.kdeplot(a, ax=ax, shade = True, bw = 0.01)
# ax.set_xlim([-0.04, 0.04])
fig, ax = plt.subplots()
for a in frames1:
    sns.distplot(a, ax=ax, kde=False, bins=50)
ax.set_xlim([-0.03, 0.03])

labels = ['1_0_0', '0.95_0.05_0', '0.90_0.10_0',
          '0.85_0.10_0.05', '0.75_0.15_0.10',
          '0.50_0.30_0.20', '0.25_0.50_0.25']
labels1 = ['0.95_0.05_0', '0.90_0.10_0.0', '0.70_0.20_0.10']
# ax.set_xticklabels(['1_0_0', '0.95_0.05_0', '0.90_0.10_0',
#                     '0.85_0.10_0.05', '0.75_0.15_0.10',
#                     '0.50_0.30_0.20', '0.25_0.50_0.25'], fontsize = 12)
ax.set_title("DIP Distributions by Initial Composition", fontsize = 20)
ax.set_ylabel("Density", fontsize = 16)
ax.set_xlabel("DIP Rate", fontsize = 16)
# ax.set(title = "DIP Distributions by Initial Composition",
#        ylabel = "Density",
#        xlabel = "DIP Rate")
ax.legend(labels = labels1)
# ax.set_xticklabels(labels1)
# ax.set_xticklabels(['1_0_0','0.95_0.05_0','0.90_0.10_0',
#                     '0.90_0.08_0.02','0.85_0.10_0.05','0.80_0.13_0.07',
#                     '0.75_0.15_0.10','0.70_0.20_0.10','0.55_0.30_0.15',
#                     '0.50_0.30_0.20','0.25_0.50_0.25'])
# ax.set_axis_labels('Distributions', 'DIP Rates')

#
# data = np.concatenate(frames)
#
# fig, axes = plt.subplots()
#
# axes.violinplot(data, showmeans=True, showextrema=True,
#                 showmedians=False)
# axes.set_title('violin plot')
# for dat in frames:
#     axes.violinplot(frames, dat,
#                     showmeans=True,
#                     showmedians=False)
#     axes.set_title('violin plot')

q75, q25 = np.percentile(data2, [75 ,25])
iqr = q75 - q25

q75, q25 = np.percentile(data3, [75 ,25])
iqr1 = q75 - q25

q75, q25 = np.percentile(data8, [75 ,25])
iqr2 = q75 - q25

print(iqr, iqr1, iqr2)

plt.show()

