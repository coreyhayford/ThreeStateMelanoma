import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import scipy.stats as sp
import plotly.graph_objs as go
import plotly.plotly as py

import os
import re

sns.set(font_scale = 2)
sns.set_style("whitegrid")

dict = {}

bd_nums = []
path = "/Users/Corey/git/ThreeStateMelanoma"
file_match = "barcoding_data_\d+bar\d+exp\d+states.npy"
for filename in os.listdir(path):
    # print filename
    if re.match(file_match, filename):
        print filename
        print re.findall("\d+", filename)
        bd_nums.append(re.findall("\d+", filename))
        # dict["bd%dB%dE%dS" % re.findall("\d+", filename)] = "Hello"
print(bd_nums)

for bd_num in bd_nums:
    dict["bd%sB%sE%sS" % tuple(bd_num)] = np.load("barcoding_data_%sbar%sexp%sstates.npy" % tuple(bd_num))
        # with np.load(os.path.join(path, filename))


bd_1000B100E100S = np.load("barcoding_data_1000bar100exp100states.npy")

keys = []
# CS_stats = []
# p_vals = []
KS_stats = []
p_vals = []

for key, value in dict.iteritems():
    # print key, value
    # bdf = pd.DataFrame(value)
    # sns.boxplot(data = bdf)
    # plt.show()

    bdf1 = pd.DataFrame.transpose(bdf)
    df_means = bdf1.mean(axis=1)
    # df_sort = pd.DataFrame(df_means.sort_values(ascending=False))
    # df_sort.columns = ['Mean Count']
    # df_sort.plot(y='Mean Count', kind='bar', title="", legend=False)
    # plt.show()
    df_sum = bdf1.sum(axis = 1)
    print(df_sum)
    df_sort = pd.DataFrame(df_sum.sort_values(ascending=False))
    df_sort.columns = ['Count']
    df_sort['Relative Count'] = df_sort['Count'] / df_sort['Count'].sum()
    df_sort['Theoretical Relative Count'] = [(1./len(df_sort)) for i in range(len(df_sort))]
    # df_sort['Theoretical Relative Count'] = [(df_sort['Count'].sum()/len(df_sort)) for i in range(len(df_sort))]
    # print(df_sort)
    # quit()
    # KS_stat = sp.stats.kstest(df_sort['Relative Count'], 'uniform')[0]
    # p_val = sp.stats.kstest(df_sort['Count'], 'uniform')[1]
    KS_stat = sp.stats.ks_2samp(df_sort['Relative Count'], df_sort['Theoretical Relative Count'])[0]
    p_val = sp.stats.ks_2samp(df_sort['Count'], df_sort['Theoretical Relative Count'])[1]
    # print(KS_stat, KS_stat_2s)
    # CS_stat = sp.stats.chisquare(f_obs = df_sort['Relative Count'], f_exp= df_sort['Theoretical Relative Count'])[0]
    # p_val = sp.stats.chisquare(f_obs = df_sort['Relative Count'], f_exp= df_sort['Theoretical Relative Count'])[1]
    # print(CS_stat, p_val)
    keys.append(key)
    KS_stats.append(KS_stat)
    # CS_stats.append(CS_stat)
    p_vals.append(p_val)
    # quit()

KS_data = np.column_stack((KS_stats, p_vals))
print(KS_data)
KS_df = pd.DataFrame(KS_data, index = keys)
KS_df.columns = ['KS Stat', 'p_value']
print(KS_df)
# CS_data = np.column_stack((CS_stats, p_vals))
# print(CS_data)
# CS_df = pd.DataFrame(CS_data, index = keys)
# CS_df.columns = ['ChiSq Stat', 'p_value']
# print(CS_df)
KS_df.to_csv('KSTable_2samp_relative', sep='\t')


quit()


means = []
COVs = []
FanoFactors = []
for key, value in dict.iteritems():
    # print key, value
    bdf = pd.DataFrame(value)
    # print(bdf.shape)
    bdf1 = pd.DataFrame.transpose(bdf)
    df_sum = bdf1.sum(axis = 1)
    # print(bdf.shape, np.mean(df_sum), sp.variation(df_sum))
    means.append(np.mean(df_sum))
    COVs.append(sp.variation(df_sum))
    FanoFactors.append((np.std(df_sum)**2)/np.mean(df_sum))


means_COV = np.column_stack((means, COVs, FanoFactors))
# print(means_COV)
means_COV_df = pd.DataFrame(means_COV, index = dict.keys())
# means_COV_df['Barcodes', 'Experiments', 'States'] = [re.findall("\d+",key) for key in dict.keys()]
means_COV_df['Barcodes'] = [int(re.findall("\d+",key)[0]) for key in dict.keys()]
means_COV_df['Experiments'] = [int(re.findall("\d+",key)[1]) for key in dict.keys()]
means_COV_df['States'] = [int(re.findall("\d+",key)[2]) for key in dict.keys()]
# means_COV_df['Barcodes'] = [int(i[0]) for i in bd_nums]
# means_COV_df['Experiments'] = [int(j[1]) for j in bd_nums]
# means_COV_df['States'] = [int(k[2]) for k in bd_nums]
means_COV_df.columns = ['Mean', 'COV', 'FanoFactor', 'Barcodes', 'Experiments', 'States']
# means_COV_df.sort_values(by='COV')
print(means_COV_df)
print(means_COV_df.dtypes)

test_df1 = means_COV_df.loc[(means_COV_df['Barcodes'] == 100) & (means_COV_df['States'] == 15)]
test_df2 = means_COV_df.loc[(means_COV_df['Barcodes'] == 1000) & (means_COV_df['States'] == 15)]
test_df3 = means_COV_df.loc[(means_COV_df['Barcodes'] == 5) & (means_COV_df['States'] == 15)]
print(test_df1)
print(test_df2)
print(test_df3)
quit()

test_df = means_COV_df.loc[means_COV_df['Barcodes'] == 100]
test_df = test_df[['FanoFactor', 'States', 'Experiments']]
test_df = test_df.pivot('States', 'Experiments', 'FanoFactor')
# print(test_df)
sns.heatmap(test_df, cbar_kws={'label': 'Fano Factor'})
plt.show()
quit()
sns.heatmap(test_df, annot=True, vmin = 0, vmax=1.3, cbar_kws={'label': 'COV'})
plt.title("Barcoding Model - %d States" %sv)
plt.savefig("Barcoding_HM_comparison_%dStates.pdf" %sv)
plt.show()

# bd_10B10E16S = np.load("barcoding_data_10bar10exp16states.npy")
# bd_10B50E16S = np.load("barcoding_data_10bar50exp16states.npy")
# bd_10B100E16S = np.load("barcoding_data_10bar100exp16states.npy")
# bd_100B10E16S = np.load("barcoding_data_100bar10exp16states.npy")
# bd_100B50E16S = np.load("barcoding_data_100bar50exp16states.npy")
# bd_100B100E16S = np.load("barcoding_data_100bar100exp16states.npy")
# bd_1000B10E16S = np.load("barcoding_data_1000bar10exp16states.npy")
# bd_1000B50E16S = np.load("barcoding_data_1000bar50exp16states.npy")
# bd_1000B100E16S = np.load("barcoding_data_1000bar100exp16states.npy")
# bd_10B10E100S = np.load("barcoding_data_10bar10exp100states.npy")
# bd_10B50E100S = np.load("barcoding_data_10bar50exp100states.npy")
# bd_10B100E100S = np.load("barcoding_data_10bar100exp100states.npy")
# bd_100B10E100S = np.load("barcoding_data_100bar10exp100states.npy")
# bd_100B50E100S = np.load("barcoding_data_100bar50exp100states.npy")
# bd_100B100E100S = np.load("barcoding_data_100bar100exp100states.npy")
# bd_1000B10E100S = np.load("barcoding_data_1000bar10exp100states.npy")
# bd_1000B50E100S = np.load("barcoding_data_1000bar50exp100states.npy")
# bd_1000B100E100S = np.load("barcoding_data_1000bar100exp100states.npy")
#
# bd_lib = [bd_10B10E16S,bd_10B50E16S,bd_10B100E16S,bd_100B10E16S,bd_100B50E16S,bd_100B100E16S,
#           bd_1000B10E16S,bd_1000B50E16S,bd_1000B100E16S,bd_10B10E100S,bd_10B50E100S,bd_10B100E100S,
#           bd_100B10E100S,bd_100B50E100S,bd_100B100E100S,bd_1000B10E100S,bd_1000B50E100S,bd_1000B100E100S]

# means = []
# COVs = []
#
# for dat in bd_lib:
#     bdf = pd.DataFrame(dat)
#     # print(bdf.shape)
#     bdf1 = pd.DataFrame.transpose(bdf)
#     df_sum = bdf1.sum(axis = 1)
#     # print(bdf.shape, np.mean(df_sum), sp.variation(df_sum))
#     means.append(np.mean(df_sum))
#     COVs.append(sp.variation(df_sum))

# means_COV = np.column_stack((means, COVs))
# # print(means_COV)
# means_COV_df = pd.DataFrame(means_COV)
# print(means_COV_df)
# quit()
# print(means_COV_df.loc[(means_COV_df['States'] == 5) & (means_COV_df['Barcodes'] >= 100)])

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), sp.stats.sem(a)
    h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
    # return m, m-h, m+h
    return h
    # print(m, m-h, m+h)

sub1 = means_COV_df.loc[(means_COV_df['States'] == 15) & (means_COV_df['Barcodes'] == 100)]
sub2 = means_COV_df.loc[(means_COV_df['States'] == 15) & (means_COV_df['Barcodes'] == 1000)]
sub3 = means_COV_df.loc[(means_COV_df['States'] == 15) & (means_COV_df['Barcodes'] == 5)]

print(sub1)
print(sub2)
print(sub3)
#
# dd[ order(-dd[,4], dd[,1]), ]
# dd[with(dd, order(-z, b)), ]
# sub1_o = sub1[with(sub1, order("Experients")), ]

sub1_s = sub1.sort_values(by = "Experiments")
sub2_s = sub2.sort_values(by = "Experiments")
sub3_s = sub3.sort_values(by = "Experiments")

plt.plot(sub1_s['Experiments'], sub1_s['COV'], color = "red", lw = 5)
plt.plot(sub2_s['Experiments'], sub2_s['COV'], color = "blue", lw = 5)
plt.plot(sub3_s['Experiments'], sub3_s['COV'], color = "green", lw = 5)

sns.plt.title("DNA Barcoding Variability", weight = "bold")
plt.ylabel("COV", weight = "bold")
plt.xlabel("Number of Experiments", weight = "bold")
# plt.xlim(-0.02,0.02)
# plt.ylim(0,200)
labels = ['15 States 100 Barcodes', '15 States 1000 Barcodes', '15 States 5 Barcodes']
plt.legend(labels = labels)



# sub1 = means_COV_df[means_COV_df[]]
# fig,ax = plt.subplots()
# plt.plot(
#     [np.mean(expData_p11),np.mean(expData_p15),np.mean(expData_p19),np.mean(expData_p28)])
#
# quit()

state_val = [5,10,15,20,25,16,100]
for sv in state_val:
    # print sv
    plt.figure("%d States" % sv)
    test_df = means_COV_df.loc[(means_COV_df['States'] == sv) & (means_COV_df['Barcodes'] >= 100)]
    test_df = test_df[['COV', 'Barcodes', 'Experiments']]
    test_df = test_df.pivot('Barcodes', 'Experiments', 'COV')
    # print(test_df)
    sns.heatmap(test_df, annot=True, vmin = 0, vmax=1.3, cbar_kws={'label': 'COV'})
    plt.title("Barcoding Model - %d States" %sv)
    plt.savefig("Barcoding_HM_comparison_%dStates.pdf" %sv)
plt.show()

    # sns.heatmap(means_COV_df, annot=True, fmt="d", linewidths=.5)
# plt.show()




