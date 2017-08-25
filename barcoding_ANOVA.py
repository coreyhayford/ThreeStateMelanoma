import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import scipy.stats as sp

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

keys = []
F_vals = []
p_vals = []
DOF_btwns = []
DOF_withins = []

for key, value in dict.iteritems():
    # print key, value
    bdf = pd.DataFrame(value)
    bdf_t = pd.DataFrame.transpose(bdf)
    F = sp.stats.f_oneway(*[d[1] for d in bdf_t.iterrows()])[0]
    p = sp.stats.f_oneway(*[d[1] for d in bdf_t.iterrows()])[1]
    k = len(bdf_t)
    DOF_btwn = k-1
    N = k*len(bdf_t.columns)
    DOF_within = N-k
    keys.append(key)
    F_vals.append(F)
    p_vals.append(p)
    DOF_btwns.append(DOF_btwn)
    DOF_withins.append(DOF_within)

ANOVA_data = np.column_stack((F_vals, p_vals, DOF_btwns, DOF_withins))
print(ANOVA_data)
ANOVA_df = pd.DataFrame(ANOVA_data, index=keys)
ANOVA_df['Barcodes'] = [int(re.findall("\d+",key)[0]) for key in dict.keys()]
ANOVA_df['Experiments'] = [int(re.findall("\d+",key)[1]) for key in dict.keys()]
ANOVA_df['States'] = [int(re.findall("\d+",key)[2]) for key in dict.keys()]
ANOVA_df.columns = ['F_value', 'p_value', 'DOF_between', 'DOF_within', 'Barcodes', 'Experiments', 'States']
print(ANOVA_df)
ANOVA_df_sub = ANOVA_df[ANOVA_df['p_value'] < 0.05]
print(ANOVA_df_sub)
ANOVA_df.to_csv('ANOVA_barcoding', sep='\t')


