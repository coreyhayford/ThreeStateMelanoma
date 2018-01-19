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

seg_nums = []
path = "/Users/Corey/git/ThreeStateMelanoma/stoch_expgrowth/"
file_match = "stoch_expgrowth_"

for filename in os.listdir(path):
    # print filename
    if re.match(file_match, filename):
        # print filename
        # print re.findall("\d+\.?\d+", filename)
        seg_nums.append(re.findall("\d*\.?\d+", filename))
        # dict["bd%dB%dE%dS" % re.findall("\d+", filename)] = "Hello"
# print(seg_nums)

for seg_num in seg_nums:
    # print(tuple(seg_num))
    dict["Stochastic Exponential Growth - %sCells %sDivRate %sDthRate %sSims" % tuple(seg_num)] = \
        np.load("/Users/Corey/git/ThreeStateMelanoma/stoch_expgrowth/stoch_expgrowth_%scells_%sdivrate_%sdthrate_%ssims.npy" % tuple(seg_num))

distr_values = []
# fig,ax = plt.subplots(nrows = 10)
for key,value in dict.iteritems():
    # print key, value
    # plt.figure()
    # print(key)
    # print(value)
    if value.size:
        # plt.figure()
        # sns.distplot(value, hist = False, color = "k", kde_kws={"shade": True})
        # plt.xlabel("DIP Rate", weight = "bold")
        # plt.ylabel("Density", weight = "bold")
        # plt.title(key, weight = "bold")
        distr_values.append(value)
    # plt.show()
    # quit()
# plt.show()

print(distr_values)
print(distr_values[0])
print(type(distr_values[0]))
quit()
np.savetxt('seg_distr_values.txt', distr_values, delimiter='\t')

# import cPickle as pickle
#
# with open('seg_dict.txt', 'w') as file:
#      file.write(pickle.dumps(dict))