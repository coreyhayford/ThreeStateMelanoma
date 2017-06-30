import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import scipy.stats as sp

# import os
# import re
# path = "/Users/Corey/git/ThreeStateModel/"
# for filename in os.listdir(path):
#     if re.match("barcoding_data\d+bar\d+exp\d+states.npy", filename):
#         with np.load(os.path.join(path, filename))


bd_10B10E16S = np.load("barcoding_data_10bar10exp16states.npy")
bd_10B50E16S = np.load("barcoding_data_10bar50exp16states.npy")
bd_10B100E16S = np.load("barcoding_data_10bar100exp16states.npy")
bd_100B10E16S = np.load("barcoding_data_100bar10exp16states.npy")
bd_100B50E16S = np.load("barcoding_data_100bar50exp16states.npy")
bd_100B100E16S = np.load("barcoding_data_100bar100exp16states.npy")
bd_1000B10E16S = np.load("barcoding_data_1000bar10exp16states.npy")
bd_1000B50E16S = np.load("barcoding_data_1000bar50exp16states.npy")
bd_1000B100E16S = np.load("barcoding_data_1000bar100exp16states.npy")
bd_10B10E100S = np.load("barcoding_data_10bar10exp100states.npy")
bd_10B50E100S = np.load("barcoding_data_10bar50exp100states.npy")
bd_10B100E100S = np.load("barcoding_data_10bar100exp100states.npy")
bd_100B10E100S = np.load("barcoding_data_100bar10exp100states.npy")
bd_100B50E100S = np.load("barcoding_data_100bar50exp100states.npy")
bd_100B100E100S = np.load("barcoding_data_100bar100exp100states.npy")
bd_1000B10E100S = np.load("barcoding_data_1000bar10exp100states.npy")
bd_1000B50E100S = np.load("barcoding_data_1000bar50exp100states.npy")
bd_1000B100E100S = np.load("barcoding_data_1000bar100exp100states.npy")

bd_lib = [bd_10B10E16S,bd_10B50E16S,bd_10B100E16S,bd_100B10E16S,bd_100B50E16S,bd_100B100E16S,
          bd_1000B10E16S,bd_1000B50E16S,bd_1000B100E16S,bd_10B10E100S,bd_10B50E100S,bd_10B100E100S,
          bd_100B10E100S,bd_100B50E100S,bd_100B100E100S,bd_1000B10E100S,bd_1000B50E100S,bd_1000B100E100S]

for dat in bd_lib:
    bdf = pd.DataFrame(dat)
    # print(bdf.shape)
    bdf1 = pd.DataFrame.transpose(bdf)
    df_sum = bdf1.sum(axis = 1)
    print(bdf.shape, sp.variation(df_sum))
