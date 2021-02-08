import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

data1 = np.load("1_0_0_data.npy")
data2 = np.load("0.95_0.05_0_data.npy")
data3 = np.load("0.90_0.10_0_data.npy")
data4 = np.load("0.90_0.08_0.02_data.npy")
data5 = np.load("0.85_0.10_0.05_data.npy")
data6 = np.load("0.80_0.13_0.07_data.npy")
data7 = np.load("0.75_0.15_0.10_data.npy")
data8 = np.load("0.70_0.20_0.10_data.npy")
data9 = np.load("0.55_0.30_0.15_data.npy")
data10 = np.load("0.50_0.30_0.20_data.npy")
data11 = np.load("0.25_0.50_0.25_data.npy")

data10k_1 = np.load("0.25_0.5_0.25_10kData.npy")

frames = [data1,data2,data3, data5,
          data7, data10, data11]

frames1 = [data2, data3, data8]

mydict = {'1_0_0':data1,
          '0.95_0.05_0':data2,
          '0.90_0.10_0':data3,
          '0.90_0.08_0.02':data4,
          '0.85_0.10_0.05':data5,
          '0.80_0.13_0.07':data6,
          '0.75_0.15_0.10':data7,
          '0.70_0.20_0.10':data8,
          '0.55_0.30_0.15':data9,
          '0.50_0.30_0.20':data10,
          '0.25_0.50_0.25':data11}

sns.set_style("white")
sns.despine(offset=10, trim=True)

col_list = ["red", "green", "blue"]
col_list_palette = sns.xkcd_palette(col_list)
sns.set_palette(col_list_palette)
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

