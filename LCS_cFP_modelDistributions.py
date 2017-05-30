import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Load model produced data
LCSdata1_trans = np.load("1.0_0.0_0.0_10kData.npy")
LCSdata2_trans = np.load("0.95_0.05_0.0_10kData.npy")
LCSdata3_trans = np.load("0.85_0.1_0.05_10kData.npy")
LCSdata4_trans = np.load("0.7_0.2_0.1_10kData.npy")
LCSdata5_trans = np.load("0.5_0.3_0.2_10kData.npy")
LCSdata6_trans = np.load("0.25_0.5_0.25_10kData.npy")

cFPdata1_trans = np.load("1.0_0.0_0.0_10kData_cFP.npy")
cFPdata2_trans = np.load("0.95_0.05_0.0_10kData_cFP.npy")
cFPdata3_trans = np.load("0.85_0.1_0.05_10kData_cFP.npy")
cFPdata4_trans = np.load("0.7_0.2_0.1_10kData_cFP.npy")
cFPdata5_trans = np.load("0.5_0.3_0.2_10kData_cFP.npy")
cFPdata6_trans = np.load("0.25_0.5_0.25_10kData_cFP.npy")

LCSdata1_notrans = np.load("1.0_0.0_0.0_10kData_LCS_noTrans.npy")
LCSdata2_notrans = np.load("0.95_0.05_0.0_10kData_LCS_noTrans.npy")
LCSdata3_notrans = np.load("0.85_0.1_0.05_10kData_LCS_noTrans.npy")
LCSdata4_notrans = np.load("0.7_0.2_0.1_10kData_LCS_noTrans.npy")
LCSdata5_notrans = np.load("0.5_0.3_0.2_10kData_LCS_noTrans.npy")
LCSdata6_notrans = np.load("0.25_0.5_0.25_10kData_LCS_noTrans.npy")

cFPdata1_notrans = np.load("1.0_0.0_0.0_10kData_cFP_noTrans.npy")
cFPdata2_notrans = np.load("0.95_0.05_0.0_10kData_cFP_noTrans.npy")
cFPdata3_notrans = np.load("0.85_0.1_0.05_10kData_cFP_noTrans.npy")
cFPdata4_notrans = np.load("0.7_0.2_0.1_10kData_cFP_noTrans.npy")
cFPdata5_notrans = np.load("0.5_0.3_0.2_10kData_cFP_noTrans.npy")
cFPdata6_notrans = np.load("0.25_0.5_0.25_10kData_cFP_noTrans.npy")

expData_p11 = np.loadtxt("experimental_psg11.txt", delimiter='\t')
expData_p15 = np.loadtxt("experimental_psg15.txt", delimiter='\t')
expData_p19 = np.loadtxt("experimental_psg19.txt", delimiter='\t')
expData_par = np.loadtxt("experimental_parental.txt", delimiter='\t')


frames_LCStrans = [LCSdata1_trans,LCSdata2_trans,LCSdata3_trans,LCSdata4_trans,LCSdata5_trans,LCSdata6_trans]
frames_LCStrans_mod = [LCSdata1_trans,LCSdata2_trans,LCSdata3_trans,LCSdata4_trans]
frames_cFPtrans = [cFPdata1_trans,cFPdata2_trans,cFPdata3_trans,cFPdata4_trans,cFPdata5_trans,cFPdata6_trans]
frames_LCSnotrans = [LCSdata1_notrans,LCSdata2_notrans,LCSdata3_notrans,LCSdata4_notrans,LCSdata5_notrans,LCSdata6_notrans]
frames_cFPnotrans = [cFPdata1_notrans,cFPdata2_notrans,cFPdata3_notrans,cFPdata4_notrans,cFPdata5_notrans,cFPdata6_notrans]
expData = [expData_p11, expData_p15, expData_p19, expData_par]

sns.set(font_scale = 2)
sns.set_style("whitegrid")
col_list = ["red", "green", "blue", "purple"]
col_list_palette = sns.xkcd_palette(col_list)
sns.set_palette(col_list_palette)
sns.despine(offset=10, trim=True)
# labels = ['1_0_0', '0.95_0.05_0', '0.85_0.10_0.05',
#           '0.70_0.20_0.10','0.50_0.30_0.20','0.25_0.50_0.25']
labels = ['1_0_0', '0.95_0.05_0', '0.85_0.10_0.05', '0.25_0.50_0.25']
expDataLabels = ['Passage 11', 'Passage 15', 'Passage 19', 'Parental']

ax = sns.violinplot(data = frames_LCStrans_mod, cut = 0, inner = 'box')
# ax = sns.violinplot(data = expData, cut = 0, inner = 'box')

# ax = sns.boxplot(data = frames_cFPnotrans)
# fig, ax = plt.subplots()
# for a in frames_LCStrans:
#     sns.distplot(a, ax=ax, kde=False, bins=50)
# ax.set_xlim([-0.02, 0.02])

ax.set_xticklabels(labels)
# ax.set_title("DIP Distributions by Initial Composition")
ax.set_ylabel("DIP Rate", weight = "bold")
ax.set_xlabel("Subpopulations", weight = "bold")
plt.ylim(-0.022,0.022)
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

