import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp
import pandas as pd
import re
from numpy import median
from scipy.stats import gaussian_kde

sns.set(font_scale = 1.25)
sns.set_style("whitegrid")

LSD_100_0_0_1cell = np.random.choice(np.load("1.0_0.0_0.0_10kData_LSD_1cell.npy"), 1000)
LSD_100_0_0_5cells = np.load("1.0_0.0_0.0_1kData_LSD_5cells.npy")
LSD_100_0_0_10cells = np.load("1.0_0.0_0.0_1kData_LSD_10cells.npy")
LSD_100_0_0_25cells = np.load("1.0_0.0_0.0_1kData_LSD_25cells.npy")

print(np.median(LSD_100_0_0_1cell))
print(np.median(LSD_100_0_0_5cells))
print(np.median(LSD_100_0_0_10cells))
print(np.median(LSD_100_0_0_25cells))


bins = np.histogram(np.hstack((LSD_100_0_0_1cell, LSD_100_0_0_5cells,
                               LSD_100_0_0_10cells, LSD_100_0_0_25cells)),
                    bins=100)[1]

sns.distplot(LSD_100_0_0_1cell, color="brown", hist=False, bins=bins, kde_kws={"shade": True}, label="1 cell")
sns.distplot(LSD_100_0_0_5cells, color="magenta", hist=False, bins=bins, kde_kws={"shade": True}, label="5 cells")
sns.distplot(LSD_100_0_0_10cells, color="yellow", hist=False, bins=bins, kde_kws={"shade": True}, label="10 cells")
sns.distplot(LSD_100_0_0_25cells, color="black", hist=False, bins=bins, kde_kws={"shade": True}, label="25 cells")

# LSD_100_0_0_25cells = np.load("1.0_0.0_0.0_1kData_LSD_25cells.npy")
# LSD_85_10_5_25cells = np.load("0.85_0.1_0.05_1kData_LSD_25cells.npy")
# LSD_70_20_10_25cells = np.load("0.7_0.2_0.1_1kData_LSD_25cells.npy")
# LSD_25_50_25_25cells = np.load("0.25_0.5_0.25_1kData_LSD_25cells.npy")
#
# bins = np.histogram(np.hstack((LSD_100_0_0_25cells, LSD_85_10_5_25cells,
#                                LSD_70_20_10_25cells, LSD_25_50_25_25cells)),
#                     bins=100)[1]
#
# sns.distplot(LSD_100_0_0_25cells, color="red", hist=False, bins=bins, kde_kws={"shade": True}, label="Low")
# sns.distplot(LSD_85_10_5_25cells, color="green", hist=False, bins=bins, kde_kws={"shade": True}, label="Medium-Low")
# sns.distplot(LSD_70_20_10_25cells, color="blue", hist=False, bins=bins, kde_kws={"shade": True}, label="Medium-High")
# sns.distplot(LSD_25_50_25_25cells, color="magenta", hist=False, bins=bins, kde_kws={"shade": True}, label="High")

# LSD_100_0_0 = np.load("1.0_0.0_0.0_10kData_LSD_5cells.npy")
# LSD_85_10_05 = np.load("0.85_0.1_0.05_10kData_LSD_5cells.npy")
# LSD_70_20_10 = np.load("0.7_0.2_0.1_10kData_LSD_5cells.npy")
# LSD_25_50_25 = np.load("0.25_0.5_0.25_10kData_LSD_5cells.npy")

# bins = np.histogram(np.hstack((LSD_100_0_0, LSD_85_10_05, LSD_70_20_10, LSD_25_50_25)),
#                     bins=100)[1]

# sns.distplot(LSD_100_0_0, color="red", kde=False, bins=bins, hist_kws={"alpha":0.35}, label="Low")
# sns.distplot(LSD_85_10_05, color="green", kde=False, bins=bins, hist_kws={"alpha":0.35}, label="Medium-Low")
# sns.distplot(LSD_70_20_10, color="blue", kde=False, bins=bins, hist_kws={"alpha":0.35}, label="Medium-High")
# sns.distplot(LSD_25_50_25, color="magenta", kde=False, bins=bins, hist_kws={"alpha":0.35}, label="High")

# sns.distplot(LSD_100_0_0, color="red", hist=False, bins=bins, kde_kws={"shade": True}, label="Low")
# sns.distplot(LSD_85_10_05, color="green", hist=False, bins=bins, kde_kws={"shade": True}, label="Medium-Low")
# sns.distplot(LSD_70_20_10, color="blue", hist=False, bins=bins, kde_kws={"shade": True}, label="Medium-High")
# sns.distplot(LSD_25_50_25, color="magenta", hist=False, bins=bins, kde_kws={"shade": True}, label="High")


# cFP_100_0_0 = np.load("1.0_0.0_0.0_10kData_LSD_1cell.npy")
# cFP_85_10_05 = np.load("0.85_0.1_0.05_10kData_LSD_1cell.npy")
# cFP_70_20_10 = np.load("0.7_0.2_0.1_10kData_LSD_1cell.npy")
# cFP_25_50_25 = np.load("0.25_0.5_0.25_10kData_LSD_1cell.npy")
#
# bins = np.histogram(np.hstack((cFP_100_0_0, cFP_85_10_05, cFP_70_20_10, cFP_25_50_25)),
#                     bins=100)[1]

# sns.distplot(cFP_100_0_0, color="red", kde=False, bins=bins, hist_kws={"alpha":0.35}, label="Low")
# sns.distplot(cFP_85_10_05, color="green", kde=False, bins=bins, hist_kws={"alpha":0.35}, label="Medium-Low")
# sns.distplot(cFP_70_20_10, color="blue", kde=False, bins=bins, hist_kws={"alpha":0.35}, label="Medium-High")
# sns.distplot(cFP_25_50_25, color="magenta", kde=False, bins=bins, hist_kws={"alpha":0.35}, label="High")

# sns.distplot(cFP_100_0_0, color="red", hist=False, bins=bins, kde_kws={"shade": True}, label="Low")
# sns.distplot(cFP_85_10_05, color="green", hist=False, bins=bins, kde_kws={"shade": True}, label="Medium-Low")
# sns.distplot(cFP_70_20_10, color="blue", hist=False, bins=bins, kde_kws={"shade": True}, label="Medium-High")
# sns.distplot(cFP_25_50_25, color="magenta", hist=False, bins=bins, kde_kws={"shade": True}, label="High")


plt.title("DIP Rate Distributions", weight="bold")
plt.xlabel("DIP Rate", weight="bold")
# plt.ylabel("Frequency", weight="bold")
plt.ylabel("Density", weight="bold")
plt.tick_params(labelsize = 10)
plt.xlim(-0.03, 0.03)
plt.ylim(0,600)
# plt.legend(title="Diversity", loc="upper right")
plt.legend(title="Number of Cells Seeded", loc="upper right")
plt.savefig("LSDSimulatedComparisons_100_0_0_1to25cells_75to125hr_kde.pdf")
# plt.savefig("LSDSimulatedComparisons_first25hr_hist.pdf")
# plt.savefig("cFPSimulatedComparisons_first25hr_kde.pdf")
plt.show()