import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

experiments = ["Experiment_%s" %i for i in range(1,11)]
barcodes = ["Barcode_%s" %i for i in range(1,1001)]

sns.set(font_scale = 2)
sns.set_style("whitegrid")
barcoding_data = np.load("barcoding_data_1000.npy")
print(barcoding_data)
print(len(barcoding_data))
barcoding_df = pd.DataFrame(barcoding_data)
# barcoding_df = pd.DataFrame.transpose(barcoding_df)

# barcoding_df.index.name = "Experiment"
print(barcoding_df)
df_sum = barcoding_df.sum()
df_sort = df_sum.sort_values(ascending=False)
print(df_sort)
df_sort.plot(kind='bar', stacked = True, title = "1000 Barcodes over 10 Experiments")
plt.xlabel("Barcodes")
plt.ylabel("Count")
plt.xticks([])
plt.show()
quit()
    #groupby(level = 0) #.sum().sort_values(ascending=False)
print(df_sort)
quit()

# quit()
barcoding_df.plot(kind='bar', stacked = True, title = "1000 Barcodes over 10 Experiments")
plt.xlabel("Barcodes")
plt.ylabel("Count")

plt.show()
quit()
# colors = ["red", "blue", "green", "black", "gold", "lightsalmon", "saddlebrown", "coral", "cyan", "aqua"]
# for i in range(len(barcoding_data)):
#     sns.barplot(range(len(barcoding_data)), [pt[i] for pt in barcoding_data])

# plt.legend(0)
    # plt.bar(range(10), [pt[i] for pt in barcoding_data], color = colors[i])
# barcoding_df = pd.DataFrame(barcoding_data, columns = barcodes)
# print(barcoding_df)
# test = pd.melt(barcoding_df, value_vars=barcodes)
# print(test)
# sns.barplot(data = test, x = "variable")
plt.show()