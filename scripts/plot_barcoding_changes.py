import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
#
# experiments = ["Experiment_%s" %i for i in range(1,11)]
# barcodes = ["Barcode_%s" %i for i in range(1,1001)]

sns.set(font_scale = 2)
sns.set_style("whitegrid")

barcoding_data_10B_10E_16S = np.load("barcoding_data_10bar10exp16states.npy")
barcoding_data_10B_50E_16S = np.load("barcoding_data_10bar50exp16states.npy")
barcoding_data_10B_100E_16S = np.load("barcoding_data_10bar100exp16states.npy")
barcoding_data_100B_10E_16S = np.load("barcoding_data_100bar10exp16states.npy")
barcoding_data_100B_50E_16S = np.load("barcoding_data_100bar50exp16states.npy")
barcoding_data_100B_100E_16S = np.load("barcoding_data_100bar100exp16states.npy")
barcoding_data_1000B_10E_16S = np.load("barcoding_data_1000bar10exp16states.npy")
barcoding_data_1000B_50E_16S = np.load("barcoding_data_1000bar50exp16states.npy")
barcoding_data_1000B_100E_16S = np.load("barcoding_data_1000bar100exp16states.npy")

barcoding_data_10B_10E_100S = np.load("barcoding_data_10bar10exp100states.npy")
barcoding_data_10B_50E_100S = np.load("barcoding_data_10bar50exp100states.npy")
barcoding_data_10B_100E_100S = np.load("barcoding_data_10bar100exp100states.npy")
barcoding_data_100B_10E_100S = np.load("barcoding_data_100bar10exp100states.npy")
barcoding_data_100B_50E_100S = np.load("barcoding_data_100bar50exp100states.npy")
barcoding_data_100B_100E_100S = np.load("barcoding_data_100bar100exp100states.npy")
barcoding_data_1000B_10E_100S = np.load("barcoding_data_1000bar10exp100states.npy")
# barcoding_data_1000B_50E_100S = np.load("barcoding_data_1000bar50exp100states.npy")
# barcoding_data_1000B_100E_100S = np.load("barcoding_data_1000bar100exp100states.npy")


df_list = [barcoding_data_10B_10E_16S,barcoding_data_10B_50E_16S,barcoding_data_10B_100E_16S,
           barcoding_data_100B_10E_16S,barcoding_data_100B_50E_16S,barcoding_data_100B_100E_16S,
           barcoding_data_1000B_10E_16S,barcoding_data_1000B_50E_16S,barcoding_data_1000B_100E_16S]

df_list_100 = [barcoding_data_10B_10E_100S, barcoding_data_10B_50E_100S, barcoding_data_10B_100E_100S,
           barcoding_data_100B_10E_100S, barcoding_data_100B_50E_100S, barcoding_data_100B_100E_100S,
           barcoding_data_1000B_10E_100S]#, barcoding_data_1000B_50E_100S, barcoding_data_1000B_100E_100S]

list_10B_10E = [barcoding_data_10B_10E_16S, barcoding_data_10B_10E_100S]
list_100B_10E = [barcoding_data_100B_10E_16S, barcoding_data_100B_10E_100S]
list_1000B_10E = [barcoding_data_1000B_10E_16S, barcoding_data_1000B_10E_100S]


list_10exp = [barcoding_data_10B_10E_16S,barcoding_data_100B_10E_16S,barcoding_data_1000B_10E_16S]
list_50exp = [barcoding_data_10B_50E_16S,barcoding_data_100B_50E_16S,barcoding_data_1000B_50E_16S]
list_100exp = [barcoding_data_10B_100E_16S,barcoding_data_100B_100E_16S,barcoding_data_1000B_100E_16S]

df_sum_list = []
# plt.figure()
for num,df in enumerate(list_1000B_10E):
    barcoding_df = pd.DataFrame.transpose(pd.DataFrame(df))
    df_sum = barcoding_df.sum(axis = 1)
    df_sum_list.append(df_sum)

    # plt.figure(num)
    # df_sort = df_sum.sort_values(ascending=False)
    # df_sort.plot(kind='bar')
    # plt.xlabel("Barcodes")
    # plt.ylabel("Count")
    # sns.distplot(df_sum, kde=False, rug=False)
    # ax = sns.violinplot(data=df_sum, cut=0, inner='box')


ax = sns.violinplot(data=df_sum_list, cut=0, inner='box')
# ax = sns.boxplot(data=df_sum_list)
plt.show()
quit()


print(barcoding_data)
print(len(barcoding_data))
barcoding_df = pd.DataFrame(barcoding_data)
barcoding_df = pd.DataFrame.transpose(barcoding_df)

# barcoding_df.index.name = "Experiment"
print(barcoding_df)
df_sum = barcoding_df.sum(axis = 1)
df_sort = df_sum.sort_values(ascending=False)
print(df_sort)
df_sort.plot(kind='bar', title = "10 Barcodes over 10 Experiments")
sns.rugplot(df_sort, color="r", axis = 'y', height = 0.01)
plt.xlabel("Barcodes")
plt.ylabel("Count")
# plt.xticks([])
# plt.show()
# quit()

barcoding_df.plot(kind='bar', stacked = True, title = "10 Barcodes over 10 Experiments", legend = False)
plt.xlabel("Barcodes")
plt.ylabel("Count")



#groupby(level = 0) #.sum().sort_values(ascending=False)

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