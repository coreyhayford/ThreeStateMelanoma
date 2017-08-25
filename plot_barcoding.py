import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

experiments = ["Experiment_%s" %i for i in range(1,11)]
barcodes = ["Barcode_%s" %i for i in range(1,1001)]

sns.set(font_scale = 2)
sns.set_style("whitegrid")
barcoding_data = np.load("barcoding_data_100bar5exp15states.npy")
# print(barcoding_data)
# print(len(barcoding_data))
barcoding_df = pd.DataFrame(barcoding_data)
barcoding_df = pd.DataFrame.transpose(barcoding_df)

# barcoding_df.index.name = "Experiment"
# print(barcoding_df)

barcoding_df[barcoding_df >0] = 1.0
print(type(barcoding_df))
barcoding_df.sum(axis = 1)
# print(test)
# test = test.sum(axis = 1)
# print(test)
# test_df = barcoding_df[barcoding_df > 0] = 1
# print test_df
# test_df = re.sub()
barcoding_df.columns = ["Experiment 1","Experiment 2","Experiment 3","Experiment 4","Experiment 5"]
barcoding_df['Sum'] = barcoding_df.sum(axis = 1)
# test = barcoding_df.sum()
# barcoding_df.replace('^[1-9][0-9]*$', '0')
print(barcoding_df)

test = barcoding_df

for exp in test:
    test[exp] = np.where(test[exp] == 0, np.NaN, test[exp])
    test[exp] = np.where(test[exp] == 1, test['Sum']-test[exp], test[exp])

print(test)
test1 = test.apply(pd.Series.value_counts)
del test1["Sum"]

test1 = pd.DataFrame.transpose(test1)

test1['Sum'] = test1.sum(axis = 1)
del test1[5.0]
test1.index = [1,2,3,4,5]
print(test1)
a = test1[0]/test1['Sum'] * 100
b = test1[1]/test1['Sum'] * 100
c = test1[2]/test1['Sum'] * 100
d = test1[3]/test1['Sum'] * 100
e = test1[4]/test1['Sum'] * 100

dat = {'0':a,'1':b,'2':c,'3':d,'4':e}
index = [1,2,3,4,5]
test2 = pd.DataFrame(data=dat, index=index)

print(test2)

test2.plot(kind='bar', stacked = True, legend = True)
plt.legend(loc = 9, ncol = 5, title = "Other Replicates Sharing Unique Barcodes", frameon = True)
plt.xlabel("Experimental Replicate", weight = 'bold')
plt.ylabel("Percentage of Unique Barcodes", weight = 'bold')
plt.title("Percent Unique Barcodes by Replicate", weight = 'bold')
plt.show()
quit()
test1 = pd.DataFrame.transpose(test1)
test1.index = [1,2,3,4,5,'Sum']
# test1.loc[[1,2,3,4,5]].div(test1.loc['Sum'], axis = 0)
# test1.iloc[0:4,:].div(test1.iloc['Sum',:], axis=1)

print(test1)

quit()
test1[[0.0,1.0,2.0,3.0,4.0]].div(test1['Sum'], axis=0)
# test1 = pd.DataFrame.transpose(test1)
print(test1)
quit()
for row in test1.rows:
    np.where(test1[row] > 1, test1['Sum']/test1[row], test1[row])
print(test1)
quit()
# test1.index = [1,2,3,4,5]
for i in test1:
    print i
    test1[i] = np.where(test[i] > 0, test[i]/test['Sum'], test[i])
print(test1)
quit()
test1.append(test1.sum(numeric_only=True), ignore_index=True)
print(test1)
test1.plot(kind='bar', stacked = True, title = "Percent Unique Barcodes by Replicate", legend = False)
plt.show()
quit()
for row in range(len(barcoding_df)):
    for col in barcoding_df.columns:
        print row
        print col
        if barcoding_df[row,col] == 0:
            print ("Yay!")

quit()


#             barcoding_df[row,col] == 1
#
# print(barcoding_df)

quit()
df_sum = barcoding_df.sum(axis = 1)
df_sort = pd.DataFrame(df_sum.sort_values(ascending=False))
df_sort.columns = ['Count']
df_sort['Relative Count'] = df_sort['Count'] / df_sort['Count'].sum()

print(df_sort)

df_sort.plot(y = 'Relative Count', kind='bar', title = "100 Barcodes over 5 Experiments - 15 states", legend = False)
plt.plot([0, 100], [0.01, 0.01], "k--", lw = 5)

# sns.rugplot(df_sort, color="r", axis = 'y', height = 0.01)
# plt.xlabel("Barcodes", weight = "bold")
plt.ylabel("Relative Count", weight = "bold")
plt.xticks([])
# plt.show()
# quit()

barcoding_df.plot(kind='bar', stacked = True, title = "100 Barcodes over 100 Experiments - 15 states", legend = False)
plt.xlabel("Barcodes", weight = "bold")
plt.ylabel("Count", weight = "bold")

plt.show()
quit()

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