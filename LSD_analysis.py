import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import scipy.stats as sp

import os
import re

# sns.set(font_scale = 1.25)
sns.set_style("whitegrid")

dict = {}

lsd_nums = []
path = "/Users/Corey/git/ThreeStateMelanoma/LSD_spikeIn"
file_match = "LSDspikein_\w+DIPs_\d+\.\d+resistant_\d+cells_\d+wells.npy"
for filename in os.listdir(path):
    # print filename
    # print "..."
    if re.match(file_match, filename):
        # print filename
        # print re.findall("\d+", filename)
        dip_type = re.search('LSDspikein_(.*)DIPs_\d+\.\d+resistant_\d+cells_\d+wells.npy', filename).group(1)
        perc_cells_wells = re.findall("\d+\.?\d*", filename)
        # perc_resist = re.findall("\d+\.\d+", filename)
        # cells_wells = re.findall("\d*[.,]?\d*", filename)
        # print(dip_type)
        print(perc_cells_wells)
        perc_cells_wells = [float(i) for i in perc_cells_wells]
        # print(test)
        # quit()
        perc_cells_wells.insert(0, dip_type)
        print(perc_cells_wells)
        # quit()
        lsd_nums.append(perc_cells_wells)
# dict["bd%dB%dE%dS" % re.findall("\d+", filename)] = "Hello"
print(lsd_nums)
# quit()
for lsd_num in lsd_nums:
    os.chdir(path)
    dict["LSDspikein_%sDIPs_%.6fresistant_%dcells_%dwells" % tuple(lsd_num)] = \
        np.load("LSDspikein_%sDIPs_%.6fresistant_%dcells_%dwells.npy" % tuple(lsd_num))
# with np.load(os.path.join(path, filename))

lsd_medians = []
lsd_datpts = []
dat_type = []
perc_resist = []
num_cells = []
num_wells = []

for key, value in dict.iteritems():
    # print key, value
    dt = re.search('LSDspikein_(.*)DIPs_\d+\.\d+resistant_\d+cells_\d+wells', key).group(1)
    # print(dt)
    dat_type.append(dt)

    perc, cells, wells = re.findall("\d+\.?\d*", key)
    # print(perc)
    # print(cells)
    # print(wells)
    perc_resist.append(perc)
    num_cells.append(cells)
    num_wells.append(wells)

    lsd_dat = value
    lsd_datpts.append(lsd_dat)
    # print(lsd_dat)
    # lsd_med = lsd_dat.median()
    # print(lsd_med)
    # lsd_medians.append(lsd_med)
    # quit()

dat_mat = np.column_stack((dat_type, perc_resist, num_cells, num_wells, lsd_datpts))
dat_mat_df = pd.DataFrame(dat_mat)#, index = dict.keys())
dat_mat_df.columns = ['Type', 'Percent Resistant', 'Cells per Well', 'Wells per Plate', 'DIP Rates']
print(dat_mat_df)

dat_mat_normal = dat_mat_df[dat_mat_df['Type'] == 'normal']
dat_mat_resistant = dat_mat_df[dat_mat_df['Type'] == 'resistant']

# print(dat_mat_normal)
# print(dat_mat_resistant)

# dat_mat_normal.to_csv('testnorm.csv')
# dat_mat_resistant.to_csv('testresist.csv')

test = pd.merge(dat_mat_normal, dat_mat_resistant, how = 'left',
                on = ['Percent Resistant', 'Cells per Well', 'Wells per Plate'])
#
test = pd.DataFrame(test)
# test.to_csv('test.csv', sep = ',')
print(np.append(test['DIP Rates_x'][0], (test['DIP Rates_y'][0])))
print(len(np.append(test['DIP Rates_x'][0], (test['DIP Rates_y'][0]))))

test['DIP Rates_all'] = [np.append(test['DIP Rates_x'][i], test['DIP Rates_y'][i]) for i in range(len(test))]

test['Cells'] = test['Cells per Well'].values.astype(np.int64)
print(test['Cells per Well'].astype(np.int64))

test['Wells'] = test['Wells per Plate'].values.astype(np.int64)

print(test)
# print(type(test['Cells per Well'].values))
# quit()
#
# print(type(int(test['Cells per Well'][1])))
# # quit()
# print(type(test['Cells per Well'][1]))
# test['Cells per Well'].astype('int')
# print(type(test['Cells per Well'][1]))
#
# quit()
# test['Cells per Well'] = test['Cells per Well'].to_numeric(test['Cells per Well'], errors = 'coerce').fillna(0).astype(np.int64)

test['DIP Median'] = [np.median(test['DIP Rates_all'][j]) for j in range(len(test))]

test['DIP Maximum'] = [np.max(test['DIP Rates_all'][k]) for k in range(len(test))]

final_test = np.column_stack((test['Percent Resistant'], test['Cells'],
                             test['Wells'], test['DIP Rates_all'],
                              test['DIP Median'], test['DIP Maximum']))

final_test_df = pd.DataFrame(final_test)
final_test_df.columns = ['Percent Resistant', 'Cells', 'Wells', 'DIP Rates', 'DIP Median', 'DIP Maximum']

# final_test_df['Cells'] = int(final_test_df['Cells'])
final_test_df = final_test_df.sort_values(['Cells'], ascending=False)

# print(final_test_df)
# quit()


### Number of Wells Figure
final_test_df_wellsSorted = final_test_df.sort_values(['Wells'], ascending=False)
final_test_df_10cells = final_test_df_wellsSorted[final_test_df_wellsSorted['Cells'] == 10]

# final_test_df_10cells = final_test_df_10cells[['Cells', 'Wells', 'Percent Resistant', 'DIP Median']]

# print(final_test_df_10cells)
# quit()

res1 = final_test_df_10cells.pivot(index='Wells', columns='Percent Resistant', values='DIP Maximum')
res1.fillna(value = np.nan, inplace=True)

print(res1)

sns.set(rc={'figure.figsize': (15,8)})

ax = plt.axes()
hm = sns.heatmap(res1, fmt="g", cmap='viridis',
                 cbar_kws={'label': 'Maximum DIP Rate'},
                 annot_kws={'size': 10})
ax.set_title('Effect of Number of Wells on DIP Rate', weight = "bold")

hm_fig = hm.get_figure()
hm_fig.savefig('numWells_percentResistant_byDIP_10cells.pdf', dpi = 600)

quit()
### Seeding density Figure
final_test_df_1536 = final_test_df[final_test_df['Wells'] == "1536"]
print(final_test_df_1536)

res = final_test_df_1536.pivot(index = 'Cells', columns= 'Percent Resistant', values = 'DIP Median')
print(res)
# res.sortlevel(level = 1, ascending=False, inplace=True)
res.fillna(value = np.nan, inplace=True)
# res_test = res.dropna(axis = 1, how = 'any')

print(res)


# res.to_csv('res.csv')
# sns.heatmap(result, annot=True, fmt="g", cmap='viridis')

sns.set(rc={'figure.figsize': (15,8)})

ax = plt.axes()
hm = sns.heatmap(res, fmt="g", cmap='coolwarm',
                 cbar_kws={'label': 'Median DIP Rate'},
                 annot_kws={'size': 10})
ax.set_title('Effect of Seeding Density on DIP Rate', weight = "bold")

hm_fig = hm.get_figure()
hm_fig.savefig('seedingDensity_percentResistant_byDIP_1536.pdf', dpi = 600)
# plt.show()