from pysb import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sp
from pysb.simulator.bng import BngSimulator
from pysb.bng import generate_equations
import pandas as pd
import seaborn as sns
import math
import matplotlib

# current_cmap = matplotlib.cm.get_cmap()
# current_cmap.set_bad(color='lightgrey')
# # fig = plt.figure(figsize = (15,9)) # width x height
# # ax1 = fig.add_subplot(3, 2, 1) # row, column, position
# # ax2 = fig.add_subplot(3, 2, 2)
# # ax3 = fig.add_subplot(3, 2, 3)
# # ax4 = fig.add_subplot(3, 2, 4)
# # ax5 = fig.add_subplot(3, 2, 5)
# # ax6 = fig.add_subplot(3, 2, 6)


##### THIS IS FOR CONVERSION OF PICKLED DATA TO EASY INPUT TO R #####
# df_DS8 = pd.read_pickle('PC9-DS8_param-scan_twoState.pkl')
# # df_DS8 = df_DS8.round({'death rate': 5, 'division rate': 5, 'DIP rate': 5})
# df_DS8['cell.line'] = 'PC9.DS8'
# print(df_DS8)
# print(type(df_DS8))
#
# DS8_dipdf = pd.DataFrame(df_DS8['DIP rate'].values.tolist(), columns=['DIP1','DIP2'])
# DS8_divdf = pd.DataFrame(df_DS8['division rate'].values.tolist(), columns=['div1','div2'])
# DS8_dthdf = pd.DataFrame(df_DS8['death rate'].values.tolist(), columns=['dth1','dth2'])
# DS8_rest = df_DS8[['p-value',  'cell.line']]
#
# print(DS8_dipdf)
# print(DS8_divdf)
# print(DS8_dthdf)
# print(DS8_rest)
#
# DS8 = pd.concat([DS8_rest, DS8_dipdf, DS8_divdf,DS8_dthdf], axis = 1)
# print(DS8)
#
# DS8['cell.line'] = np.where(DS8['p-value']>0.1, 'PC9-DS8', 'not.assigned')
# DS8['param.pair'] = range(DS8.shape[0])
# print(DS8)
# DS8_sig1 = DS8[['p-value', 'cell.line', 'DIP1', 'div1', 'dth1', 'param.pair']]
# DS8_sig2 = DS8[['p-value', 'cell.line', 'DIP2', 'div2', 'dth2', 'param.pair']]
#
# DS8_sig1.rename(columns={'DIP1': 'DIP Rate',
#                          'div1': 'Division Rate',
#                          'dth1': 'Death Rate'},
#                  inplace=True)
# DS8_sig2.rename(columns={'DIP2': 'DIP Rate',
#                          'div2': 'Division Rate',
#                          'dth2': 'Death Rate'},
#                 inplace=True)
#
# DS8_sig1['Cell Line'] = np.where(DS8_sig1['cell.line'] == "PC9-DS8", 'PC9-DS8.1', 'not.assigned')
# DS8_sig2['Cell Line'] = np.where(DS8_sig2['cell.line'] == "PC9-DS8", 'PC9-DS8.2', 'not.assigned')
#
# print(DS8_sig1)
# print(DS8_sig2)
#
# DS8_sig_all = pd.concat([DS8_sig1, DS8_sig2])
# print(DS8_sig_all)
# DS8_sig_all.to_csv('DS8_twoState_tile_extra.csv')
# quit()
#######


##### PLOTTING DATA ######
# result_DS8 = df_DS8.pivot(index='division rate', columns='DIP rate', values='p-value')
# plt.figure('DS8', figsize=(3,4))
# a = sns.heatmap(result_DS8, cbar_kws={'label': 'p-value'}, cmap = 'jet')
# a.set_xticklabels(a.get_xticklabels(), rotation = 45, fontsize = 8)
# a.set_yticklabels(a.get_yticklabels(), rotation = 45, fontsize = 8)
# plt.title("DS8", weight = "bold")
# plt.xlabel("DIP Rate (doublings/hour)")
# plt.ylabel("Division Rate (per hour)")
# plt.tight_layout()
# # plt.savefig('DS8_paramscan_LR_tiny.png')
# # plt.savefig('DS8_paramscan_TR.png')
#
# df_DS6 = pd.read_pickle('PC9-DS6_param-scan_tighterRange.pkl')
# df_DS6 = df_DS6.round({'death rate': 5, 'division rate': 5, 'DIP rate': 5})
# df_DS6['cell.line'] = 'PC9.DS6'
# result_DS6 = df_DS6.pivot(index='division rate', columns='DIP rate', values='p-value')
# plt.figure('DS6', figsize=(3,4))
# a = sns.heatmap(result_DS6, cbar_kws={'label': 'p-value'}, cmap = 'jet')
# a.set_xticklabels(a.get_xticklabels(), rotation = 45, fontsize = 8)
# a.set_yticklabels(a.get_yticklabels(), rotation = 45, fontsize = 8)
# plt.title("DS6", weight = "bold")
# plt.xlabel("DIP Rate (doublings/hour)")
# plt.ylabel("Division Rate (per hour)")
# plt.tight_layout()
# # plt.savefig('DS6_paramscan_TR1_tiny.png')
# # plt.savefig('DS6_paramscan_TR1.png')
#
# df_DS9 = pd.read_pickle('PC9-DS9_param-scan_tighterRange.pkl')
# df_DS9 = df_DS9.round({'death rate': 5, 'division rate': 5, 'DIP rate': 5})
# df_DS9['cell.line'] = 'PC9.DS9'
# result_DS9 = df_DS9.pivot(index='division rate', columns='DIP rate', values='p-value')
# plt.figure('DS9', figsize=(3,4))
# a = sns.heatmap(result_DS9, cbar_kws={'label': 'p-value'}, cmap = 'jet')
# a.set_xticklabels(a.get_xticklabels(), rotation = 45, fontsize = 8)
# a.set_yticklabels(a.get_yticklabels(), rotation = 45, fontsize = 8)
# plt.title("DS9", weight = "bold")
# plt.xlabel("DIP Rate (doublings/hour)")
# plt.ylabel("Division Rate (per hour)")
# plt.tight_layout()
# # plt.savefig('DS9_paramscan_TR1_tiny.png')
# # plt.savefig('DS9_paramscan_TR1.png')
#
# df_DS1 = pd.read_pickle('PC9-DS1_param-scan_tighterRange.pkl')
# df_DS1 = df_DS1.round({'death rate': 5, 'division rate': 5, 'DIP rate': 5})
# df_DS1['cell.line'] = 'PC9.DS1'
# result_DS1 = df_DS1.pivot(index='division rate', columns='DIP rate', values='p-value')
# plt.figure('DS1', figsize=(3,4))
# a = sns.heatmap(result_DS1, cbar_kws={'label': 'p-value'}, cmap = 'jet')
# a.set_xticklabels(a.get_xticklabels(), rotation = 45, fontsize = 8)
# a.set_yticklabels(a.get_yticklabels(), rotation = 45, fontsize = 8)
# plt.title("DS1", weight = "bold")
# plt.xlabel("DIP Rate (doublings/hour)")
# plt.ylabel("Division Rate (per hour)")
# plt.tight_layout()
# # plt.savefig('DS1_paramscan_TR1_tiny.png')
# # plt.savefig('DS1_paramscan_TR1.png')
#
# df_DS3 = pd.read_pickle('PC9-DS3_param-scan_tighterRange.pkl')
# df_DS3 = df_DS3.round({'death rate': 5, 'division rate': 5, 'DIP rate': 5})
# df_DS3['cell.line'] = 'PC9.DS3'
# result_DS3 = df_DS3.pivot(index='division rate', columns='DIP rate', values='p-value')
# plt.figure('DS3', figsize=(3,4))
# a = sns.heatmap(result_DS3, cbar_kws={'label': 'p-value'}, cmap = 'jet')
# a.set_xticklabels(a.get_xticklabels(), rotation = 45, fontsize = 8)
# a.set_yticklabels(a.get_yticklabels(), rotation = 45, fontsize = 8)
# plt.title("DS3", weight = "bold")
# plt.xlabel("DIP Rate (doublings/hour)")
# plt.ylabel("Division Rate (per hour)")
# plt.tight_layout()
# # plt.savefig('DS3_paramscan_TR1_tiny.png')
# # plt.savefig('DS3_paramscan_TR1.png')
#
# df_DS4 = pd.read_pickle('PC9-DS4_param-scan_tighterRange.pkl')
# df_DS4 = df_DS4.round({'death rate': 5, 'division rate': 5, 'DIP rate': 5})
# df_DS4['cell.line'] = 'PC9.DS4'
# result_DS4 = df_DS4.pivot(index='division rate', columns='DIP rate', values='p-value')
# plt.figure('DS4', figsize=(3,4))
# a = sns.heatmap(result_DS4, cbar_kws={'label': 'p-value'}, cmap = 'jet')
# a.set_xticklabels(a.get_xticklabels(), rotation = 45, fontsize = 8)
# a.set_yticklabels(a.get_yticklabels(), rotation = 45, fontsize = 8)
# plt.title("DS4", weight = "bold")
# plt.xlabel("DIP Rate (doublings/hour)")
# plt.ylabel("Division Rate (per hour)")
# plt.tight_layout()
# # plt.savefig('DS4_paramscan_TR1_tiny.png')
# # plt.savefig('DS4_paramscan_TR1.png')
#
# df_DS7 = pd.read_pickle('PC9-DS7_param-scan_tighterRange.pkl')
# df_DS7 = df_DS7.round({'death rate': 5, 'division rate': 5, 'DIP rate': 5})
# df_DS7['cell.line'] = 'PC9.DS7'
# result_DS7 = df_DS7.pivot(index='division rate', columns='DIP rate', values='p-value')
# plt.figure('DS7', figsize=(3,4))
# a = sns.heatmap(result_DS7, cbar_kws={'label': 'p-value'}, cmap = 'jet')
# a.set_xticklabels(a.get_xticklabels(), rotation = 45, fontsize = 8)
# a.set_yticklabels(a.get_yticklabels(), rotation = 45, fontsize = 8)
# plt.title("DS7", weight = "bold")
# plt.xlabel("DIP Rate (doublings/hour)")
# plt.ylabel("Division Rate (per hour)")
# plt.tight_layout()
# # plt.savefig('DS7_paramscan_TR1_tiny.png')
# # plt.savefig('DS7_paramscan_TR1.png')
##########


##### DEVELOPMENT - FINDING DUPLICATED ROWS #####
# all_dfs = [df_DS1, df_DS3, df_DS4, df_DS6, df_DS7, df_DS8, df_DS9]
# df_all = pd.concat(all_dfs)
#
# def label_cell_line(df):
#     if df['p-value'] > 0.1:
#         return df['cell.line']
#     else:
#         return "not.assigned"
#
# df_all['cell.line.new'] = df_all.apply(lambda df_all: label_cell_line(df_all), axis = 1)
# # df_all.columns = ['DIP', 'death', 'division', 'pval', 'cellLine', 'cellLineNew']
# df_all.to_csv('all_cellLine_tile.csv')
# quit()
# print(df_all)
# duplicateRowsDF = df_all[df_all.duplicated(['DIP', 'division', 'cellLineNew'])]
#
# print("Duplicate Rows except first occurrence based on all columns are :")
# print(duplicateRowsDF)
# # quit()
# df_all_new = df_all.drop_duplicates(subset=['DIP', 'division', 'cellLineNew'], keep='first')
# # df_all_new = df_all_new.round({'death': 5, 'division': 5, 'DIP': 5})
# print(df_all_new)
# result_all = df_all_new.pivot(index="division", columns='DIP', values='cellLineNew')
# print(result_all)
# quit()
#####




##### OLD PHASE SCAN PLOTTING #####
# # fig = plt.figure(figsize = (15,9)) # width x height
# # ax1 = fig.add_subplot(2, 3, 1) # row, column, position
# # ax2 = fig.add_subplot(2, 3, 2)
# # ax3 = fig.add_subplot(2, 3, 3)
# # ax4 = fig.add_subplot(2, 3, 4)
# # ax5 = fig.add_subplot(2, 3, 5)
# # ax6 = fig.add_subplot(2, 3, 6)
# # # ax7 = fig.add_subplot(3, 3, 7)
# #
# # a = sns.heatmap(result_DS1, cbar_kws={'label': 'p-value'}, cmap = 'jet', ax = ax1)
# # a.set_xticklabels(a.get_xticklabels(), rotation = 45, fontsize = 8)
# # a.set_yticklabels(a.get_yticklabels(), rotation = 45, fontsize = 8)
# #
# # b = sns.heatmap(result_DS3, cbar_kws={'label': 'p-value'}, cmap = 'jet', ax = ax2)
# # b.set_xticklabels(a.get_xticklabels(), rotation = 45, fontsize = 8)
# # b.set_yticklabels(a.get_yticklabels(), rotation = 45, fontsize = 8)
# #
# # c = sns.heatmap(result_DS4, cbar_kws={'label': 'p-value'}, cmap = 'jet', ax = ax3)
# # c.set_xticklabels(a.get_xticklabels(), rotation = 45, fontsize = 8)
# # c.set_yticklabels(a.get_yticklabels(), rotation = 45, fontsize = 8)
# #
# # d = sns.heatmap(result_DS6, cbar_kws={'label': 'p-value'}, cmap = 'jet', ax = ax4)
# # d.set_xticklabels(a.get_xticklabels(), rotation = 45, fontsize = 8)
# # d.set_yticklabels(a.get_yticklabels(), rotation = 45, fontsize = 8)
# # +
# # e = sns.heatmap(result_DS7, cbar_kws={'label': 'p-value'}, cmap = 'jet', ax = ax5)
# # e.set_xticklabels(a.get_xticklabels(), rotation = 45, fontsize = 8)
# # e.set_yticklabels(a.get_yticklabels(), rotation = 45, fontsize = 8)
# #
# # f = sns.heatmap(result_DS9, cbar_kws={'label': 'p-value'}, cmap = 'jet', ax = ax6)
# # f.set_xticklabels(a.get_xticklabels(), rotation = 45, fontsize = 8)
# # f.set_yticklabels(a.get_yticklabels(), rotation = 45, fontsize = 8)
# #
# # plt.tight_layout()
# # plt.show()
# quit()
#####




##### PHASED OUT FUNCTIONS FOR DRAWING HEATMAPS #####
# def draw_heatmap(*args, **kwargs):
#     data = kwargs.pop('data')
#     d = data.pivot(index=args[1], columns=args[0], values=args[2])
#     z = sns.heatmap(d, cmap = current_cmap, **kwargs)
#     z.set_xticklabels(z.get_xticklabels(), rotation=45, fontsize=8)
#     z.set_yticklabels(z.get_yticklabels(), rotation=45, fontsize=8)
#     plt.title("Discrete Subline Parameter Scan", weight="bold")
#
# fg = sns.FacetGrid(data, col='cell.line')
# fg.map_dataframe(draw_heatmap, 'DIP rate', 'division rate', 'p-value')
# plt.show()
#
# quit()

# df = pd.read_pickle('PC9_param-scan.pkl')
# df['DIP rate'] = df['division rate'] - df['death rate']
# df = df.round({'death rate': 3, 'division rate': 3, 'DIP rate': 3})
# df.rename(columns = {"dataset": "cell.line"}, inplace = True)
# df['cell.line'] = [item[0] for item in df['cell.line']]
#
# def draw_heatmap(*args, **kwargs):
#     data = kwargs.pop('data')
#     d = data.pivot(index=args[1], columns=args[0], values=args[2])
#     z = sns.heatmap(d, cbar_kws={'label': 'p-value'}, cmap = current_cmap, **kwargs)
#     z.set_xticklabels(z.get_xticklabels(), rotation=45, fontsize=8)
#     z.set_yticklabels(z.get_yticklabels(), rotation=45, fontsize=8)
#     plt.title("Discrete Subline Parameter Scan", weight="bold")
#
# fg = sns.FacetGrid(df, col='cell.line')
# fg.map_dataframe(draw_heatmap, 'DIP rate', 'division rate', 'p-value')
# plt.savefig('allDS_paramscan.png')
# plt.show()
# quit()
#####



##### SECOND PHASE OF SCAN PLOTTING #####
# df_DS8 = df[df['dataset']=='PC9.DS8']
# result_DS8 = df_DS8.pivot(index='division rate', columns='DIP rate', values='p-value')
# # print(result_DS8)
# # quit()
# plt.figure('DS8')
# a = sns.heatmap(result_DS8, cbar_kws={'label': 'p-value'}, cmap = current_cmap)
# a.set_xticklabels(a.get_xticklabels(), rotation = 45, fontsize = 8)
# a.set_yticklabels(a.get_yticklabels(), rotation = 45, fontsize = 8)
# plt.title("DS8 Parameter Scan", weight = "bold")
# plt.savefig('DS8_paramscan_secondpass.pdf')
#
# df_DS6 = df[df['dataset']=='PC9.DS6']
# result_DS6 = df_DS6.pivot(index='division rate', columns='DIP rate', values='p-value')
# plt.figure('DS6')
# b = sns.heatmap(result_DS6, cbar_kws={'label': 'p-value'}, cmap = current_cmap)
# b.set_xticklabels(b.get_xticklabels(), rotation = 45, fontsize = 8)
# b.set_yticklabels(b.get_yticklabels(), rotation = 45, fontsize = 8)
# plt.title("DS6 Parameter Scan", weight = "bold")
# plt.savefig('DS6_paramscan_secondpass.pdf')
#
# df_DS9 = df[df['dataset']=='PC9.DS9']
# result_DS9 = df_DS9.pivot(index='division rate', columns='DIP rate', values='p-value')
# plt.figure('DS9')
# c = sns.heatmap(result_DS9, cbar_kws={'label': 'p-value'}, cmap = current_cmap)
# c.set_xticklabels(c.get_xticklabels(), rotation = 45, fontsize = 8)
# c.set_yticklabels(c.get_yticklabels(), rotation = 45, fontsize = 8)
# plt.title("DS9 Parameter Scan", weight = "bold")
# plt.savefig('DS9_paramscan_secondpass.pdf')
# plt.show()
# quit()
#####



#### OPENING DATA AND PLOTTING FRICK'S OLD cFP DATA #####
# cFP_rates = pd.read_csv("Frick_cFP_rates.csv")
# cell_lines = ['PC9', 'PC9.BR1', 'HCC4006', 'PC9.MGH', 'HCC827', 'II-18']
# cFP_rates_CL = cFP_rates[cFP_rates['cell.line'].isin(cell_lines)]
# discrete_sublines = ['PC9.DS9', 'PC9.DS6', 'PC9.DS8']
# cFP_rates_DS = cFP_rates[cFP_rates['Cell_Line'].isin(cell_lines)]
# sns.distplot(cFP_rates[(cFP_rates['cell.line'] == 'PC9') & (cFP_rates['drug1'] == 'erlotinib')]['rate'],
#              hist=False, kde_kws={"color": "b", "lw": 3, "label": "PC9-VU", "shade": True, "ls": "--"})
# sns.distplot(cFP_rates[(cFP_rates['cell.line'] == 'PC9.DS8') & (cFP_rates['drug1'] == 'erlotinib')]['rate'],
#              hist=False, kde_kws={"color": "y", "lw": 3, "label": "PC9-DS8", "shade": True})
# sns.distplot(cFP_rates[(cFP_rates['cell.line'] == 'PC9.DS9') & (cFP_rates['drug1'] == 'erlotinib')]['rate'],
#              hist=False, kde_kws={"color": "k", "lw": 3, "label": "PC9-DS9", "shade": True})
# sns.distplot(cFP_rates[(cFP_rates['cell.line'] == 'PC9.DS6') & (cFP_rates['drug1'] == 'erlotinib')]['rate'],
#              hist=False, kde_kws={"color": "c", "lw": 3, "label": "PC9-DS6", "shade": True})
# plt.axvline(x=0, c='k', ls='--', lw=3)
# plt.ylim(0,150)
# plt.xlim(-0.065, 0.045)
# plt.legend()
# plt.xlabel("Proliferation Rate")
# plt.ylabel("Density")
# plt.title("Single-cell Rate Distribution", weight = "bold")
# plt.savefig('PC9-SubLines_DIPdist_kde_DS8.pdf')
# plt.show()
# quit()
# dsPlot = sns.FacetGrid(cFP_rates_DS, col = "Cell_Line", hue = 'Cell_Line',
#                        col_order=cell_lines)
# dsPlot = (dsPlot.map(sns.distplot, 'DIP_Rate', hist=False, kde_kws={"shade":True, "lw":3}))
# plt.savefig("all_DSLines.pdf")
# quit()

# g = sns.FacetGrid(cFP_rates_CL, col = "cell.line", hue = 'drug1', col_wrap=3)
# g = (g.map(sns.distplot, 'rate'))
# # g.set(ylabel='Frequency')
# plt.ylim(0,150)
# plt.legend(loc='lower right', bbox_to_anchor=(1,0.75))
# plt.savefig('Fig1a_firstpass.jpeg')
# plt.show()
# quit()

# cFP_PC9_DMSO = cFP_rates[(cFP_rates['cell.line'] == 'PC9') & (cFP_rates['drug1'] == 'control')]
# cFP_PC9_Erl = cFP_rates[(cFP_rates['cell.line'] == 'PC9') & (cFP_rates['drug1'] == 'erlotinib')]

# cFP_DS9_DMSO = cFP_rates[(cFP_rates['cell.line'] == 'PC9.DS9') & (cFP_rates['drug1'] == 'control')]
# cFP_DS9_Erl = cFP_rates[(cFP_rates['cell.line'] == 'PC9.DS9') & (cFP_rates['drug1'] == 'erlotinib')]

# cFP_DS8_DMSO = cFP_rates[(cFP_rates['cell.line'] == 'PC9.DS8') & (cFP_rates['drug1'] == 'control')]
# cFP_DS8_Erl = cFP_rates[(cFP_rates['cell.line'] == 'PC9.DS8') & (cFP_rates['drug1'] == 'erlotinib')]
# print(cFP_DS8_Erl['cell.line'].unique())
# quit()
# cFP_DS6_DMSO = cFP_rates[(cFP_rates['cell.line'] == 'PC9.DS6') & (cFP_rates['drug1'] == 'control')]
# cFP_DS6_Erl = cFP_rates[(cFP_rates['cell.line'] == 'PC9.DS6') & (cFP_rates['drug1'] == 'erlotinib')]
#####



cFP_rates = pd.read_csv("cFP_rates_VUlines.csv")
cell_lines = ['PC9-DS1', 'PC9-DS3', 'PC9-DS4', 'PC9-DS7', 'PC9-DS6', 'PC9-DS9', 'PC9-DS8']



##### PHASE ONE OF PLOTTING MODEL TRAJECTORIES #####
# num_cells = 1
# div = 0.02 * np.log(2)
# dth = 0.03 * np.log(2)
# dat = cFP_DS8_Erl
# Model()
# Monomer('Cell')
# Parameter('Cell_init', num_cells)
# Initial(Cell, Cell_init)
# Observable('Obs_Cell', Cell())
# Parameter('k_div', div)
# Parameter('k_dth', dth)
# Rule('Division', Cell() >> Cell() + Cell(), k_div)
# Rule('Death', Cell() >> None, k_dth)
# sim_dist = []
# plt.figure()
# for rep in range(len(dat)):
#     t1 = np.linspace(0, 144, 145) #hrs
#     sim1 = BngSimulator(model, tspan=t1, verbose=False)  # , seed = 1095205711)
#     x1 = sim1.run(tspan=t1, param_values={'k_div': 0.04 * np.log(2),
#                                          'k_dth': 0.005 * np.log(2)},
#                  verbose=False)
#     cell_tot = x1.observables["Obs_Cell"][-1]
#     t2 = np.linspace(0, 144, 145)  # 6 days in drug
#     x2 = sim1.run(tspan=t2, param_values={'k_div': model.parameters['k_div'].value,
#                                          "k_dth": model.parameters['k_dth'].value,
#                                          "Cell_init": cell_tot},
#                  n_sim=1, verbose=False)
#     if cell_tot != 0:
#         slope, intercept, r_value, p_value, std_err = sp.stats.linregress(t2, np.log2(
#             x2.observables["Obs_Cell"] / x2.observables["Obs_Cell"][0]))
#         if math.isnan(slope) == False:
#             sim_dist = sim_dist + [slope]
# sns.distplot(cFP_DS8_Erl['rate'], kde = True, color="r")
# sns.distplot(sim_dist, kde=True, color = "g")
# D_stat, p_val = sp.ks_2samp(dat['rate'], sim_dist)
# print(p_val)
# # print "Normal", slope
# #     plt.plot(t1, np.log2(x1.observables["Obs_Cell"]), 'g', lw=2)
# #     plt.plot(t1[-1] + x2.tout[0], np.log2(x2.observables["Obs_Cell"]), 'r', lw = 2)
# # plt.axvline(x=144, linewidth=4, color='k')
# # plt.xlabel("Time (hours)", weight = "bold")
# # plt.ylabel("log2(cell count)", weight = "bold")
# plt.show()
# quit()
#####



##### PHASE TWO OF MODEL TRAJECTORIES PLOTTING #####
def dist_compare(samp, div, dth):
    dat = cFP_rates[cFP_rates['Cell_Line'] == samp]

    kdiv = div
    kdth = dth
    num_cells = 1
    Model()
    Monomer('Cell')
    Parameter('Cell_init', num_cells)
    Initial(Cell, Cell_init)
    Observable('Obs_Cell', Cell())
    Parameter('k_div', kdiv)
    Parameter('k_dth', kdth)
    Rule('Division', Cell() >> Cell() + Cell(), k_div)
    Rule('Death', Cell() >> None, k_dth)

    sim_dist = []
    # plt.figure(figsize=[4,3])
    t1 = np.linspace(0, 200, 201) #hrs
    sim1 = BngSimulator(model, tspan=t1, verbose=False)  # , seed = 1095205711)
    x1 = sim1.run(tspan=t1, param_values={'k_div': 0.04 * np.log(2),
                                          'k_dth': 0.005 * np.log(2)},
                  n_runs=len(dat), verbose=False)

    trajs = np.array(np.array([tr[:]["Obs_Cell"] for tr in np.array(x1.observables)]).T)
    # print(x1.observables[:]['Obs_Cell'])
    # print(np.array(np.array([tr[:]["Obs_Cell"] for tr in np.array(x1.observables)]).T))
    # plt.plot(t1, np.log2(x1.observables["Cell"]), 'red', lw=1)
    # plt.plot(t1, np.log2(trajs), 'red', lw=1)


    # cell_tot = trajs[-1]
    cell_pop = trajs[-1]
    print(cell_pop)
    trajects = []
    for ind, cell in enumerate(cell_pop):
        print("%d" % ind)
        print("%d" % cell)
        # cell_pop_1n = cell_pop1[ind]
        t2 = np.linspace(0, 225, 226)  # in drug
        x2 = sim1.run(tspan=t2,
                      initials={model.species[0]: cell},
                      verbose=False)
        # trajects.append(np.log2(x2.observables["Obs_Cell"] / x2.observables["Obs_Cell"][0]))
        if cell != 0:
            slope, intercept, r_value, p_value, std_err = sp.stats.linregress(t2, np.log2(
                x2.observables["Obs_Cell"] / x2.observables["Obs_Cell"][0]))
            if math.isnan(slope) == False:
                sim_dist = sim_dist + [slope]
        if cell > 50:
            trajects.append(np.log2(x2.observables["Obs_Cell"] / x2.observables["Obs_Cell"][0]))

        # plt.plot(t1[-1] + x2.tout[0],
        #          np.log2(x2.observables["Obs_Cell"] / x2.observables["Obs_Cell"][0]),
        #          color = 'grey', lw=0.5)



    # plt.xlabel("Time (hours)")
    # plt.ylabel("Norm Log2 Cell Count")
    # plt.title(r"Post-drug Dynamics $(k_{division}=%.3f, k_{death}=%.5f)$" % (kdiv,kdth), fontsize = 10)
    # plt.savefig("PreDrug_withMetrics_DS7_normalized.svg")

    D_stat, p_val = sp.ks_2samp(dat['DIP_Rate'], sim_dist)

    return trajects, sim_dist, p_val

# data_DS7 = dist_compare('PC9-DS7', 0.028, 0.02625)
# trajectories_DS7 = pd.DataFrame(data_DS7[0])
# trajectories_DS7 = trajectories_DS7.transpose()
# trajectories_DS7.to_csv('trajectories_DS7.csv')
#
# distributions_DS7 = pd.DataFrame({'DS7': data_DS7[1]})
# distributions_DS7.to_csv('distributions_DS7.csv')
#
# print(data_DS7[2])
# quit()



DSs = ['PC9-DS1', 'PC9-DS3', 'PC9-DS4', 'PC9-DS6', 'PC9-DS7', 'PC9-DS9']
divs_a = [0.032, 0.030, 0.035, 0.030, 0.028, 0.025]
dths_a = [0.03050, 0.03075, 0.03175, 0.02940, 0.02625, 0.02375]
data_a = []
for i in range(len(DSs)):
    data_a.append(dist_compare(DSs[i], divs_a[i], dths_a[i]))

trajectories_DS1 = pd.DataFrame(data_a[0][0])
trajectories_DS1 = trajectories_DS1.transpose()
trajectories_DS1.to_csv('trajectories_DS1_G50.csv')

trajectories_DS3 = pd.DataFrame(data_a[1][0])
trajectories_DS3 = trajectories_DS3.transpose()
trajectories_DS3.to_csv('trajectories_DS3_G50.csv')

trajectories_DS4 = pd.DataFrame(data_a[2][0])
trajectories_DS4 = trajectories_DS4.transpose()
trajectories_DS4.to_csv('trajectories_DS4_G50.csv')

trajectories_DS6 = pd.DataFrame(data_a[3][0])
trajectories_DS6 = trajectories_DS6.transpose()
trajectories_DS6.to_csv('trajectories_DS6_G50.csv')

trajectories_DS7 = pd.DataFrame(data_a[4][0])
trajectories_DS7 = trajectories_DS7.transpose()
trajectories_DS7.to_csv('trajectories_DS7_G50.csv')

trajectories_DS9 = pd.DataFrame(data_a[5][0])
trajectories_DS9 = trajectories_DS9.transpose()
trajectories_DS9.to_csv('trajectories_DS9_G50.csv')

# print(trajectories_DS1)

distributions_0 = pd.DataFrame({'DS1': data_a[0][1]})
distributions_1 = pd.DataFrame({'DS3': data_a[1][1]})
distributions_2 = pd.DataFrame({'DS4': data_a[2][1]})
distributions_3 = pd.DataFrame({'DS6': data_a[3][1]})
distributions_4 = pd.DataFrame({'DS7': data_a[4][1]})
distributions_5 = pd.DataFrame({'DS9': data_a[5][1]})
distributions = pd.concat([distributions_0, distributions_1,
                           distributions_2, distributions_3,
                           distributions_4, distributions_5],
                          ignore_index=True, axis = 1)
distributions.columns = ['DS1', 'DS3', 'DS4', 'DS6', 'DS7', 'DS9']

print(distributions)
distributions.to_csv('distributions_G50.csv')

pvalues = pd.DataFrame({'DS1': data_a[0][2], 'DS3': data_a[1][2],
                       'DS4': data_a[2][2], 'DS6': data_a[3][2],
                       'DS7': data_a[4][2], 'DS9': data_a[5][2]}, index=[0])
print(pvalues)
pvalues.to_csv('pvalues_G50.csv')

quit()






mean_exp = np.mean(dat['DIP_Rate'])
sd_exp = np.std(dat['DIP_Rate'])
mean_sim = np.mean(sim_dist)
sd_sim = np.std(sim_dist)

plt.figure(figsize=[4,3])
sns.distplot(sim_dist, kde = True, hist= False, color = "grey", kde_kws={"shade": True})
sns.distplot(dat['DIP_Rate'], kde = True, hist= False, color = "darkorchid", kde_kws={"shade": True})
plt.xlabel("DIP Rate")
plt.ylabel("Density")
plt.xlim(-0.025, 0.025)
plt.ylim(0,200)
plt.legend(labels=['Simulated','Experimental'], loc = 'upper left')
# plt.text(0.007, 200, "Mean (Exp)=%.3f" % round(mean_exp, 3), fontsize = 8)
# plt.text(0.007, 180, "SD (Exp)=%.4f" % round(sd_exp, 4), fontsize = 8)
# plt.text(0.007, 160, "Mean (Sim)=%.3f" % round(mean_sim, 3), fontsize = 8)
# plt.text(0.007, 140, "SD (Sim)=%.4f" % round(sd_sim, 4), fontsize = 8)
plt.text(0.007, 150, "p=%.3f" % round(p_val, 3), fontsize = 12)
plt.title("DIP Rate Distribution")
plt.savefig("DIP_norm_DS7.svg")

plt.show()
quit()


##### MODEL GENERALIZED FOR MULTIPLE STATES ####
# kdiv = [0.026, 0.036]
# kdth = [0.025, 0.03]
# Model()
# Monomer('Cell', ['type'], {'type': ["x"+str(i) for i in range(len(kdiv))]})
# [Initial(Cell(type="x"+str(i)), Parameter('Cell%d_Init' % i, num_cells)) for i in range(len(kdiv))]
# [Rule('divide_Cell%d' % i, Cell(type="x"+str(i)) >> Cell(type="x"+str(i)) + Cell(type="x"+str(i)),
#       Parameter('kdiv_%d' % i, kdiv[i])) for i in range(len(kdiv))]
# [Rule('die_Cell%d' % i, Cell(type="x"+str(i)) >> None,
#       Parameter('kdth_%d' % i, kdth[i])) for i in range(len(kdiv))]
# Observable("Cell_total", Cell())
# [Observable("Cell_t_%s" % i, Cell(type="x"+str(i))) for i in range(len(kdiv))]
#####


for ind, cell in enumerate(cell_tot):
    print("%d" % ind)

    # cell_pop1 = np.random.multinomial(int(round(cell)), [1,0])
    # cell_pop2 = np.random.multinomial(int(round(cell)), [0,1])
    # print(cell_pop1)
    t2 = np.linspace(0, 225, 226)  # in drug
    x2 = sim1.run(tspan=t2, n_runs=cell_tot_s[0],
                  initials={model.species[i]: pop for i, pop in enumerate(cell_pop1)},
                  verbose=False)
    # print(type(cell_pop))
    # print(enumerate(cell_pop))
    # print([e for e in enumerate(cell_pop)])
    # print({model.species[i]:pop for i,pop in enumerate(cell_pop)})
    # # quit()
    # print(model.species)
    # print(len(model.species))
    # print(model.parameters)
    # print(model.initial_conditions)
    # for i, pop in enumerate(cell_pop):
    #     print(model.species[i])
    #     print(pop)
    # quit()
    t2 = np.linspace(0, 225, 226)  #  in drug
    x2 = sim1.run(tspan=t2, n_runs=cell_tot_s[0],
                  initials={model.species[i]: pop for i, pop in enumerate(cell_pop1)},
                  verbose=False)
    x3 = sim1.run(tspan=t2, n_runs=cell_tot_s[1],
                  initials={model.species[i]: pop for i, pop in enumerate(cell_pop2)},
                  verbose=False)
    # x2 = sim1.run(tspan=t2, n_runs=1,
    #               initials={model.species[i]:pop for i,pop in enumerate(cell_pop)},
    #               verbose=False)
    # x2 = sim1.run(tspan=t2, param_values={'k_div': model.parameters['k_div'].value,
    #                                   "k_dth": model.parameters['k_dth'].value
    #                                   "Cell_init": cell_tot},
    #          n_sim=1, verbose=False)

    plt.plot(t1[-1] + x3.tout[0], np.log2(x3.observables["Cell_total"]), color = 'green', lw=1)

    #
    # if cell != 0:
    #     slope, intercept, r_value, p_value, std_err = sp.stats.linregress(t2, np.log2(
    #         x2.observables["Cell_total"] / x2.observables["Cell_total"][0]))
    #     if math.isnan(slope) == False:
    #         sim_dist = sim_dist + [slope]
    #         # print "Normal", slope

plt.show()
quit()

D_stat, p_val = sp.ks_2samp(dat['DIP_Rate'], sim_dist)
mean_exp = np.mean(dat['DIP_Rate'])
sd_exp = np.std(dat['DIP_Rate'])
mean_sim = np.mean(sim_dist)
sd_sim = np.std(sim_dist)

    # return p_val

# print(cell_tot)
# print(len(dat))
plt.xlabel("Time (hours)")
plt.ylabel("Normalized Log2 Count")

# a = round(np.log(2) * 0.04, 3)
# b = round(np.log(2) * 0.005, 3)
# plt.title(r"Post-drug Dynamics $(k_{division}=%.3f, k_{death}=%.5f)$" % (div,dth), fontsize = 10)
# plt.savefig("PreDrug_withMetrics_DS7bad.svg")



# cell_tot_all_log = np.log10(cell_tot_all)

# plt.figure(figsize=[4,3])
# sns.distplot(cell_tot_all_log, kde = True, hist= False,  color = "grey", kde_kws={"shade": True})
# plt.xlabel("Log 10 Count")
# plt.ylabel("Density")
# plt.title('Pre-drug Cell Count')
# plt.savefig("PreDrugCount_withTitle.svg")
#
# plt.figure()
# colors = cm.rainbow(np.linspace(0, 1, len(cell_tot)))
# for y, c in zip(cell_tot, colors):
#     plt.plot(t1[-1] + x2.tout[0], np.log2(x2.observables["Obs_Cell"]), color = c, lw = 2)
# plt.xlabel("Time (hours)")
# plt.ylabel("Normalized Log2 Count")
#
plt.figure(figsize=[4,3])
sns.distplot(sim_dist, kde = True, hist= False, color = "grey", kde_kws={"shade": True})
sns.distplot(dat['DIP_Rate'], kde = True, hist= False, color = "seagreen", kde_kws={"shade": True})
plt.xlabel("DIP Rate")
plt.ylabel("Density")
plt.xlim(-0.025, 0.025)
plt.ylim(0,200)
plt.legend(labels=['Simulated','Experimental'], loc = 'upper left')
# plt.text(0.007, 200, "Mean (Exp)=%.3f" % round(mean_exp, 3), fontsize = 8)
# plt.text(0.007, 180, "SD (Exp)=%.4f" % round(sd_exp, 4), fontsize = 8)
# plt.text(0.007, 160, "Mean (Sim)=%.3f" % round(mean_sim, 3), fontsize = 8)
# plt.text(0.007, 140, "SD (Sim)=%.4f" % round(sd_sim, 4), fontsize = 8)
plt.text(0.007, 150, "p=%.3f" % round(p_val, 3), fontsize = 12)
plt.title("DIP Rate Distribution")
# plt.savefig("DIP_DS7bad.svg")
plt.show()
quit()
#
# divs_DS9 = np.linspace(0, 0.07, 24) * np.log(2)
# dips_DS9 = np.linspace(-0.03, -0.01, 9) * np.log(2)
# p_vals_DS9 = []
# div_rates_DS9 = []
# dth_rates_DS9 = []
# dip_rates_DS9 = []
# for dip in dips_DS9:
#     for di in divs_DS9:
#         dt = di-dip
#         # print(dip, di, dt)
#         p = dist_compare(cFP_DS9_Erl, 1, di, dt)
#         p_vals_DS9.append(p)
#         dip_rates_DS9.append(dip)
#         div_rates_DS9.append(di)
#         dth_rates_DS9.append(dt)
#         print("run done")
#
# dict = {'DIP rate': dip_rates_DS9, 'division rate': div_rates_DS9,
#         'death rate': dth_rates_DS9, 'p-value': p_vals_DS9}
#
# df = pd.DataFrame(data=dict)
#
# df.to_pickle('PC9-DS9_param-scan.pkl')
#
# quit()
### Doing first pass (non DIP targeted) parameter scan ###
# dats = [cFP_DS8_Erl, cFP_DS6_Erl, cFP_DS9_Erl]
# divs = np.linspace(0, 0.1, 11) * np.log(2) # 31) #[.01 .04]
# dths = np.linspace(0, 0.1, 11) * np.log(2) # 21) #[0.04, 0.06]
#
# # print(divs)
# # print(dths)
# # quit()
# p_vals = []
# div_rates = []
# dth_rates = []
# datasets = []
# count = 0
# for da in dats:
#     for di in divs:
#         for dt in dths:
#             p = dist_compare(da, 1, di, dt)
#             p_vals.append(p)
#             div_rates.append(di)
#             dth_rates.append(dt)
#             datasets.append(da)
#             # count = count + 1
#             print("run done" + count)
#
#
# dict = {'dataset': datasets, 'division rate': div_rates,
#         'death rate': dth_rates, 'p-value': p_vals}
#
# df = pd.DataFrame(data=dict)
#
# df.to_pickle('PC9_param-scan.pkl')

# print(p_vals)
# range of numbers multiplied by np.log(2)
# return pvalues (p_vals.append(p_val))
# do for subset of one DS, then whole DS, then all DSs

# print(dist_compare(cFP_DS8_Erl, 1, 0.02*np.log(2), 0.03*np.log(2)))

# quit()
#
# all_data = list(cFP_PC9_Erl['rate']) + list(cFP_DS9_Erl['rate']) + list(cFP_DS8_Erl['rate']) + list(cFP_DS6_Erl['rate'])
# min_lim = -0.05
# max_lim = 0.05
#
# labels = ['cFP PC9','cFP DS9', 'cFP DS8', 'cFP DS6']
# bins = np.arange(min(all_data), max(all_data) + 0.001, 0.001)
# plt.figure()
# sns.distplot(cFP_PC9_Erl['rate'], kde = True, color="k", bins = bins)
# sns.distplot(cFP_DS9_Erl['rate'], kde = True, color="b", bins = bins)
# sns.distplot(cFP_DS8_Erl['rate'], kde = True, color="r", bins = bins)
# sns.distplot(cFP_DS6_Erl['rate'], kde = True, color="g", bins = bins)
# plt.xlim(-0.075,0.075)
# plt.ylim(0,200)
# plt.xlabel("DIP Rate", weight = "bold")
# plt.ylabel("Density", weight = "bold")
# plt.title("PC9 Variants", weight = "bold")
# plt.legend(title = "Experiments", labels = labels, loc = "upper left")
# plt.savefig('PC9_variants_cFP_Erl.pdf')
# # plt.show()
# quit()

# div_rates = list(np.linspace(0,0.05, 51))
# death_rates = list(np.linspace(0, 0.1, 101))
#
# for div in div_rates:
#     for death in death_rates:
