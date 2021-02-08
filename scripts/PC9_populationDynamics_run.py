from pysb.bng import generate_equations
from pysb.core import *
from pysb.simulator.bng import BngSimulator
from pysb.integrate import odesolve
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
# from ggplot import *
import scipy.stats as ss
from scipy.stats import sem
from collections import OrderedDict

sns.set(font_scale = 1.25)
sns.set_style("whitegrid")

def plot_mean_min_max(name, title=None, color="red"):
    trajectories = y1
    x = np.array([tr[:][name] for tr in trajectories]).T
    tout = x1.tout.T
    # plt.figure()
    plt.title(title)
    plt.subplot(111)
    mean = np.log2(x.mean(1) / x.mean(1)[0])
    ci_95 = 1.96 * sem(np.log2(x / x[0]).T)
    ci_min = mean - ci_95
    ci_max = mean + ci_95
    # print(mean)
    # print(ci_95)
    # print(mean.shape)
    # print(np.log2(x/x[0]))
    # print(ci_95)
    # print(ci_95.shape)
    # print(tout.shape)
    # print(ci_min.shape)
    # print(ci_max.shape)
    # min = np.log2(x.min(1)/x.min(1)[0])
    # max = np.log2(x.max(1)/x.max(1)[0])
    plt.fill_between(t, ci_min, ci_max, facecolor=color, alpha=0.25)
    plt.plot(tout, np.log2(x / x[0]), '0.5', color=color, lw=1, alpha=0.25)  # individual trajectories
    # plt.plot(tout, np.log2(x.mean(1)/x.mean(1)[0]), 'k-*', lw=3)#, label="Mean")
    # plt.plot(tout, np.log2(x.min(1)/x.min(1)[0]), 'b--', lw=3)#, label="Minimum")
    # plt.plot(tout, np.log2(x.max(1)/x.max(1)[0]), 'r--', lw=3)#, label="Maximum")
    # if model is not None:
    #     y = odesolve(model, tspan = t)
    #     plt.plot(tout, y[title], 'g--', lw=3, label="ODE")
    # plt.tight_layout()
    # if savename is not None:
    #     plt.savefig("%s.png" % savename)
    # if show:
    #     plt.show()
    # plt.close()

import PC9_populationDynamics_model as popD

num_cells = 5000
pre_drug_time = 25
post_drug_time = 750

Model()

popD.declare_monomers()
popD.declare_parameters()
popD.declare_initial_conditions()
popD.declare_observables()
popD.declare_functions()


generate_equations(model, verbose=True)
from pysb.generator.bng import BngGenerator
print BngGenerator(model).get_content()

n_sims = 6

# slopes1 = []
# for sim in range(n_sims):
#
#     t0 = np.linspace(0,pre_drug_time,pre_drug_time + 1)
#     sim0 = BngSimulator(model, tspan=t0, verbose=False)
#     x0 = sim0.run(tspan=t0, verbose=False, cleanup=False,
#                  param_values={"k_dip_PC9_VU": 0.023 * np.log(2)}, method='ssa')
#
#     PC9_VU_0.value = x0.all["Obs_PC9_VU"][-1]
#     PC9_BR1_0.value = x0.all["Obs_PC9_BR1"][-1]
#
#     t = np.linspace(0, post_drug_time, post_drug_time + 1)
#     sim = BngSimulator(model, tspan=t, verbose=False)
#     x1 = sim.run(tspan=t, verbose=False, cleanup=False,
#                  param_values={"k_dip_PC9_VU": model.parameters['k_dip_PC9_VU'].value},
#                  method='ssa')
#
#     slope, intercept, r_value, p_value, std_err = \
#         ss.linregress(x1.tout[0], np.log2(x1.all["Obs_All"] / x1.all["Obs_All"][0]))
#
#     slopes1.append(slope)
#
#     plt.plot(x0.tout[0], np.log2(x0.all["Obs_All"] / x0.all["Obs_All"][0]), 'pink', lw=2)
#     plt.plot(x0.tout[0][-1] + x1.tout[0],
#              np.log2(x1.all["Obs_All"] / x0.all["Obs_All"][0]),
#              'pink', lw=2, label = "PC9-VU")
#
# ###
#
# model.parameters["PC9_VU_0"].value = num_cells * 0.999
# model.parameters["PC9_BR1_0"].value = num_cells * 0.001
#
# slopes2 = []
# for sim in range(n_sims):
#     t0 = np.linspace(0, pre_drug_time, pre_drug_time + 1)
#     sim0 = BngSimulator(model, tspan=t0, verbose=False)
#     x0 = sim0.run(tspan=t0, verbose=False, cleanup=False,
#                   param_values={"k_dip_PC9_VU": 0.023 * np.log(2)}, method='ssa')
#
#     PC9_VU_0.value = x0.all["Obs_PC9_VU"][-1]
#     PC9_BR1_0.value = x0.all["Obs_PC9_BR1"][-1]
#
#     t = np.linspace(0, post_drug_time, post_drug_time + 1)
#     sim = BngSimulator(model, tspan=t, verbose=False)
#     x1 = sim.run(tspan=t, verbose=False, cleanup=False,
#                  param_values={"k_dip_PC9_VU": model.parameters['k_dip_PC9_VU'].value},
#                  method='ssa')
#
#     slope, intercept, r_value, p_value, std_err = \
#         ss.linregress(x1.tout[0], np.log2(x1.all["Obs_All"] / x1.all["Obs_All"][0]))
#
#     slopes2.append(slope)
#
#     plt.plot(x0.tout[0], np.log2(x0.all["Obs_All"] / x0.all["Obs_All"][0]), 'gold', lw=2)
#     plt.plot(x0.tout[0][-1] + x1.tout[0],
#              np.log2(x1.all["Obs_All"] / x0.all["Obs_All"][0]),
#              'gold', lw=2, label = "0.1% PC9-BR1")

###

model.parameters["PC9_VU_0"].value = num_cells * 0.99
model.parameters["PC9_BR1_0"].value = num_cells * 0.01

slopes3 = []
for sim in range(n_sims):
    t0 = np.linspace(0, pre_drug_time, pre_drug_time + 1)
    sim0 = BngSimulator(model, tspan=t0, verbose=False)
    x0 = sim0.run(tspan=t0, verbose=False, cleanup=False,
                  param_values={"k_dip_PC9_VU": 0.023 * np.log(2)}, method='ssa')

    PC9_VU_0.value = x0.all["Obs_PC9_VU"][-1]
    PC9_BR1_0.value = x0.all["Obs_PC9_BR1"][-1]

    t = np.linspace(0, post_drug_time, post_drug_time + 1)
    sim = BngSimulator(model, tspan=t, verbose=False)
    x1 = sim.run(tspan=t, verbose=False, cleanup=False,
                 param_values={"k_dip_PC9_VU": model.parameters['k_dip_PC9_VU'].value},
                 method='ssa')

    slope, intercept, r_value, p_value, std_err = \
        ss.linregress(x1.tout[0], np.log2(x1.all["Obs_All"] / x1.all["Obs_All"][0]))

    slopes3.append(slope)
    plt.plot(x0.tout[0], np.log2(x0.all["Obs_PC9_VU"]), 'r', lw=2)
    plt.plot(x0.tout[0][-1] + x1.tout[0],
             np.log2(x1.all["Obs_PC9_VU"]),
             'r', lw=2, label = "PC9-VU Cells")
    plt.plot(x0.tout[0], np.log2(x0.all["Obs_PC9_BR1"]), 'b', lw=2)
    plt.plot(x0.tout[0][-1] + x1.tout[0],
             np.log2(x1.all["Obs_PC9_BR1"]),
             'b', lw=2, label = "PC9-BR1 Cells")
    plt.plot(x0.tout[0], np.log2(x0.all["Obs_All"]), 'k', lw=4)
    plt.plot(x0.tout[0][-1] + x1.tout[0],
             np.log2(x1.all["Obs_All"]),
             'k', lw=4, label = "Combined")

handles, labels = plt.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys(), loc=0, title="Cell Line")
# plt.legend()
plt.xlabel("Time")
plt.ylabel("log2(Cell Count)")
plt.title("1% Spike-In Simulation", weight = "bold")
plt.show()
quit()

###

model.parameters["PC9_VU_0"].value = num_cells * 0.95
model.parameters["PC9_BR1_0"].value = num_cells * 0.05

slopes4 = []
for sim in range(n_sims):
    t0 = np.linspace(0, pre_drug_time, pre_drug_time + 1)
    sim0 = BngSimulator(model, tspan=t0, verbose=False)
    x0 = sim0.run(tspan=t0, verbose=False, cleanup=False,
                  param_values={"k_dip_PC9_VU": 0.023 * np.log(2)}, method='ssa')

    PC9_VU_0.value = x0.all["Obs_PC9_VU"][-1]
    PC9_BR1_0.value = x0.all["Obs_PC9_BR1"][-1]

    t = np.linspace(0, post_drug_time, post_drug_time + 1)
    sim = BngSimulator(model, tspan=t, verbose=False)
    x1 = sim.run(tspan=t, verbose=False, cleanup=False,
                 param_values={"k_dip_PC9_VU": model.parameters['k_dip_PC9_VU'].value},
                 method='ssa')

    slope, intercept, r_value, p_value, std_err = \
        ss.linregress(x1.tout[0], np.log2(x1.all["Obs_All"] / x1.all["Obs_All"][0]))

    slopes4.append(slope)
    plt.plot(x0.tout[0], np.log2(x0.all["Obs_All"] / x0.all["Obs_All"][0]), 'green', lw=2)
    plt.plot(x0.tout[0][-1] + x1.tout[0],
             np.log2(x1.all["Obs_All"] / x0.all["Obs_All"][0]),
             'green', lw=2, label = "5% PC9-BR1")

###

model.parameters["PC9_VU_0"].value = num_cells * 0.90
model.parameters["PC9_BR1_0"].value = num_cells * 0.10

slopes5 = []
for sim in range(n_sims):
    t0 = np.linspace(0, pre_drug_time, pre_drug_time + 1)
    sim0 = BngSimulator(model, tspan=t0, verbose=False)
    x0 = sim0.run(tspan=t0, verbose=False, cleanup=False,
                  param_values={"k_dip_PC9_VU": 0.023 * np.log(2)}, method='ssa')

    PC9_VU_0.value = x0.all["Obs_PC9_VU"][-1]
    PC9_BR1_0.value = x0.all["Obs_PC9_BR1"][-1]

    t = np.linspace(0, post_drug_time, post_drug_time + 1)
    sim = BngSimulator(model, tspan=t, verbose=False)
    x1 = sim.run(tspan=t, verbose=False, cleanup=False,
                 param_values={"k_dip_PC9_VU": model.parameters['k_dip_PC9_VU'].value},
                 method='ssa')

    slope, intercept, r_value, p_value, std_err = \
        ss.linregress(x1.tout[0], np.log2(x1.all["Obs_All"] / x1.all["Obs_All"][0]))

    slopes5.append(slope)
    plt.plot(x0.tout[0], np.log2(x0.all["Obs_All"] / x0.all["Obs_All"][0]), 'blue', lw=2)
    plt.plot(x0.tout[0][-1] + x1.tout[0],
             np.log2(x1.all["Obs_All"] / x0.all["Obs_All"][0]),
             'blue', lw=2, label = "10% PC9-BR1")

###

model.parameters["PC9_VU_0"].value = num_cells * 0.75
model.parameters["PC9_BR1_0"].value = num_cells * 0.25

slopes6 = []
for sim in range(n_sims):
    t0 = np.linspace(0, pre_drug_time, pre_drug_time + 1)
    sim0 = BngSimulator(model, tspan=t0, verbose=False)
    x0 = sim0.run(tspan=t0, verbose=False, cleanup=False,
                  param_values={"k_dip_PC9_VU": 0.023 * np.log(2)}, method='ssa')

    PC9_VU_0.value = x0.all["Obs_PC9_VU"][-1]
    PC9_BR1_0.value = x0.all["Obs_PC9_BR1"][-1]

    t = np.linspace(0, post_drug_time, post_drug_time + 1)
    sim = BngSimulator(model, tspan=t, verbose=False)
    x1 = sim.run(tspan=t, verbose=False, cleanup=False,
                 param_values={"k_dip_PC9_VU": model.parameters['k_dip_PC9_VU'].value},
                 method='ssa')

    slope, intercept, r_value, p_value, std_err = \
        ss.linregress(x1.tout[0], np.log2(x1.all["Obs_All"] / x1.all["Obs_All"][0]))

    slopes6.append(slope)
    plt.plot(x0.tout[0], np.log2(x0.all["Obs_All"] / x0.all["Obs_All"][0]), 'magenta', lw=2)
    plt.plot(x0.tout[0][-1] + x1.tout[0],
             np.log2(x1.all["Obs_All"] / x0.all["Obs_All"][0]),
             'magenta', lw=2, label = "25% PC9-BR1")

###

model.parameters["PC9_VU_0"].value = num_cells * 0.50
model.parameters["PC9_BR1_0"].value = num_cells * 0.50

slopes7 = []
for sim in range(n_sims):
    t0 = np.linspace(0, pre_drug_time, pre_drug_time + 1)
    sim0 = BngSimulator(model, tspan=t0, verbose=False)
    x0 = sim0.run(tspan=t0, verbose=False, cleanup=False,
                  param_values={"k_dip_PC9_VU": 0.023 * np.log(2)}, method='ssa')

    PC9_VU_0.value = x0.all["Obs_PC9_VU"][-1]
    PC9_BR1_0.value = x0.all["Obs_PC9_BR1"][-1]

    t = np.linspace(0, post_drug_time, post_drug_time + 1)
    sim = BngSimulator(model, tspan=t, verbose=False)
    x1 = sim.run(tspan=t, verbose=False, cleanup=False,
                 param_values={"k_dip_PC9_VU": model.parameters['k_dip_PC9_VU'].value},
                 method='ssa')

    slope, intercept, r_value, p_value, std_err = \
        ss.linregress(x1.tout[0], np.log2(x1.all["Obs_All"] / x1.all["Obs_All"][0]))

    slopes7.append(slope)
    plt.plot(x0.tout[0], np.log2(x0.all["Obs_All"] / x0.all["Obs_All"][0]), 'red', lw=2)
    plt.plot(x0.tout[0][-1] + x1.tout[0],
             np.log2(x1.all["Obs_All"] / x0.all["Obs_All"][0]),
             'red', lw=2, label = "50% PC9-BR1")

###

model.parameters["PC9_VU_0"].value = num_cells * 0
model.parameters["PC9_BR1_0"].value = num_cells * 1

slopes8 = []
for sim in range(n_sims):
    t0 = np.linspace(0, pre_drug_time, pre_drug_time + 1)
    sim0 = BngSimulator(model, tspan=t0, verbose=False)
    x0 = sim0.run(tspan=t0, verbose=False, cleanup=False,
                  param_values={"k_dip_PC9_VU": 0.023 * np.log(2)}, method='ssa')

    PC9_VU_0.value = x0.all["Obs_PC9_VU"][-1]
    PC9_BR1_0.value = x0.all["Obs_PC9_BR1"][-1]

    t = np.linspace(0, post_drug_time, post_drug_time + 1)
    sim = BngSimulator(model, tspan=t, verbose=False)
    x1 = sim.run(tspan=t, verbose=False, cleanup=False,
                 param_values={"k_dip_PC9_VU": model.parameters['k_dip_PC9_VU'].value},
                 method='ssa')

    slope, intercept, r_value, p_value, std_err = \
        ss.linregress(x1.tout[0], np.log2(x1.all["Obs_All"] / x1.all["Obs_All"][0]))

    slopes8.append(slope)
    plt.plot(x0.tout[0], np.log2(x0.all["Obs_All"] / x0.all["Obs_All"][0]), 'red', lw=2)
    plt.plot(x0.tout[0][-1] + x1.tout[0],
             np.log2(x1.all["Obs_All"] / x0.all["Obs_All"][0]),
             'red', lw=2, label = "PC9-BR1")


# cell_tot = [x[-1]["Obs_All"] for x in x0.all]
# cell_PC9_VU = [x[-1]["Obs_PC9_VU"] for x in x0.all]
# cell_PC9_BR1 = [x[-1]["Obs_PC9_BR1"] for x in x0.all]
#
# print(cell_tot, cell_PC9_VU, cell_PC9_BR1)
# quit()
#
# for ind, cell in enumerate(cell_tot):
#     t = np.linspace(0, 125, 126)
#     sim = BngSimulator(model, tspan=t, verbose=False)
#     x1 = sim.run(tspan=t, verbose=False, cleanup=False,
#                  initials={model.species[0]:cell},
#                  param_values={"k_dip_PC9_VU": model.parameters['k_dip_PC9_VU'].value},
#                  method='ssa')
#     plt.plot(x0.tout[ind], np.log2(x0.all[ind]["Obs_All"] / x0.all[ind]["Obs_All"][0]), 'red', lw=2)
#     plt.plot(x0.tout[ind][-1] + x1.tout[0], np.log2(x1.all["Obs_All"]/x0.all[ind]["Obs_All"][0]), 'red', lw=2)

# plot_mean_min_max('Obs_All', color = 'red')
# y = odesolve(model = model,tspan = t, verbose = True)
# plt.plot(t, np.log2(y["Obs_All"]/y["Obs_All"][0]), 'r-', lw=4, label = "PC9-VU")


# model.parameters["PC9_VU_0"].value = num_cells * 0.999
# model.parameters["PC9_BR1_0"].value = num_cells * 0.001


#
# t = np.linspace(0, 150, 151)
# sims = 10
# sim = BngSimulator(model, tspan=t, verbose=False)
# x1 = sim.run(tspan=t, verbose=False, n_runs=sims, cleanup=False, method='ssa')
# # tout = x1.tout
# y1 = np.array(x1.observables)
#
# plot_mean_min_max('Obs_All', color = 'green')
# y = odesolve(model = model,tspan = t, verbose = True)
# plt.plot(t, np.log2(y["Obs_All"]/y["Obs_All"][0]), 'g-', lw=4, label = "0.1% PC9-BR1")
#
#
# model.parameters["PC9_VU_0"].value = num_cells * 0.99
# model.parameters["PC9_BR1_0"].value = num_cells * 0.01
#
# t = np.linspace(0, 150, 151)
# sims = 10
# sim = BngSimulator(model, tspan=t, verbose=False)
# x1 = sim.run(tspan=t, verbose=False, n_runs=sims, cleanup=False, method='ssa')
# # tout = x1.tout
# y1 = np.array(x1.observables)
#
# plot_mean_min_max('Obs_All', color = 'blue')
# y = odesolve(model = model,tspan = t, verbose = True)
# plt.plot(t, np.log2(y["Obs_All"]/y["Obs_All"][0]), 'b-', lw=4, label = "1% PC9-BR1")
#
#
# model.parameters["PC9_VU_0"].value = num_cells * 0.95
# model.parameters["PC9_BR1_0"].value = num_cells * 0.05
#
# t = np.linspace(0, 150, 151)
# sims = 10
# sim = BngSimulator(model, tspan=t, verbose=False)
# x1 = sim.run(tspan=t, verbose=False, n_runs=sims, cleanup=False, method='ssa')
# # tout = x1.tout
# y1 = np.array(x1.observables)
#
# plot_mean_min_max('Obs_All', color = 'magenta')
# y = odesolve(model = model,tspan = t, verbose = True)
# plt.plot(t, np.log2(y["Obs_All"]/y["Obs_All"][0]), 'm-', lw=4, label = "5% PC9-BR1")
#
# model.parameters["PC9_VU_0"].value = num_cells * 0.90
# model.parameters["PC9_BR1_0"].value = num_cells * 0.10
#
# t = np.linspace(0, 150, 151)
# sims = 10
# sim = BngSimulator(model, tspan=t, verbose=False)
# x1 = sim.run(tspan=t, verbose=False, n_runs=sims, cleanup=False, method='ssa')
# # tout = x1.tout
# y1 = np.array(x1.observables)
#
# plot_mean_min_max('Obs_All', color = 'yellow')
# y = odesolve(model = model,tspan = t, verbose = True)
# plt.plot(t, np.log2(y["Obs_All"]/y["Obs_All"][0]), 'y-', lw=4, label = "10% PC9-BR1")
#
#
# model.parameters["PC9_VU_0"].value = num_cells * 0.75
# model.parameters["PC9_BR1_0"].value = num_cells * 0.25
#
# t = np.linspace(0, 150, 151)
# sims = 10
# sim = BngSimulator(model, tspan=t, verbose=False)
# x1 = sim.run(tspan=t, verbose=False, n_runs=sims, cleanup=False, method='ssa')
# # tout = x1.tout
# y1 = np.array(x1.observables)
#
# plot_mean_min_max('Obs_All', color = 'black')
# y = odesolve(model = model,tspan = t, verbose = True)
# plt.plot(t, np.log2(y["Obs_All"]/y["Obs_All"][0]), 'k-', lw=4, label = "25% PC9-BR1")
#
#
# model.parameters["PC9_VU_0"].value = num_cells * 0.50
# model.parameters["PC9_BR1_0"].value = num_cells * 0.50
#
# t = np.linspace(0, 150, 151)
# sims = 10
# sim = BngSimulator(model, tspan=t, verbose=False)
# x1 = sim.run(tspan=t, verbose=False, n_runs=sims, cleanup=False, method='ssa')
# # tout = x1.tout
# y1 = np.array(x1.observables)
#
# plot_mean_min_max('Obs_All', color = 'cyan')
# y = odesolve(model = model,tspan = t, verbose = True)
# plt.plot(t, np.log2(y["Obs_All"]/y["Obs_All"][0]), 'c-', lw=4, label = "50% PC9-BR1")
#
# model.parameters["PC9_VU_0"].value = num_cells * 0
# model.parameters["PC9_BR1_0"].value = num_cells * 1
#
# t = np.linspace(0, 150, 151)
# sims = 10
# sim = BngSimulator(model, tspan=t, verbose=False)
# x1 = sim.run(tspan=t, verbose=False, n_runs=sims, cleanup=False, method='ssa')
# # tout = x1.tout
# y1 = np.array(x1.observables)
#
# plot_mean_min_max('Obs_All', color = 'orange')
# y = odesolve(model = model,tspan = t, verbose = True)
# plt.plot(t, np.log2(y["Obs_All"]/y["Obs_All"][0]), color = 'orange', ls = '-', lw=4, label = "PC9-BR1")

# plt.xlim(0,pre_drug_time+post_drug_time+(0.05*(pre_drug_time+post_drug_time)))
# plt.ylim(-0.5,5)
plt.xlabel("Time (hours)")
plt.ylabel("Population Doublings")
plt.title("PC9 Spike-In Model", weight = "bold")

handles, labels = plt.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys(), loc=0, title="Cell Line")
# plt.legend(loc=0, title="Cell Line")
plt.savefig("PC9_SpikeInModel_25hrdelay_750hrtotal.pdf")

print(slopes1)
print(slopes2)
print(slopes3)
print(slopes4)
print(slopes5)
print(slopes6)
print(slopes7)
print(slopes8)

plt.show()