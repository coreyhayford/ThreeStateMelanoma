# Uses LSD_model
# Prints distribution of DIP rates multi-well plate.

from pysb.core import *
from pysb.simulator.bng import BngSimulator
from pylab import *
import numpy as np
from numpy import linspace
import scipy.stats as ss
import matplotlib.pyplot as plt
import scipy as sp
import pandas as pd
import math
import seaborn as sns

sns.set(font_scale = 1.25)
sns.set_style("whitegrid")

# from pysb import *  # import everything from pysb
# from pysb.util import *  # for alias_model_components()
#
# def declare_monomers():  # Monomer("name_of_monomer")
#     Monomer("A")
#     Monomer("B")
#     Monomer("C")
#
#
# def declare_parameters():  # Parameter("name_of_par", value)
#     Parameter("A0")
#     Parameter("B0")
#     Parameter("C0")
#
#     Parameter("k_div_A")#, .033)  # np.random.normal(.03, .005) )
#     Parameter("k_div_B")#, .033)  # np.random.normal(.05,.005) )
#     Parameter("k_div_C")#, .033)  # np.random.normal(.07,.005) )
#
#     Parameter('k_death_A')#, .005)
#     Parameter('k_death_B')#, .005)
#     Parameter('k_death_C')#, .005)
#
#     Parameter("k_AB",     0.025)
#     Parameter("k_BA",     0.00004)
#     Parameter("k_CB",     0.025)
#     Parameter("k_BC",     0.00004)
#
#     alias_model_components()  # "Make all model components visible as symbols in the caller's global namespace"
#
# def declare_initial_conditions():  # Initial(call monomer's name, call parameter's name for IC value = number of initial cells)
#     Initial(A, A0)
#     Initial(B, B0)
#     Initial(C, C0)
#
#
# def declare_observables():  # Observables("name_of_obs", call monomer name you are observing)
#     Observable("Obs_A", A)
#     Observable("Obs_B", B)
#     Observable("Obs_C", C)
#     Observable("Obs_AB", A() + B())
#     Observable("Obs_BC", B() + C())
#     Observable("Obs_AC", A() + C())
#     Observable("Obs_All", A() + B() + C())  ## Don't forget parentheses here
#
#
# def declare_rules():  # Rule("rule_name", chemical reaction, call name of DIP rate/parameter)
#     Rule("DIV_A", A() >> A() + A(), k_div_A)  ##Don't forget parentheses here
#     Rule("DIV_B", B() >> B() + B(), k_div_B)
#     Rule("DIV_C", C() >> C() + C(), k_div_C)
#     Rule('Death_A', A() >> None, k_death_A)
#     Rule('Death_B', B() >> None, k_death_B)
#     Rule("Death_C", C() >> None, k_death_C)
#
#
# def declare_transition_rules():
#     Rule("A_to_B", A() >> B(), k_AB)
#     Rule("B_to_A", B() >> A(), k_BA)
#     Rule("B_to_C", B() >> C(), k_BC)
#     Rule("C_to_B", C() >> B(), k_CB)

import LSD_model as pre

def Transition_model(prob, n_sim, n_cells):

    distr_all = []
    count = 0
    for i in range(n_sim):
        # plt.figure()
        Model()
        pre.declare_monomers()
        pre.declare_parameters()
        pre.declare_initial_conditions()
        pre.declare_observables()
        pre.declare_rules()

        # plt.figure()
        ## LCS
        num_cells = np.random.poisson(n_cells)
        ## cFP
        # num_cells = 1

        picks = np.random.multinomial(num_cells, prob)
        picks_total = np.sum(picks)
        if picks_total != 0:
            A_init = picks[0]
            B_init = picks[1]
            C_init = picks[2]

            A0.value = A_init
            B0.value = B_init
            C0.value = C_init

            k_div_A.value = .033
            k_div_B.value = .033
            k_div_C.value = .033

            k_death_A.value = .005
            k_death_B.value = .005
            k_death_C.value = .005
            print "Pre-Growth"
            print A0, B0, C0
            # Set up space
            t = linspace(0, 200, 201)  # set up the time axis from 0 to 100 in 101 increments
            sim = BngSimulator(model, tspan=t, verbose=False)
            y = sim.run(nsims=1, verbose=False)

            # plt.plot(y.tout[0], y.all["Obs_All"])
            plt.plot(y.tout[0], np.log2(y.all["Obs_All"] / y.all["Obs_All"][0]), '0.5', lw=2)
    # plt.show()

            #########
            A0.value = y.all["Obs_A"][-1]
            B0.value = y.all["Obs_B"][-1]
            C0.value = y.all["Obs_C"][-1]
            print "Post-Growth"
            print A0, B0, C0

            k_div_A.value = .033
            k_div_B.value = .033
            k_div_C.value = .033
            #

            k_death_A.value = .043
            k_death_B.value = .033
            k_death_C.value = .023

            pre.declare_transition_rules()
            k_AB.value = 0.025
            k_BA.value = 0.00004
            k_CB.value = 0.025
            k_BC.value = 0.00004

            #########
            t1 = linspace(0, 200, 201)
            sim1 = BngSimulator(model, tspan=t1, verbose=False)
            y1 = sim1.run(nsims=1, verbose=False)

            print "Post-Drug"
            print(y1.all["Obs_All"][-1])
            print(y1.all["Obs_A"][-1], y1.all["Obs_B"][-1], y1.all["Obs_C"][-1])

            count = count + 1
            print "Run #"
            print count

            # plt.plot(y.tout[0][-1] + y1.tout[0],
            #          np.log2(y1.all["Obs_All"] / y.all["Obs_All"][0]),
            #          '0.5', lw=2)
            # plt.xlabel("Time (hours)")
            # plt.ylabel("Population Doublings")
            # plt.title("Model Trajectories")

    #         plt.plot(y.tout[0][-1] + y1.tout[0][:25],
    #                  np.log2(y1.all["Obs_All"][:25] / y1.all["Obs_All"][0]),
    #                  lw=2)
    # plt.show()
    # quit()

            if picks_total != 0:
                slope_all, intercept_all, r_value_all, p_value_all, std_err_all = \
                    ss.linregress(y1.tout[0][75:125], np.log2(y1.all["Obs_All"][75:125] / y1.all["Obs_All"][0]))
                if math.isnan(slope_all) == False:
                    distr_all = distr_all + [slope_all]
                print "All", slope_all

    data = np.array(distr_all)

    plt.figure("DIP Rate Distribution")
    sns.distplot([distr_all], kde=False, color='k', hist_kws={"alpha": 0.5})
    # plt.hist([distr_all], color='k', alpha=0.5)
    plt.xlabel("DIP Rate", fontsize = 16, weight="bold")
    plt.ylabel("Frequency", fontsize = 16, weight="bold")
    plt.tick_params(labelsize = 10)
    plt.xlim(-0.03, 0.03)
    # plt.ylim(0, 125)
    plt.title("Three-State LSD Model with Probabilities %s, %s, %s" % (str(prob[0]), str(prob[1]), str(prob[2])))
    np.save('%s_%s_%s_1kData_LSD_%dcells.npy' % (str(prob[0]), str(prob[1]), str(prob[2]), n_cells), data)
    plt.savefig('3S_LSDModel1k_%s_%s_%s_%dcells.pdf' % (str(prob[0]), str(prob[1]), str(prob[2]), n_cells))
    plt.close()

print "Model #1"
Transition_model(prob = [1.0, .0, .0], n_sim=1000, n_cells=5)
Transition_model(prob = [1.0, .0, .0], n_sim=1000, n_cells=10)
Transition_model(prob = [1.0, .0, .0], n_sim=1000, n_cells=25)
print "Model #2"
Transition_model(prob = [0.85, .10, .05], n_sim=1000, n_cells=5)
Transition_model(prob = [0.85, .10, .05], n_sim=1000, n_cells=10)
Transition_model(prob = [0.85, .10, .05], n_sim=1000, n_cells=25)
print "Model #3"
Transition_model(prob = [0.70, .20, .10], n_sim=1000, n_cells=5)
Transition_model(prob = [0.70, .20, .10], n_sim=1000, n_cells=10)
Transition_model(prob = [0.70, .20, .10], n_sim=1000, n_cells=25)
print "Model #4"
Transition_model(prob = [0.25, .50, .25], n_sim=1000, n_cells=5)
Transition_model(prob = [0.25, .50, .25], n_sim=1000, n_cells=10)
Transition_model(prob = [0.25, .50, .25], n_sim=1000, n_cells=25)