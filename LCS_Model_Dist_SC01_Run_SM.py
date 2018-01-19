# Uses LCS_stochastic_model
# Prints color coded distribution of DIP rates for 384 well plate. DIP rate = highest proliferative cell in well DIP rate

from pysb import *
from pysb.util import *
from pysb.macros import *
from pysb.bng import *
from pysb.core import *
from pysb.integrate import odesolve
from pylab import *
import pylab as pl
import numpy as np
from numpy import linspace
from sympy import sympify
from scipy import constants
import scipy.stats as ss
import matplotlib.pyplot as plt
import scipy as sp
import pandas as pd
import math

# Import model from other file
# import Postdrug_Stochastic_Model as post
# import Predrug_stochastic_model as pre
import LCS_Stochastic_Model_SM as pre


# def plot(t, t_start=0, label=False):
#     plt.plot(t + t_start, np.log2(y["Obs_A"] / A_init), 'c-', lw=2,
#              label=(
#              'A DIP = %g, IC = %d (%d)' % (k_div_A.value - k_death_A.value, A0.value, picks[0])) if label else None)
#     plt.plot(t + t_start, np.log2(y["Obs_B"] / B_init), 'r-', lw=2,
#              label=(
#              'B DIP = %g, IC = %d (%d)' % (k_div_B.value - k_death_B.value, B0.value, picks[1])) if label else None)
#     plt.plot(t + t_start, np.log2(y["Obs_C"] / C_init), 'g-', lw=2,
#              label=(
#              'C DIP = %g, IC = %d (%d)' % (k_div_C.value - k_death_C.value, C0.value, picks[2])) if label else None)
#     plt.plot(t + t_start, np.log2(y["Obs_All"] / (A_init + B_init + C_init)), 'k:', lw=4,
#              label=('All' if label else None))
#     plt.xlabel("Time")
#     plt.ylabel("Population Doublings")
#     plt.legend(loc=0)  # , fontsize = 6)


# Use these picks as IC for predrug model and run but don't print
# Use y["Obs_A"][200] etc. as IC for postdrug model

# seed = 20
# np.random.seed(seed)


def Transition_model(prob):
    distr_a = []
    distr_b = []
    distr_c = []
    distr_all = []
    count = 0
    # Model()
    # pre.declare_monomers()
    # pre.declare_parameters()
    # pre.declare_initial_conditions()
    # pre.declare_observables()
    # pre.declare_rules()
    for i in range(10000):
        # plt.figure()
        Model()
        pre.declare_monomers()
        pre.declare_parameters()
        pre.declare_initial_conditions()
        pre.declare_observables()
        pre.declare_rules()


        # plt.figure()
        ## LCS
        num_cells = np.random.poisson(3)
        ## cFP
        # num_cells = 1

        picks = np.random.multinomial(num_cells, prob)
        print picks
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
            #

            k_death_A.value = .005
            k_death_B.value = .005
            k_death_C.value = .005
            print A0, B0, C0
            # Set up space
            t = linspace(0, 200, 201)  # set up the time axis from 0 to 200 in 201 increments
            y = run_ssa(model, t_end=t[-1], n_steps=len(t) - 1, verbose=False)
            # y = run_ssa(model,t_end = t[-1], n_steps = len(t)-1, verbose=False, seed=seed)
            # y = odesolve(model, t, verbose=False)
            picks_total = np.sum(picks)
            #         print "A:", y["Obs_A"][-1]/y["Obs_All"][-1], "(%g)" % (1.*picks[0]/picks_total)
            #         print "B:", y["Obs_B"][-1]/y["Obs_All"][-1], "(%g)" % (1.*picks[1]/picks_total)
            #         print "C:", y["Obs_C"][-1]/y["Obs_All"][-1], "(%g)" % (1.*picks[2]/picks_total)

            # plot(t, label=False)

            # plt.axvline(x=t[-1], color='red', ls='--', lw=2)

            # print(model.parameters)

            #########
            A0.value = y["Obs_A"][-1]
            B0.value = y["Obs_B"][-1]
            C0.value = y["Obs_C"][-1]
            # print A0, B0, C0

            k_div_A.value = .033
            k_div_B.value = .033
            k_div_C.value = .033
            #

            k_death_A.value = .043
            k_death_B.value = .033
            k_death_C.value = .023

            # k_AB.value = 0.025
            # k_BA.value = 0.00004
            # k_BC.value = 0.025
            # k_CB.value = 0.00004

            # print(model.parameters)

            #########
            t = linspace(0, 75, 76)
            y = run_ssa(model, t_end=t[-1], n_steps=len(t) - 1, verbose=False)
            # y = run_ssa(model,t_end = t[-1], n_steps = len(t)-1, verbose=False, seed=seed)
            # y = odesolve(model, t, verbose=False)

            # plot(t, t_start=t[-1], label=True)

            # plot(t, t_start=t+200, label=True)

            if picks[2] != 0:
                slope, intercept, r_value, p_value, std_err = ss.linregress(t, np.log2(y["Obs_C"] / y["Obs_C"][0]))
                if math.isnan(slope) == False:
                    distr_c = distr_c + [slope]
                print "C", slope
            elif picks[1] != 0:
                slope, intercept, r_value, p_value, std_err = ss.linregress(t, np.log2(y["Obs_B"] / y["Obs_B"][0]))
                if math.isnan(slope) == False:
                    distr_b = distr_b + [slope]
                print "B", slope
            else:
                slope, intercept, r_value, p_value, std_err = ss.linregress(t, np.log2(y["Obs_A"] / y["Obs_A"][0]))
                if math.isnan(slope) == False:
                    distr_a = distr_a + [slope]
                print "A", slope
            count = count + 1
            print count

            slope_all, intercept_all, r_value_all, p_value_all, std_err_all = ss.linregress(t, np.log2(y["Obs_All"] / y["Obs_All"][0]))
            if math.isnan(slope_all) == False:
                distr_all = distr_all + [slope_all]
            print "All", slope_all

            print(np.log2(y["Obs_A"][-1]))
            print(np.log2(y["Obs_All"][-1]))

            print(np.log2(y["Obs_A"][-1] / y["Obs_A"][0]))
            print(np.log2(y["Obs_All"][-1] / y["Obs_All"][0]))
    # plt.show()
    # quit()
    fig = plt.figure("DIP Rate Distribution")
    print distr_a, distr_b, distr_c
    print distr_all
    # n, bin_edges = np.histogram(distr_c, bins = 10, range = (-.05, .05)) #bins -.01, 0, and .01
    # # Normalize it, so that every bins value gives the probability of that bin
    # bin_probability = n/float(n.sum())
    # # Get the mid points of every bin
    # #bin_middles = [-.01, 0, .01]
    # bin_middles = (bin_edges[1:]+bin_edges[:-1])/2.
    # # Compute the bin-width
    # #bin_width = [.01,.01,.01]
    # bin_width = bin_edges[1]-bin_edges[0]
    # # Plot the histogram as a bar plot
    # plt.bar(bin_middles, bin_probability, width=bin_width, align='center')
    # plt.grid(True)
    # plt.xlim(-0.05,0.05)
    print np.sort(distr_a)
    print np.sort(distr_b)
    print np.sort(distr_c)
    data = np.array(distr_all)
    print(data)

    # plt.hist([distr_a, distr_b, distr_c], bins=50, stacked=True,  # color = ["c", "r", "g"],
    #          label=["A: %d (%g)" % (len(distr_a), float(len(distr_a)) / count),
    #                 "B: %d (%g)" % (len(distr_b), float(len(distr_b)) / count),
    #                 "C: %d (%g)" % (len(distr_c), float(len(distr_c)) / count)])
    plt.hist([distr_all], bins = 50, alpha = 0.5, color = 'k')
             # label="Mean = %g; SD = %g" % (np.mean(distr_all), np.std(distr_all)))
    # plt.legend(loc="upper left")
    # # plt.grid(True)
    plt.axis([-0.05, 0.05, 0, 3000])
    plt.text(-0.04, 2800, r'$\mu=%g$' % (np.mean(distr_all)))
    plt.text(-0.04, 2700, r'$\sigma=%g$' % (np.std(distr_all)))
    plt.xlabel("DIP Rate", fontsize = 16)
    plt.ylabel("Frequency", fontsize = 16)
    plt.tick_params(labelsize = 10)
    # plt.xlim(-0.05, 0.05)
    # plt.ylim(0, 125)
    plt.title("LCS SC01 Model with Probabilities %s, %s, %s" % (str(prob[0]), str(prob[1]), str(prob[2])), fontsize = 18)
    np.save('%s_%s_%s_10kData_LCS_forselfdistance.npy' % (str(prob[0]), str(prob[1]), str(prob[2])), data)
    fig.savefig('LCSmodelSC01diversify10k_selfdist_%s_%s_%s.pdf' % (str(prob[0]), str(prob[1]), str(prob[2])))
    # plt.show()
    plt.close()

Transition_model(prob = [1.0, .0, .0])
# noTransition_model(prob = [1.0, .0, .0])
# noTransition_model(prob = [0.95, .05, .0])
# noTransition_model(prob = [0.85, .10, .05])
# noTransition_model(prob = [0.70, .20, .10])
# noTransition_model(prob = [0.50, .30, .20])
# noTransition_model(prob = [0.25, .50, .25])