import pandas as pd
from pylab import *
from numpy import *
from numpy.random import rand
import numpy as np
from scipy.optimize import minimize, fmin
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy import integrate
from pysb import *
from pysb.integrate import Solver
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

mixing_p9 = pd.read_csv('/Users/Corey/Documents/QuarantaLab/Mixing/mixing_Psg9_timenl2', sep = "\t")
mixing_p13 = pd.read_csv('/Users/Corey/Documents/QuarantaLab/Mixing/mixing_Psg9_timenl2', sep = "\t")
mixing_p22 = pd.read_csv('/Users/Corey/Documents/QuarantaLab/Mixing/mixing_Psg9_timenl2', sep = "\t")
mixing_p26 = pd.read_csv('/Users/Corey/Documents/QuarantaLab/Mixing/mixing_Psg9_timenl2', sep = "\t")
#print mixing_p9
#print mixing_p13
print mixing_p22

quit()

mixing_p9_list = []


mixing_p9_list.append(mixing_p9.loc[mixing_p9['CellLine'] == "SKMEL5 Subclone01"])
print mixing_p9_list



# Open an empty list of averaged (from 3 replicates) normalized log2 cell counts
nl2_average_p9 = []
time_points_p9 = 19 # Replace with however many time points in experimental data
#time_points = 17 # for mixing_all_rep

# Loop over the mixing list and create sublists corresponding to the experimental replicates, and
# take the average, which will be appended to the empty list above. Loop for all conditions, and 
# plot the experimental replicates and average.
for a,b in enumerate(mixing_p9_list):
    sublist = [b[x:x+time_points_p9] for x in range(0, len(b), time_points_p9)]
    rep1_p9 = sublist[0] # first replicate data points
    rep2_p9 = sublist[1]
    rep3_p9 = sublist[2]
    print rep1_p9
    print type(rep1_p9)
    print rep1_p9.keys()
#     quit()
    print rep2_p9
    print rep3_p9
#     quit()
    # Create a time list for experimental and computational results -- same for all
    time_p9 = [x for x in rep1_p9["Time"]] 
    # Make lists numpy arrays so can easily do math on them
    nl2_rep1_p9 = np.array([x for x in rep1_p9["nl2"]])
    nl2_rep2_p9 = np.array([x for x in rep2_p9["nl2"]])
    nl2_rep3_p9 = np.array([x for x in rep3_p9["nl2"]])

    nl2_average_p9.append((nl2_rep1_p9 + nl2_rep2_p9 + nl2_rep3_p9) / 3.0)
    
    plt.figure(a)
    
    plt.plot(time_p9, nl2_average_p9[-1], ms = 0.50, mfc = "0.25")
    
#     plt.plot(time_p9, nl2_rep1_p9, '*', ms = 12, mfc = "0.75")
#     plt.plot(time_p9, nl2_rep2_p9, '*', ms = 12, mfc = "0.75")
#     plt.plot(time_p9, nl2_rep3_p9, '*', ms = 12, mfc = "0.75")  




mixing_p13_list = []


mixing_p13_list.append(mixing_p13.loc[mixing_p13['CellLine'] == "SKMEL5 Subclone01"])
print mixing_p13_list



# Open an empty list of averaged (from 3 replicates) normalized log2 cell counts
nl2_average_p13 = []
time_points_p13 = 19 # Replace with however many time points in experimental data
#time_points = 17 # for mixing_all_rep

# Loop over the mixing list and create sublists corresponding to the experimental replicates, and
# take the average, which will be appended to the empty list above. Loop for all conditions, and 
# plot the experimental replicates and average.
for c,d in enumerate(mixing_p13_list):
    sublist = [d[x:x+time_points_p13] for x in range(0, len(d), time_points_p13)]
    rep1_p13 = sublist[0] # first replicate data points
    rep2_p13 = sublist[1]
    rep3_p13 = sublist[2]
    print rep1_p13
    print type(rep1_p13)
    print rep1_p13.keys()
#     quit()
    print rep2_p13
    print rep3_p13
#     quit()
    # Create a time list for experimental and computational results -- same for all
    time_p13 = [x for x in rep1_p13["Time"]] 
    # Make lists numpy arrays so can easily do math on them
    nl2_rep1_p13 = np.array([x for x in rep1_p13["nl2"]])
    nl2_rep2_p13 = np.array([x for x in rep2_p13["nl2"]])
    nl2_rep3_p13 = np.array([x for x in rep3_p13["nl2"]])

    nl2_average_p13.append((nl2_rep1_p13 + nl2_rep2_p13 + nl2_rep3_p13) / 3.0)
    
    plt.figure(c)
    
    plt.plot(time_p13, nl2_average_p13[-1], ms = 0.50, mfc = "0.25")
    
#     plt.plot(time_p13, nl2_rep1_p13, '*', ms = 12, mfc = "0.75")
#     plt.plot(time_p13, nl2_rep2_p13, '*', ms = 12, mfc = "0.75")
#     plt.plot(time_p13, nl2_rep3_p13, '*', ms = 12, mfc = "0.75")  


mixing_p22_list = []


mixing_p22_list.append(mixing_p22.loc[mixing_p22['CellLine'] == "SKMEL5 Subclone01"])
print mixing_p22_list



# Open an empty list of averaged (from 3 replicates) normalized log2 cell counts
nl2_average_p22 = []
time_points_p22 = 19 # Replace with however many time points in experimental data
#time_points = 17 # for mixing_all_rep

# Loop over the mixing list and create sublists corresponding to the experimental replicates, and
# take the average, which will be appended to the empty list above. Loop for all conditions, and 
# plot the experimental replicates and average.
for e,f in enumerate(mixing_p22_list):
    sublist = [f[x:x+time_points_p22] for x in range(0, len(f), time_points_p22)]
    rep1_p22 = sublist[0] # first replicate data points
    rep2_p22 = sublist[1]
    rep3_p22 = sublist[2]
    print rep1_p22
    print type(rep1_p22)
    print rep1_p22.keys()
#     quit()
    print rep2_p22
    print rep3_p22
#     quit()
    # Create a time list for experimental and computational results -- same for all
    time_p22 = [x for x in rep1_p22["Time"]] 
    # Make lists numpy arrays so can easily do math on them
    nl2_rep1_p22 = np.array([x for x in rep1_p22["nl2"]])
    nl2_rep2_p22 = np.array([x for x in rep2_p22["nl2"]])
    nl2_rep3_p22 = np.array([x for x in rep3_p22["nl2"]])

    nl2_average_p22.append((nl2_rep1_p22 + nl2_rep2_p22 + nl2_rep3_p22) / 3.0)
    
    plt.figure(e)
    
    plt.plot(time_p22, nl2_average_p22[-1], ms = 0.50, mfc = "0.25")
    
#     plt.plot(time_p22, nl2_rep1_p22, '*', ms = 12, mfc = "0.75")
#     plt.plot(time_p22, nl2_rep2_p22, '*', ms = 12, mfc = "0.75")
#     plt.plot(time_p22, nl2_rep3_p22, '*', ms = 12, mfc = "0.75")  


mixing_p26_list = []


mixing_p26_list.append(mixing_p26.loc[mixing_p26['CellLine'] == "SKMEL5 Subclone01"])
print mixing_p26_list



# Open an empty list of averaged (from 3 replicates) normalized log2 cell counts
nl2_average_p26 = []
time_points_p26 = 19 # Replace with however many time points in experimental data
#time_points = 17 # for mixing_all_rep

# Loop over the mixing list and create sublists corresponding to the experimental replicates, and
# take the average, which will be appended to the empty list above. Loop for all conditions, and 
# plot the experimental replicates and average.
for g,h in enumerate(mixing_p26_list):
    sublist = [h[x:x+time_points_p26] for x in range(0, len(h), time_points_p26)]
    rep1_p26 = sublist[0] # first replicate data points
    rep2_p26 = sublist[1]
    rep3_p26 = sublist[2]
    print rep1_p26
    print type(rep1_p26)
    print rep1_p26.keys()
#     quit()
    print rep2_p26
    print rep3_p26
#     quit()
    # Create a time list for experimental and computational results -- same for all
    time_p26 = [x for x in rep1_p26["Time"]] 
    # Make lists numpy arrays so can easily do math on them
    nl2_rep1_p26 = np.array([x for x in rep1_p26["nl2"]])
    nl2_rep2_p26 = np.array([x for x in rep2_p26["nl2"]])
    nl2_rep3_p26 = np.array([x for x in rep3_p26["nl2"]])

    nl2_average_p26.append((nl2_rep1_p26 + nl2_rep2_p26 + nl2_rep3_p26) / 3.0)
    
    plt.figure(g)
    
    plt.plot(time_p26, nl2_average_p26[-1], ms = 0.50, mfc = "0.25")
    
#     plt.plot(time_p26, nl2_rep1_p26, '*', ms = 12, mfc = "0.75")
#     plt.plot(time_p26, nl2_rep2_p26, '*', ms = 12, mfc = "0.75")
#     plt.plot(time_p26, nl2_rep3_p26, '*', ms = 12, mfc = "0.75")  


plt.ylim([-1.0,1.5])
plt.show()