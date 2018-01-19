from pysb import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sp
from pysb.bng import run_ssa
from pysb.simulator.bng import BngSimulator
from pysb.bng import generate_equations
import pandas as pd
import seaborn as sns
import math

sns.set(font_scale = 1.5)
sns.set_style("white")

LSD_spikein_75resistant_resistant = np.load("LSDspikein-rep2_resistantDIPs_0.750resistant.npy")
LSD_spikein_75resistant_normal = np.load("LSDspikein-rep2_normalDIPs_0.750resistant.npy")

cFP_spikein_75resistant_resistant = np.load("cFPspikein-rep2_resistantDIPs_0.750resistant.npy")
cFP_spikein_75resistant_normal = np.load("cFPspikein-rep2_normalDIPs_0.750resistant.npy")

bins = np.histogram(np.hstack((LSD_spikein_75resistant_resistant, LSD_spikein_75resistant_normal,
                               cFP_spikein_75resistant_resistant, cFP_spikein_75resistant_normal)),
                    bins = 100)[1]

plt.figure()
plt.hist(LSD_spikein_75resistant_resistant, bins, color= 'r', alpha = 0.5, label = 'LSD Resistant')
plt.hist(LSD_spikein_75resistant_normal, bins, color = 'b', alpha = 0.5, label = 'LSD Normal')
plt.hist(cFP_spikein_75resistant_resistant, bins, color = 'g', alpha = 0.5, label = 'cFP Resistant')
plt.hist(cFP_spikein_75resistant_normal, bins, color = 'm', alpha = 0.5, label = 'cFP Normal')
plt.xlabel('DIP Rate', weight = "bold")
plt.ylabel('Frequency', weight = "bold")
plt.xlim(-0.04, 0.06)
plt.ylim(0,240)
plt.legend(loc = 'upper center')
plt.title('25% PC9-VU | 75% PC9-BR1', weight = "bold")
plt.savefig('SpikeInComparison_25percentPC9-VU_75percentPC9-BR1.pdf')

LSD_spikein_50resistant_resistant = np.load("LSDspikein-rep2_resistantDIPs_0.500resistant.npy")
LSD_spikein_50resistant_normal = np.load("LSDspikein-rep2_normalDIPs_0.500resistant.npy")

cFP_spikein_50resistant_resistant = np.load("cFPspikein-rep2_resistantDIPs_0.500resistant.npy")
cFP_spikein_50resistant_normal = np.load("cFPspikein-rep2_normalDIPs_0.500resistant.npy")

bins = np.histogram(np.hstack((LSD_spikein_50resistant_resistant, LSD_spikein_50resistant_normal,
                               cFP_spikein_50resistant_resistant, cFP_spikein_50resistant_normal)),
                    bins = 100)[1]

plt.figure()
plt.hist(LSD_spikein_50resistant_resistant, bins, color= 'r', alpha = 0.5, label = 'LSD Resistant')
plt.hist(LSD_spikein_50resistant_normal, bins, color = 'b', alpha = 0.5, label = 'LSD Normal')
plt.hist(cFP_spikein_50resistant_resistant, bins, color = 'g', alpha = 0.5, label = 'cFP Resistant')
plt.hist(cFP_spikein_50resistant_normal, bins, color = 'm', alpha = 0.5, label = 'cFP Normal')
plt.xlabel('DIP Rate', weight = "bold")
plt.ylabel('Frequency', weight = "bold")
plt.xlim(-0.04, 0.06)
plt.ylim(0,240)
plt.legend(loc = 'upper center')
plt.title('50% PC9-VU | 50% PC9-BR1', weight = "bold")
plt.savefig('SpikeInComparison_50percentPC9-VU_50percentPC9-BR1.pdf')

LSD_spikein_25resistant_resistant = np.load("LSDspikein-rep2_resistantDIPs_0.250resistant.npy")
LSD_spikein_25resistant_normal = np.load("LSDspikein-rep2_normalDIPs_0.250resistant.npy")

cFP_spikein_25resistant_resistant = np.load("cFPspikein-rep2_resistantDIPs_0.250resistant.npy")
cFP_spikein_25resistant_normal = np.load("cFPspikein-rep2_normalDIPs_0.250resistant.npy")

bins = np.histogram(np.hstack((LSD_spikein_25resistant_resistant, LSD_spikein_25resistant_normal,
                               cFP_spikein_25resistant_resistant, cFP_spikein_25resistant_normal)),
                    bins = 100)[1]

plt.figure()
plt.hist(LSD_spikein_25resistant_resistant, bins, color= 'r', alpha = 0.5, label = 'LSD Resistant')
plt.hist(LSD_spikein_25resistant_normal, bins, color = 'b', alpha = 0.5, label = 'LSD Normal')
plt.hist(cFP_spikein_25resistant_resistant, bins, color = 'g', alpha = 0.5, label = 'cFP Resistant')
plt.hist(cFP_spikein_25resistant_normal, bins, color = 'm', alpha = 0.5, label = 'cFP Normal')
plt.xlabel('DIP Rate', weight = "bold")
plt.ylabel('Frequency', weight = "bold")
plt.xlim(-0.04, 0.06)
plt.ylim(0,240)
plt.legend(loc = 'upper center')
plt.title('75% PC9-VU | 25% PC9-BR1', weight = "bold")
plt.savefig('SpikeInComparison_75percentPC9-VU_25percentPC9-BR1.pdf')

LSD_spikein_20resistant_resistant = np.load("LSDspikein-rep2_resistantDIPs_0.200resistant.npy")
LSD_spikein_20resistant_normal = np.load("LSDspikein-rep2_normalDIPs_0.200resistant.npy")

cFP_spikein_20resistant_resistant = np.load("cFPspikein-rep2_resistantDIPs_0.200resistant.npy")
cFP_spikein_20resistant_normal = np.load("cFPspikein-rep2_normalDIPs_0.200resistant.npy")

bins = np.histogram(np.hstack((LSD_spikein_20resistant_resistant, LSD_spikein_20resistant_normal,
                               cFP_spikein_20resistant_resistant, cFP_spikein_20resistant_normal)),
                    bins = 100)[1]

plt.figure()
plt.hist(LSD_spikein_20resistant_resistant, bins, color= 'r', alpha = 0.5, label = 'LSD Resistant')
plt.hist(LSD_spikein_20resistant_normal, bins, color = 'b', alpha = 0.5, label = 'LSD Normal')
plt.hist(cFP_spikein_20resistant_resistant, bins, color = 'g', alpha = 0.5, label = 'cFP Resistant')
plt.hist(cFP_spikein_20resistant_normal, bins, color = 'm', alpha = 0.5, label = 'cFP Normal')
plt.xlabel('DIP Rate', weight = "bold")
plt.ylabel('Frequency', weight = "bold")
plt.xlim(-0.04, 0.06)
plt.ylim(0,240)
plt.legend(loc = 'upper center')
plt.title('80% PC9-VU | 20% PC9-BR1', weight = "bold")
plt.savefig('SpikeInComparison_80percentPC9-VU_20percentPC9-BR1.pdf')

LSD_spikein_10resistant_resistant = np.load("LSDspikein-rep2_resistantDIPs_0.100resistant.npy")
LSD_spikein_10resistant_normal = np.load("LSDspikein-rep2_normalDIPs_0.100resistant.npy")

cFP_spikein_10resistant_resistant = np.load("cFPspikein-rep2_resistantDIPs_0.100resistant.npy")
cFP_spikein_10resistant_normal = np.load("cFPspikein-rep2_normalDIPs_0.100resistant.npy")

bins = np.histogram(np.hstack((LSD_spikein_10resistant_resistant, LSD_spikein_10resistant_normal,
                               cFP_spikein_10resistant_resistant, cFP_spikein_10resistant_normal)),
                    bins = 100)[1]

plt.figure()
plt.hist(LSD_spikein_10resistant_resistant, bins, color= 'r', alpha = 0.5, label = 'LSD Resistant')
plt.hist(LSD_spikein_10resistant_normal, bins, color = 'b', alpha = 0.5, label = 'LSD Normal')
plt.hist(cFP_spikein_10resistant_resistant, bins, color = 'g', alpha = 0.5, label = 'cFP Resistant')
plt.hist(cFP_spikein_10resistant_normal, bins, color = 'm', alpha = 0.5, label = 'cFP Normal')
plt.xlabel('DIP Rate', weight = "bold")
plt.ylabel('Frequency', weight = "bold")
plt.xlim(-0.04, 0.06)
plt.ylim(0,240)
plt.legend(loc = 'upper center')
plt.title('90% PC9-VU | 10% PC9-BR1', weight = "bold")
plt.savefig('SpikeInComparison_90percentPC9-VU_10percentPC9-BR1.pdf')

LSD_spikein_05resistant_resistant = np.load("LSDspikein-rep2_resistantDIPs_0.050resistant.npy")
LSD_spikein_05resistant_normal = np.load("LSDspikein-rep2_normalDIPs_0.050resistant.npy")

cFP_spikein_05resistant_resistant = np.load("cFPspikein-rep2_resistantDIPs_0.050resistant.npy")
cFP_spikein_05resistant_normal = np.load("cFPspikein-rep2_normalDIPs_0.050resistant.npy")

bins = np.histogram(np.hstack((LSD_spikein_05resistant_resistant, LSD_spikein_05resistant_normal,
                               cFP_spikein_05resistant_resistant, cFP_spikein_05resistant_normal)),
                    bins = 100)[1]

plt.figure()
plt.hist(LSD_spikein_05resistant_resistant, bins, color= 'r', alpha = 0.5, label = 'LSD Resistant')
plt.hist(LSD_spikein_05resistant_normal, bins, color = 'b', alpha = 0.5, label = 'LSD Normal')
plt.hist(cFP_spikein_05resistant_resistant, bins, color = 'g', alpha = 0.5, label = 'cFP Resistant')
plt.hist(cFP_spikein_05resistant_normal, bins, color = 'm', alpha = 0.5, label = 'cFP Normal')
plt.xlabel('DIP Rate', weight = "bold")
plt.ylabel('Frequency', weight = "bold")
plt.xlim(-0.04, 0.06)
plt.ylim(0,240)
plt.legend(loc = 'upper center')
plt.title('95% PC9-VU | 5% PC9-BR1', weight = "bold")
plt.savefig('SpikeInComparison_95percentPC9-VU_5percentPC9-BR1.pdf')

LSD_spikein_01resistant_resistant = np.load("LSDspikein-rep2_resistantDIPs_0.010resistant.npy")
LSD_spikein_01resistant_normal = np.load("LSDspikein-rep2_normalDIPs_0.010resistant.npy")

cFP_spikein_01resistant_resistant = np.load("cFPspikein-rep2_resistantDIPs_0.010resistant.npy")
cFP_spikein_01resistant_normal = np.load("cFPspikein-rep2_normalDIPs_0.010resistant.npy")

bins = np.histogram(np.hstack((LSD_spikein_01resistant_resistant, LSD_spikein_01resistant_normal,
                               cFP_spikein_01resistant_resistant, cFP_spikein_01resistant_normal)),
                    bins = 100)[1]

plt.figure()
plt.hist(LSD_spikein_01resistant_resistant, bins, color= 'r', alpha = 0.5, label = 'LSD Resistant')
plt.hist(LSD_spikein_01resistant_normal, bins, color = 'b', alpha = 0.5, label = 'LSD Normal')
plt.hist(cFP_spikein_01resistant_resistant, bins, color = 'g', alpha = 0.5, label = 'cFP Resistant')
plt.hist(cFP_spikein_01resistant_normal, bins, color = 'm', alpha = 0.5, label = 'cFP Normal')
plt.xlabel('DIP Rate', weight = "bold")
plt.ylabel('Frequency', weight = "bold")
plt.xlim(-0.04, 0.06)
plt.ylim(0,240)
plt.legend(loc = 'upper center')
plt.title('99% PC9-VU | 01% PC9-BR1', weight = "bold")
plt.savefig('SpikeInComparison_99percentPC9-VU_01percentPC9-BR1.pdf')

LSD_spikein_001resistant_resistant = np.load("LSDspikein-rep2_resistantDIPs_0.001resistant.npy")
LSD_spikein_001resistant_normal = np.load("LSDspikein-rep2_normalDIPs_0.001resistant.npy")

cFP_spikein_001resistant_resistant = np.load("cFPspikein-rep2_resistantDIPs_0.001resistant.npy")
cFP_spikein_001resistant_normal = np.load("cFPspikein-rep2_normalDIPs_0.001resistant.npy")

bins = np.histogram(np.hstack((LSD_spikein_001resistant_resistant, LSD_spikein_001resistant_normal,
                               cFP_spikein_001resistant_resistant, cFP_spikein_001resistant_normal)),
                    bins = 100)[1]

plt.figure()
plt.hist(LSD_spikein_001resistant_resistant, bins, color= 'r', alpha = 0.5, label = 'LSD Resistant')
plt.hist(LSD_spikein_001resistant_normal, bins, color = 'b', alpha = 0.5, label = 'LSD Normal')
plt.hist(cFP_spikein_001resistant_resistant, bins, color = 'g', alpha = 0.5, label = 'cFP Resistant')
plt.hist(cFP_spikein_001resistant_normal, bins, color = 'm', alpha = 0.5, label = 'cFP Normal')
plt.xlabel('DIP Rate', weight = "bold")
plt.ylabel('Frequency', weight = "bold")
plt.xlim(-0.04, 0.06)
plt.ylim(0,240)
plt.legend(loc = 'upper center')
plt.title('99.9% PC9-VU | 0.1% PC9-BR1', weight = "bold")
plt.savefig('SpikeInComparison_99.9percentPC9-VU_0.1percentPC9-BR1.pdf')