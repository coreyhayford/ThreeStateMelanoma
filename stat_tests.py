import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from numpy import median

# Load model produced data
LCSdata1_trans = np.load("1.0_0.0_0.0_10kData.npy")
LCSdata2_trans = np.load("0.95_0.05_0.0_10kData.npy")
LCSdata3_trans = np.load("0.85_0.1_0.05_10kData.npy")
LCSdata4_trans = np.load("0.7_0.2_0.1_10kData.npy")
LCSdata5_trans = np.load("0.5_0.3_0.2_10kData.npy")
LCSdata6_trans = np.load("0.25_0.5_0.25_10kData.npy")

cFPdata1_trans = np.load("1.0_0.0_0.0_10kData_cFP.npy")
cFPdata2_trans = np.load("0.95_0.05_0.0_10kData_cFP.npy")
cFPdata3_trans = np.load("0.85_0.1_0.05_10kData_cFP.npy")
cFPdata4_trans = np.load("0.7_0.2_0.1_10kData_cFP.npy")
cFPdata5_trans = np.load("0.5_0.3_0.2_10kData_cFP.npy")
cFPdata6_trans = np.load("0.25_0.5_0.25_10kData_cFP.npy")

LCSdata1_notrans = np.load("1.0_0.0_0.0_10kData_LCS_noTrans.npy")
LCSdata2_notrans = np.load("0.95_0.05_0.0_10kData_LCS_noTrans.npy")
LCSdata3_notrans = np.load("0.85_0.1_0.05_10kData_LCS_noTrans.npy")
LCSdata4_notrans = np.load("0.7_0.2_0.1_10kData_LCS_noTrans.npy")
LCSdata5_notrans = np.load("0.5_0.3_0.2_10kData_LCS_noTrans.npy")
LCSdata6_notrans = np.load("0.25_0.5_0.25_10kData_LCS_noTrans.npy")

cFPdata1_notrans = np.load("1.0_0.0_0.0_10kData_cFP_noTrans.npy")
cFPdata2_notrans = np.load("0.95_0.05_0.0_10kData_cFP_noTrans.npy")
cFPdata3_notrans = np.load("0.85_0.1_0.05_10kData_cFP_noTrans.npy")
cFPdata4_notrans = np.load("0.7_0.2_0.1_10kData_cFP_noTrans.npy")
cFPdata5_notrans = np.load("0.5_0.3_0.2_10kData_cFP_noTrans.npy")
cFPdata6_notrans = np.load("0.25_0.5_0.25_10kData_cFP_noTrans.npy")

expData_p11 = np.loadtxt("experimental_psg11.2.txt", delimiter='\t')
expData_p15 = np.loadtxt("experimental_psg15.2.txt", delimiter='\t')
expData_p19 = np.loadtxt("experimental_psg19.2.txt", delimiter='\t')
expData_p28 = np.loadtxt("experimental_psg28.2.txt", delimiter='\t')
expData_par = np.loadtxt("experimental_parental2.txt", delimiter='\t')

