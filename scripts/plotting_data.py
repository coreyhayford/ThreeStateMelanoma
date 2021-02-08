import numpy as np
import scipy.stats as ss
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.cm
import pandas as pd

# colors = matplotlib.colors.cnames.keys()
# colors = ['blue', 'green', 'red', 'cyan', 'magenta', 'DarkOrange', 'indigo', 'maroon', 'black']
# colors = ['0.0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8']
#colors = [matplotlib.cm.spectral(i) for i in np.linspace(0.9, 0.2, 9)]
#index = 0

data = pd.read_csv('/Users/Corey/Documents/QuarantaLab/mixing_mar24_processed.csv', 
                         delimiter='\t', dtype=None, names=True)

print data