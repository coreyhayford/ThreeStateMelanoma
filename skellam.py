from scipy.stats import skellam
import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots()

mu1, mu2 = 0.03, 0.02
mean, var, skew, kurt = skellam.stats(mu1, mu2, moments='mvsk')

print(mean, var, skew, kurt)

x = np.arange(skellam.ppf(0.01, mu1, mu2),
              skellam.ppf(0.99, mu1, mu2))

ax.plot(x, skellam.pmf(x, mu1, mu2), 'bo', ms=8, label = 'skellam pmf')
ax.vlines(x, 0, skellam.pmf(x, mu1, mu2), colors='b', lw = 5, alpha=0.5)
plt.show()