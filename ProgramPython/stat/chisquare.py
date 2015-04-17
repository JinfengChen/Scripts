import numpy as np
from scipy.stats.mstats import chisquare
from scipy.stats import fisher_exact

observed = np.array([2, 188])
expected = np.array([1, 80])
a = chisquare(observed, expected)
b = fisher_exact([observed, expected])

print a
print b
