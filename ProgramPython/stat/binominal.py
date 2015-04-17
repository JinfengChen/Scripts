from scipy.stats import binom_test

d = binom_test(51, 235, 1.0/6)
print d
