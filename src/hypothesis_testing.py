
import numpy as np
from scipy.stats import kstest, norm
import matplotlib.pyplot as plt
from scipy import stats
from functools import partial

file_path = '../build/objLimitingDistribution.txt' 
with open(file_path, 'r') as file:
    numbers_str = file.read().strip()

# Convert the string of numbers to a list of floats
objectives = np.array( [float(num) for num in numbers_str.split(', ')] )
mu = objectives.mean()
sigma = objectives.std()
print("mu =",mu)
print("sigma =",sigma)

# Kolmogorov-Smirnov test
ks_statistic, p_value = stats.kstest(objectives, stats.norm.cdf, args=(mu,sigma))
# Set significance level
alpha = 0.05

print(f"K-S Statistic: {ks_statistic}")
print(f"P-value: {p_value}")

# Make a decision based on the p-value and significance level
if p_value < alpha:
    print("Reject the null hypothesis: The sample does not follow a normal distribution.")
else:
    print("Fail to reject the null hypothesis: The sample follows a normal distribution.")


# Plotting the histogram
plt.hist(objectives, bins=19, color='skyblue', edgecolor='black')
plt.xlabel('Optimal Objective Values')
plt.ylabel('Frequency')
plt.show()

# Plotting ECDF
objectives.sort()
ecdf = np.arange(1, len(objectives) + 1) / len(objectives)
plt.step(objectives, ecdf, where='post')
plt.xlabel('Optimal objective values')
plt.ylabel('Empirical Cumulative Distribution Function')
plt.show()


