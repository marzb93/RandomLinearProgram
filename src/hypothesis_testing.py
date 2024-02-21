
import numpy as np
from scipy.stats import kstest, norm
import matplotlib.pyplot as plt
from scipy import stats
from functools import partial

file_path = '../build/obj.txt' 
with open(file_path, 'r') as file:
    numbers_str = file.read().strip()

# Convert the string of numbers to a list of floats
objectives = np.array( [float(num) for num in numbers_str.split(', ')] )
mu = objectives.mean()
sigma = objectives.std()
print("mu =",mu)
print("sigma =",sigma)


# Kolmogorov-Smirnov test
#shifting mean, and scaling by sigma to get standard normal 
# normed_objectives = [x/sigma for x in [obj-mu for obj in objectives]]
# ks_statistic, p_value = stats.kstest(normed_objectives, stats.norm.cdf)
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
#plt.title('Histogram of Optimal Objective Values')
#plt.grid(True)
plt.show()

ecdf = np.arange(1, len(objectives) + 1) / len(objectives)
# Plotting ECDF
plt.step(objectives, ecdf, where='post')
plt.xlabel('Optimal objective values')
plt.ylabel('Empirical Cumulative Distribution Function')
# plt.title('Empirical Cumulative Distribution Function')
plt.show()


"""
loc, scale = stats.gumbel_r.fit(objectives)
gumbel_dist = stats.gumbel_r(loc = loc, scale=scale )
# ks_statistic, p_value = stats.kstest(objectives, gumbel_dist.cdf )
x = np.linspace(gumbel_dist.ppf(0.01), gumbel_dist.ppf(0.99), 100)  
pdf = gumbel_dist.pdf(x)

# Plot the Gumbel distribution
plt.figure(figsize=(8, 6))
plt.hist(objectives, density=True, bins=19, alpha=0.5, label='Data Histogram')  
plt.plot(x, pdf, label='Gumbel Distribution PDF')
plt.title('Gumbel Distribution')
plt.xlabel('Values')
plt.ylabel('Probability Density')
plt.legend()
plt.grid(True)
plt.show()

"""
