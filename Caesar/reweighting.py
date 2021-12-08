# little program to test the reweighting of the integral F(T) = int f(x)g(x, T)
# f(x) is a simple function, g(x,T) is a Gaussian whose width changes as a function of temperature T

# Gaussian is that of Monserrat - Electron-phonon coupling from finite differences, with the frequency set to 1, and in units where k_B = 1


import math
import random
import numpy as np
from mpmath import coth, inf

def sigma(T):
    if T == 0:
        return np.sqrt(0.5)
    else:
        return float(np.sqrt(0.5*coth(1/(2*T))))

def gaussian(x, T):
    if T == 0:
        return float(math.exp(-x**2/(coth(inf)))/(np.sqrt(np.pi*coth(inf)))) # coth(inf) = 1
    else:
        return float(math.exp(-x**2/(coth(1/(2*T))))/(np.sqrt(np.pi*coth(1/(2*T)))))

def f(l):
    """returns an array of f(x), for all x in list l"""

    out_list = []
    for x in l:
        out_list.append(x**3+2*x**2-3*x)

    return np.array(out_list)

num_samples = int(input("How many points would you like to generate? "))

# generate num_samples values according to a Gaussian distribution at 50 K
points_300K = []
for i in range(num_samples):
    points_300K.append(random.gauss(0, sigma(300)))

f_x_300K = f(points_300K)

# estimate the integral at 300 K
F_300K = np.mean(f_x_300K)

# estimate the uncertainty of F_300K
squared_differences_300K = (f_x_300K - F_300K)**2
uncertainty_300K = np.sqrt(np.sum(squared_differences_300K)/(num_samples*(num_samples-1)))

print("F(300): {}".format(F_300K))
print("Uncertainty in estimate: {}".format(uncertainty_300K))

# generate num_samples values according to a Gaussian distribution at 50 K
points_50K = []
for i in range(num_samples):
    points_50K.append(random.gauss(0, sigma(50)))

f_x_50K = f(points_50K)

# estimate the integral at 50 K
F_50K = np.mean(f_x_50K)

# estimate the uncertainty of F_50K
squared_differences_50K = (f_x_50K - F_50K)**2
uncertainty_50K = np.sqrt(np.sum(squared_differences_50K)/(num_samples*(num_samples-1)))

print("F(50): {}".format(F_50K))
print("Uncertainty in estimate: {}".format(uncertainty_50K))

# estimate the integral at 300K using the points at 50K
reweighted_points_300K = []

for i in range(len(f_x_50K)):
    reweighted_points_300K.append(f_x_50K[i] * gaussian(points_50K[i], 300) / gaussian(points_50K[i], 50))

reweighted_points_300K = np.array(reweighted_points_300K)
reweighted_F_300K = np.mean(reweighted_points_300K)

# estimate the uncertainty of reweighted_F_300K
reweighted_squared_differences_300K = (reweighted_points_300K - reweighted_F_300K)**2
reweighted_uncertainty_300K = np.sqrt(np.sum(reweighted_squared_differences_300K)/(num_samples*(num_samples-1)))

print("F(300), estimated by reweighting the points generated according to the Gaussian distribution at 50K: {}".format(reweighted_F_300K))
print("Uncertainty in this reweighted estimate for F(300): {}".format(reweighted_uncertainty_300K))

# estimate the integral at 50K using the points at 300K
reweighted_points_50K = []

for i in range(len(f_x_300K)):
    reweighted_points_50K.append(f_x_300K[i] * gaussian(points_300K[i], 50) / gaussian(points_300K[i], 300))

reweighted_points_50K = np.array(reweighted_points_50K)
reweighted_F_50K = np.mean(reweighted_points_50K)

# estimate the uncertainty of reweighted_F_50K
reweighted_squared_differences_50K = (reweighted_points_50K - reweighted_F_50K)**2
reweighted_uncertainty_50K = np.sqrt(np.sum(reweighted_squared_differences_50K)/(num_samples*(num_samples-1)))

print("F(50), estimated by reweighting the points generated according to the Gaussian distribution at 300K: {}".format(reweighted_F_50K))
print("Uncertainty in this reweighted estimate for F(50): {}".format(reweighted_uncertainty_50K))
