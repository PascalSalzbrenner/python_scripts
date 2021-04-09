# Script to calculate the difference of the electronic observable at every point from the mean at that point

# inputs:
# bs.i.j.dat - giving the values for the band j at k-point i for all N/2 configuration pairs, in order
# mean.i.j.dat - giving the mean for the electronic observable at band j at k-point i

# output:
# points_ordered.dat: the points, in ascending order according to their difference from the mean - the odd point of the pair is given

# run in directory containing a bs directory
# written by Pascal Salzbrenner, pts28

import sys
import numpy as np

# read in i and j for bs.i.j.dat and mean.i.j.dat
k_point = sys.argv[1]
lower_band_index = sys.argv[2]
higher_band_index = sys.argv[3]

points_observable = [] # list for the values of the observable at the points, in order

# read the values of the observable at the different points
with open ("bs_difference_{}_{}_{}.dat".format(k_point, lower_band_index, higher_band_index), "r") as points_file:
    for line in points_file:
        if not "#" in line:
            points_observable.append(float(line))

points_observable = np.array(points_observable)

# read the observable's mean
with open("renormalisation_{}_{}_{}.dat".format(k_point, lower_band_index, higher_band_index), "r") as mean_file:
    # read past the two comment lines
    for i in range(2):
            mean_file.readline()

    # read mean
    mean = float(mean_file.readline().split()[0])

# find the difference at every point
points_observable_difference = np.abs(points_observable-mean)

# sort in ascending order
points_observable_difference_sorted = np.sort(points_observable_difference)

order_list = [] # list of configuration indexes, starting from that with the lowest difference from the mean

for i in range(len(points_observable_difference_sorted)):
    # if there are several points with equivalent distance from the mean, they will all be found
    # treat them all the first time they are found, and skip points that have the same difference as their predecessor
    if i != 0 and np.isclose(points_observable_difference_sorted[i], points_observable_difference_sorted[i-1]):
        continue

    indices = np.where(np.isclose(points_observable_difference_sorted[i], points_observable_difference))[0]

    for index in indices:
        order_list.append(str(2*index+1)) # such that always the odd (ie lower) element of the pair is given

with open("points_order.dat", "w") as outfile:
    outfile.write("\n".join(order_list))
