# code to generate a list of q-points in the IBZ as required as input by the Voronoi tiling module of Caesar
# detects and removes duplicates

# input: paths to directories containing the desired Caesar-generated ibz.dat files

import sys
import numpy as np

from collections import Counter

def array_to_string(array):
    # convert an array into a neat string. Eg if the array is [1,2,3], the string will be "1 2 3"
    raw_string=str(array)
    nice_string=raw_string.lstrip("[")
    nice_string=nice_string.rstrip("]")

    return nice_string

# set up list to ultimately contain all points from all grid sizes
points = []

# set up list to contain multiplicity of q-point in the same order
multiplicities = []

# read in all input files
for path in sys.argv[1:]:
    ibz_file = open("{}/ibz.dat".format(path.rstrip("/")), "r")

    for line in ibz_file:
        split_line=line.split()
        point = np.array(split_line[:-1], dtype=np.float64)

        # check if a certain element is already in the list

        if (not points) or (not len(point) in list(Counter(np.where(np.isclose(point, points))[0]).values())):
            # first iteration, where no points have been added yet, or
            # if any element in the first list returned by np.where is present D (number of dimensions) number of times, that means all
            # elements of point have been found in a single element in points - ie point is already in it. Only add it if this isn't the case

            points.append(point)
            multiplicities.append(split_line[-1])

# convert points to strings for writing
string_points = []

for point in points:
    string_points.append(array_to_string(point))

# write output
with open("multiplicities.dat", "w") as outfile:
    outfile.write("\n".join(multiplicities))

with open("q_points.dat", "w") as outfile:
    outfile.write("#q-points from calculations in the following directories: {}\n".format("; ".join(sys.argv[1:])))
    outfile.write("\n".join(string_points))
