# visualise in the vibrational BZ the q-points connecting points k and k' as given in the ep_kpoint_pairs block in CASTEP .cell files
# written by Pascal Salzbrenner, pts28@cam.ac.uk

import sys
import numpy as np
import matplotlib.pyplot as plt

from itertools import product, combinations

# read fileroot of .cell file
fileroot = sys.argv[1]

# open cell file for reading
cell_file = open("{}.cell".format(fileroot), "r")

# initialise list of q-points
q_points = []

# read the lines in the file, skipping all that are unrelated to the ep_kpoint_pairs
for line in cell_file:
    if not "ep_kpoint_pairs" in line.lower():
        continue
    else:
        # we have found the beginning of the ep_kpoint_pairs block
        # all following lines (until the next line containing "ep_kpoint_pairs", ie the endblock statement) contain a pair of k and k'
        # format: k_x k_y k_z k'_x k'_y k'_z

        for kpoint_pair in cell_file:
            if "ep_kpoint_pairs" in kpoint_pair.lower():
                # end of the block with relevant information
                break
            else:
                # data point
                k_k_dash = kpoint_pair.split()
                k = np.array(k_k_dash[0:3], dtype=float)
                k_dash = np.array(k_k_dash[3:6], dtype=float)
                q = k_dash - k
                q_points.append(q)

    break

# map q-points to the first BZ
for i in range(len(q_points)):
    greater_point_five_array = np.array((np.abs(q_points[i]) > 0.5), dtype=int)
    q_points[i] = q_points[i] - np.sign(q_points[i])*greater_point_five_array

q_points = np.array(q_points)


# plot the output together with the BZ
# as we are in fractional coordinates, the BZ will just be a cube stretching from -0.5 to 0.5 in each dimension
fig = plt.figure()
ax = fig.add_subplot(projection="3d")

# plot the BZ
edges = [-0.5, 0.5]
for s, e in combinations(np.array(list(product(edges, edges, edges))), 2):
   if np.sum(np.abs(s-e)) == edges[1]-edges[0]:
      ax.plot3D(*zip(s, e))

# plot the q-points

# scatter requires lists / arrays containing the x, y, and z coordinates respectively so we have to access the array columns
ax.scatter(q_points[:,0], q_points[:,1], q_points[:,2])

plt.show()
