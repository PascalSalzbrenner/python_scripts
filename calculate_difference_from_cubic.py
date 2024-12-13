# script to calculate the absolute value of the average difference of all angles from the cubic angle of 90°

import sys
import numpy as np
import ase.neighborlist
import matplotlib.pyplot as plt

from ase.io import read

filename = sys.argv[1]
# all atoms with a distance less than this from the central atom will be considered its first-nearest-neighbours:
nn_dist_max = float(sys.argv[2])

time_step = None

# time step, if given, must be the actual step between images, in ps
if len(sys.argv) > 3:
    time_step = float(sys.argv[3])

# set up list to store the angle
angle_diff = []

# read trajectory
trajectory = read(filename, index=":", format="extxyz")

# define neighbourlist
nl = ase.neighborlist.NeighborList(nn_dist_max, self_interaction = False,
                                   primitive=ase.neighborlist.NewPrimitiveNeighborList)

# iterate over every frame
for i in range(len(trajectory)):

    print("Frame {} out of {}".format(i+1, len(trajectory)))

    nl.update(trajectory[i])

    # list to store angles:
    angles = []

    # iterate over all atoms
    for atom_index in range(len(trajectory[i])):

        atom_nl = nl.get_neighbors(atom_index)[0]
        # iterate over all neighbours for the atom
        for j in range(len(atom_nl)-1):
            for k in range(j+1, len(atom_nl)):
                angle = trajectory[i].get_angle(atom_nl[j], atom_index, atom_nl[k], mic=True)

                # only consider angles which can lie within the trigonal deformation path
                if 60 <= angle < 110:
                    angles.append(abs(90-angle))

    angle_diff.append(np.mean(angles))

# write data
with open("angle_diff.dat", "w") as f:
    if not time_step:
        f.write("# MD Step; Average difference from 90° [°]\n")
        multiply_factor = 1
    else:
        f.write("# Time [ps]; Average difference from 90° [°]\n")
        multiply_factor = time_step

    for i, angle in enumerate(angle_diff):
        f.write("{} {}\n".format(multiply_factor*i, angle))

# plotting
if not time_step:
    plt.plot(angle_diff, color="#E6AB02", linestyle="solid")
    plt.xlabel("MD Step")
else:
    times = np.arange(len(angle_diff)) * time_step
    plt.plot(times, angle_diff, color="#E6AB02", linestyle="solid")
    plt.xlabel("Time [ps]")

plt.ylabel("Angle Difference from 90°")
plt.ylim(bottom=0)
plt.minorticks_on()
plt.savefig("angle_diff.png", dpi=300)
