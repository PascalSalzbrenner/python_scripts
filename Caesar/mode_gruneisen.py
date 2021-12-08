# script to estimate the Grüneisen parameters for every mode given phonon calculations for two volumes
# Grüneisen parameters defined as per Ashcroft and Mermin, Solid State Physics, equation (25.18)
# written by Pascal Salzbrenner, pts28

import numpy as np

lattice_param_1 = input("Give the first lattice parameter at which the frequencies are read: ")
lattice_param_2 = input("Give the second lattice parameter at which the frequencies are read: ")

# set up lists for different parameters which will be read in
q_points = []
branch_indices = []

frequencies_1 = []
frequencies_2 = []

# determine the two volumes by reading the second, third, and fourth lines of Caesar structure.dat file
# the format of this is very strict, so these lines will always contain the lattice parameter
volumes = []

for lattice_param in [lattice_param_1, lattice_param_2]:
    lattice_vectors = []
    with open("lattice_parameter_{}/structure.dat".format(lattice_param), "r") as structfile:
        # skip first line
        structfile.readline()

        for i in range(3):
            lattice_vectors.append(np.array(structfile.readline().split(), dtype=float))

        basis = np.array(lattice_vectors)
        volumes.append(np.linalg.det(basis))

        # determine reciprocal lattice
        if lattice_param == lattice_param_1:
            reciprocal_lattice = 2 * np.pi * np.linalg.inv(basis).T

        # determine number of atoms
        # read past line containing the "Atoms" tag
        structfile.readline()
        no_atoms = 0
        for line in structfile:
            if not line.startswith("Symmetry"):
                # all lines preceding this contain an atom
                no_atoms += 1
                continue
            else:
                break

volumes_log_difference = np.log(volumes[1]) - np.log(volumes[0])

# there are three branches per atom
no_branches = 3 * no_atoms

branch_index = 1

# read frequencies at first volume
disp_patterns_1 = open("lattice_parameter_{}/lte/disp_patterns.dat".format(lattice_param_1), "r")

for line in disp_patterns_1:
    if line.lstrip().startswith("Frequency"):
        # line contains a Frequency
        # also means the next line contains a q-point
        frequencies_1.append(float(line.split()[2]))
        q_points.append(np.array(disp_patterns_1.readline().split(), dtype=float))
        branch_indices.append(branch_index)
        if branch_index == no_branches:
            # reset branch index, the next frequency will belong to a different q-point
            branch_index=1
        else:
            branch_index += 1
    else:
        continue

disp_patterns_1.close()

# read frequencies at second volume
disp_patterns_2 = open("lattice_parameter_{}/lte/disp_patterns.dat".format(lattice_param_2), "r")

for line in disp_patterns_2:
    if line.lstrip().startswith("Frequency"):
        # line contains a Frequency
        frequencies_2.append(float(line.split()[2]))
        # q-points and branch indices only have to be read in once
    else:
        continue

disp_patterns_2.close()

# calculate logarithmic difference between frequencies
frequencies_log_difference = []

for i in range(len(frequencies_1)):
    if frequencies_1[i] > 0 and frequencies_2[i] > 0:
        # non-zero frequencies
        frequencies_log_difference.append(np.log(frequencies_2[i])-np.log(frequencies_1[i]))
    else:
        # the frequencies are 0 - this will be constant regardless of the volume, so the difference is set to 0
        frequencies_log_difference.append(0)

frequencies_log_difference = np.array(frequencies_log_difference)

# calculate Gruneisen parameters
gruneisen_parameters = - frequencies_log_difference/volumes_log_difference

# write output
outlines = ["frequency [Hartree], q-point [reciprocal coordinates], branch index, Gruneisen parameter"]

for i in range(len(gruneisen_parameters)):
    outlines.append("{} {} {} {}".format(frequencies_1[i], np.dot(np.linalg.inv(reciprocal_lattice.T), q_points[i]), branch_indices[i],
                    gruneisen_parameters[i]))

with open("grueneisen_parameters.dat", "w") as outfile:
    outfile.write("\n".join(outlines))
