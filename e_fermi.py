#postprocessing utility to easily read the energy of the highest occupied state at the band edge
#written by Pascal Salzbrenner - pts28@cam.ac.uk

import os
import numpy as np

# before reading, stitch together one big EIGENVAL file containing the entire path if the calculations have been split up

# determine which (if any) splits exist

ls = os.listdir()
ls_temp = []

for el in ls:
    if el.startswith("split-"):
        ls_temp.append(el)

if ls_temp:
    # ls contains one or more elements, ie the band structure calculations have been split up
    ls = ls_temp[:]
    ls.sort()

    # combine EIGENVAL files in different directories into one EIGENVAL file
    full_eigenval = open("EIGENVAL", "w")

    for split_dir in ls:
        with open("{}/EIGENVAL".format(split_dir), "r") as split_eigenval:

            initial_line_counter = 0

            for line in split_eigenval:

                if initial_line_counter < 6:
                    if split_dir == "split-01":
                        # write the header (because it is read past in the band gap calculation below, so if the file doesn't have a header this might cause problems)
                        full_eigenval.write(line)
                        initial_line_counter += 1
                        continue

                    else:
                        # read past the header and ignore it
                        initial_line_counter += 1
                        continue

                if not line.split():
                    # this line is empty, so the next line is a k-point plus weight
                    k_point = split_eigenval.readline()
                    weight = float(k_point.split()[-1])

                    if np.isclose(weight, 0):
                        write_line = True # only the zero-weight k-points (those that are used in the band structure plot) are written
                        full_eigenval.write("\n")
                        full_eigenval.write(k_point)
                    else:
                        write_line = False

                else:
                    if write_line:
                        full_eigenval.write(line)
                    else:
                        continue

    full_eigenval.close()

eigenval = open("EIGENVAL", "r")

# read the full EIGENVAL file and calculate the band gap
#read past the six header lines
for line in eigenval:
    if not line.split():
        break

#the next line gives the direct coordinates of the first k-point in its first three elements

line = eigenval.readline()
k_point = line.split()

#the k-points associated with the highest energy of an occupied band, and the lowest energy of an unoccupied band
k_point_E_max = k_point[0:3]

#the first k_point is treated differently, because we use it to determine the highest occupied band, and to set the baselines

is_occupied = True
line = eigenval.readline()

while line.split():

    if is_occupied:
        if np.isclose(float(line.split()[2]), 1):
            max_occupied = line.split()[0]
            E_max_occupied = float(line.split()[1])
        else:
            is_occupied = False

    line = eigenval.readline()

#read the next k-point

line = eigenval.readline()
k_point = line.split()

#loop through the rest of the file, knowing the signatures of the highest occupied band

for l in eigenval:

    #if the line is empty, we know the next line is a k-point, so we go straight to that
    if not l.split():
        line = eigenval.readline()
        k_point = line.split()

    elif l.split()[0] == max_occupied:
        if float(l.split()[1]) > E_max_occupied:
            E_max_occupied = float(l.split()[1])
            k_point_E_max = k_point[0:3]

eigenval.close()

#write out Fermi energy
with open("fermi_energy.dat", "w") as outfile:
    outfile.write("Fermi energy = {} [eV]\n".format(E_max_occupied))
    outfile.write("In units of the reciprocal lattice vectors, the valence band edge is at: {} {} {}\n".format(k_point_E_max[0],
    k_point_E_max[1], k_point_E_max[2]))
