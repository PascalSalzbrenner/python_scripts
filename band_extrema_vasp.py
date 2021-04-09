#postprocessing utility to easily read the highest/lowest value of a specified band from VASP EIGENVAL output of a band structure calculation
#written by Pascal Salzbrenner - pts28@cam.ac.uk

import os
import numpy as np

band_index = input("Please enter the index of the band whose maximum and minimum you want to evaluate: ")

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

#read past the six header lines

line = eigenval.readline()

while line.split():
    line = eigenval.readline()
    continue

#the next line gives the direct coordinates of the first k-point in its first three elements

line = eigenval.readline()
k_point = line.split()

#the k-points associated with the highest and lowest energies of the band
k_point_E_max = k_point[0:3]
k_point_E_min = k_point[0:3]

#the first k_point is treated differently, because we set the energy baselines

while line.split():

    line_list = line.split()

    if line_list[0] == band_index:
        E_max = float(line_list[1])
        E_min = float(line_list[1])

    line = eigenval.readline()

#read the next k-point

line = eigenval.readline()
k_point = line.split()

#loop through the rest of the file

for l in eigenval:

    #if the line is empty, we know the next line is a k-point, so we go straight to that
    if not l.split():
        line = eigenval.readline()
        k_point = line.split()

    elif l.split()[0] == band_index:
        line_list = l.split()
        if float(line_list[1]) > E_max:
            E_max = float(line_list[1])
            k_point_E_max = k_point[0:3]
        elif float(line_list[1]) < E_min:
            E_min = float(line_list[1])
            k_point_E_min = k_point[0:3]

eigenval.close()

#write output

with open("band_{}_max_min.dat".format(band_index), "w") as outfile:
    outfile.write("The band's highest energy is {} eV. In units of the reciprocal lattice vectors, it is at: {} {} {}\n".format(E_max, k_point_E_max[0],
    k_point_E_max[1], k_point_E_max[2]))
    outfile.write("The band's lowest energy is {} eV. In units of the reciprocal lattice vectors, it is at: {} {} {}\n".format(E_min, k_point_E_min[0],
    k_point_E_min[1], k_point_E_min[2]))
