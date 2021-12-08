#postprocessing utility to easily calculate the value of a specific band at a specific k-point
#written by Pascal Salzbrenner - pts28@cam.ac.uk

import os
import sys
import numpy as np

band_index = sys.argv[1] # the index of the band whose energy we want to evaluate
k_point = np.array(sys.argv[2].split(), dtype=np.double) # The k-point at which you want to evaluate the band energy, in the format kx ky kz

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
                    current_k_point = split_eigenval.readline()
                    weight = float(current_k_point.split()[-1])

                    if np.isclose(weight, 0):
                        write_line = True # only the zero-weight k-points (those that are used in the band structure plot) are written
                        full_eigenval.write("\n")
                        full_eigenval.write(current_k_point)
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
current_k_point = np.array(line.split()[0:3], dtype=np.double)

for l in eigenval:

    line_list = l.split()

    #if the line is empty, we know the next line is a k-point, so we go straight to that
    if not line_list:
        line = eigenval.readline()
        current_k_point = np.array(line.split()[0:3], dtype=np.double)
    
    elif line_list[0] == band_index and np.isclose(k_point, current_k_point).all():
        energy = line_list[1]

eigenval.close()

#calculate gap energy, and determine whether the band gap is direct or indirect

with open("band_{}_k_{}_{}_{}_energy.dat".format(band_index, str(k_point[0]), str(k_point[1]), str(k_point[2])), "w") as outfile:
    outfile.write("Energy = {} eV\n".format(energy))
