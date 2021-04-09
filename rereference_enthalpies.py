# script to change the reference for the relative enthalpies given in an input file in the format

# header
#
# # structure name
# pressure enthalpy
# ...

# assumes that the blocks for each structure contain the same number of data points

# written by Pascal Salzbrenner, pts28@cam.ac.uk

import numpy as np

# energy units are assumed meV/atom, but this can easily be made an input variable
units = "meV/atom"

# read filename, current reference structure, new reference structure
filename = input("What is the name of the file containing the enthalpies? ")
current_reference = input("What is the name of the current reference structure for the enthalpies? ")
new_reference = input("What is the name of the structure you want the enthalpies to be measured relative to? ")

filename_components = filename.split(".")

# set up data containers
# will contain the data in the form "structure_name": [[pressures], [enthalpies]]
enthalpies_in_dict = {}
enthalpies_out_dict = {}

# read enthalpies
enthalpies_file = open("{}".format(filename), "r")

# read past header
while enthalpies_file.readline().lstrip().startswith("#"):
    continue

# when this is done, we have reached the empty line indicating the end of the header - first structure name in the next line
for line in enthalpies_file:
    if len(line.split()) == 0:
        # empty line
        continue
    elif line.lstrip().startswith("#"):
        # structure name
        structure_name = line.lstrip().lstrip("#").lstrip().rstrip("\n")

        enthalpies_in_dict[structure_name] = [[],[]]
    else:
        # data point
        data = line.split()
        enthalpies_in_dict[structure_name][0].append(float(data[0]))
        enthalpies_in_dict[structure_name][1].append(float(data[1]))

enthalpies_file.close()

# the enthalpies for the old references are -1 * the old enthalpies for the new reference
enthalpies_out_dict[current_reference] = [enthalpies_in_dict[new_reference][0], -1*np.array(enthalpies_in_dict[new_reference][1])]

# reference
for key, value in enthalpies_in_dict.items():

    enthalpies_out_dict[key] = [enthalpies_in_dict[key][0],
                                np.array(enthalpies_in_dict[key][1])-np.array(enthalpies_in_dict[new_reference][1])]

# write out data
outfile = open("{}_{}_reference.{}".format(filename_components[0], new_reference, filename_components[1]), "w")

outfile.write("# enthalpies relative to {}\n".format(new_reference))
outfile.write("# Pressure [GPa]; relative enthalpy [{}]\n".format(units))
outfile.write("\n")

for key, value in enthalpies_out_dict.items():
    outfile.write("# {}\n".format(key))

    for i in range(len(value[0])):
        outfile.write("{} {}\n".format(value[0][i], value[1][i]))

    outfile.write("\n")

outfile.close()
