# Python script to calculate the volume and density of a unit cell from a .xyz file giving the structure
# this assumes that the lattice is given in abcABC format
# I believe there are other formats but am overall unfamiliar with the .xyz format
# so I am sticking with this for now
# it is always possible - and should be easy - to extend this later

import sys
import numpy as np
from utilities import construct_lattice_from_abc

# read filename
filename = sys.argv[1]

# open input file and read the required data from the first and second lines
xyz_file = open("{}".format(filename), "r")
num_atoms = int(xyz_file.readline().split()[0])
abc_line = xyz_file.readline().split()
xyz_file.close()

vector_lengths = [float(length) for length in abc_line[2:5]]
vector_angles = [float(angle) for angle in abc_line[5:8]]

lattice = construct_lattice_from_abc(vector_lengths, vector_angles)

volume = np.linalg.det(np.dstack(lattice))[0]
density = num_atoms / volume

with open("volume_density.txt", "w") as outfile:
    outfile.write("# Input file: {}\n".format(filename))
    outfile.write("# Number of atoms; Volume [A**3]; Density [N/A**3]\n")
    outfile.write("{} {} {}".format(num_atoms, volume, density))
