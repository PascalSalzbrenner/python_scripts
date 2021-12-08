# Python script to calculate the number density from any filetype for which a parser has been defined in my structure.py module

import sys
from structure import Structure

filename = sys.argv[1]

structure = Structure(filename)

with open("volume_density.txt", "w") as outfile:
    outfile.write("# Input file: {}\n".format(filename))
    outfile.write("# Number of atoms; Volume [{0}**3]; Density [N/{0}**3]\n".format(structure.length_units))
    outfile.write("{} {} {}".format(structure.num_atoms, structure.volume, structure.density))
