# script to calculate the phonon zero-point energies from the eigenvalues of the phonon modes as returned from i-pi
# the frequencies in Hartree are the square roots of the eigenvalues
# from the frequencies we can calculate the ZPE via equation 15 in https://iopscience.iop.org/article/10.1088/1361-648X/aaa737/meta
# for the ZPE, s=0 for all modes

import sys
import numpy as np

# define conversion factor from Hartree to meV
conversion_factor = 27211.3825435

filename = sys.argv[1] # the name of the file containing the eigenvalues
num_atoms = float(sys.argv[2]) # the number of atoms in the unit cell

infile = open("{}".format(filename), "r")

# skip comment line
infile.readline()

# sum over all frequencies
frequency_sum = 0

for line in infile:

    eigval = float(line.split()[0])

    # the Gamma point acoustic frequencies, which are 0, can be numerically slightly negative
    # to avoid an imaginary outcome to the square root, we set them to exactly 0 here
    if np.isclose(eigval, 0, atol=1e-12):
        eigval = 0

    frequency_sum += float(line.split()[0])**0.5

zpe = 0.5 * conversion_factor * frequency_sum / num_atoms

with open("zpe.txt", "w") as outfile:
    outfile.write("# ZPE [meV/atom]\n")
    outfile.write("{}\n".format(zpe))
