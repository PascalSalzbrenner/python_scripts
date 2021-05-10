# script to take the output of a Caesar phonon calculation and add the zero point effect (T=0 value in interpolated_free_energy.dat) to the
# enthalpy in the corresponding .res file
# written by Pascal Salzbrenner, pts28@cam.ac.uk

import sys

# parse input
path = sys.argv[1].rstrip("/") # the path to the directory containing the original .res file
# the .res file naming convention is <seed>_<pressure>p0-<filename_rest>.res
seed = sys.argv[2]
filename_rest = sys.argv[3]
pressure = sys.argv[4]

# generate path to the parent directory of the .res file directory
parent_path = "/".join(path.rstrip("/").split("/")[:-1])

# generate .res filename
res_filename = "{}_{}p0-{}.res".format(seed, pressure, filename_rest)

# read the ZPE - first line of the file - the second column is in Hartree, the third in eV
with open("lte/interpolated_free_energy.dat", "r") as zpe_file:
        zpe = float(zpe_file.readline().split()[2])

# open input and output files
infile = open("{}/{}".format(path, res_filename), "r")
outfile = open("{}/enthalpy_plot_ZPE/{}".format(parent_path, res_filename), "w")

# first line of the input res file is the only one requiring special handling
# the total enthalpy (enthalpy per cell) is the 5th element of the first line

# note that the static lattice enthalpy calculations are sometimes done in the conventional cell, which makes geometry relaxations easier
# phonon calculations, on the other hand, are almost always done in the primitive cell
# thus, it is possible that the number of atoms differs
# In order to get the right energy correction per cell, the Caesar ZPE is multiplied by n_atoms_static / n_atoms_phonon
# n_atoms_static is the 8th element of the .res file's first line

res_data = infile.readline().split()
static_enthalpy = float(res_data[4])
n_atoms_static = float(res_data[7])

# read n_atoms_phonon from the structure.dat file
structure_file = open("structure.dat", "r")

n_atoms_phonon = 0

for line in structure_file:

    if line.lower().startswith("atoms"):
        # next line is the first line containing an atom - we count the number of lines to get the number of atoms

        for atoms_line in structure_file:

            if line.lower().startswith("symmetry"):
                # the last line containing an atom has been passed
                break
            else:
                n_atoms_phonon += 1

        break

structure_file.close()

atoms_factor = n_atoms_static / n_atoms_phonon

res_data[4] = str(static_enthalpy + zpe)

outfile.write(" ".join(res_data))
outfile.write("\n")

# rest of the lines are just copied over
for line in infile:
    outfile.write(line)

infile.close()
outfile.close()
