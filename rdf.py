# script to calculate the radial distribution function (RDF) of a structure or a set of structures
# supports all file types supported by my structure class

# RDF: proper RDF analysis - constructs nxnxn supercell from the structure
#      and calculates proper RDF as defined on slide 11.10 of the atomistic modelling lecture notes (note that file is actually Lecture 10)
#      by iterating over all atoms in the central cell and averaging their individual distribution functions
#      (see equation 8 in https://scripts.iucr.org/cgi-bin/paper?S0021889800019993)
#      note that the normalisation by the density is not included in my implementation
#      the reason is that the density can differ between different structures in the case where we have several
#      and if there is only one, then normalising by the density just scales the RDF by a constant factor

# written by Pascal Salzbrenner, pts28@cam.ac.uk

import os
import re
import sys
import numpy as np
from operator import itemgetter
from structure import Structure
from exceptions import InputError

# read inputs
shell_width = float(sys.argv[1])
max_radius = float(sys.argv[2])
supercell_size = int(sys.argv[3])
file_type = sys.argv[4] # will read all files of this type in the current directory

half_shell_width = shell_width/2
num_decimals = len(str(half_shell_width)) - 2

# set up list of files to read
input_files = [file for file in os.listdir() if re.match(r".*\.{}".format(file_type), file)]

# set up a list of the different shells as well as a list of shell volumes
# the volume of a shell is 4*pi*shell_width*r**2, where r is taken to be the distance of the shell's centre from the origin
shell_radius = 0
shell_list = []
shell_volumes = []
# set up a list where each row will contain the shell occupations for one atom
shell_occupations =  []

# first bin value
shell_list.append(shell_radius)
shell_volumes.append(4*np.pi*shell_width*(shell_radius+half_shell_width)**2)

# generate bin values
while shell_radius <= max_radius:

    shell_radius += shell_width

    shell_list.append(shell_radius)
    shell_volumes.append(4*np.pi*shell_width*(shell_radius+half_shell_width)**2)

# the final value in shell_list will be the upper bound of the final bin (no radius larger than it will be included by definition)

shell_volumes = np.array(shell_volumes)

# iterate over all files
for input_file in input_files:

    structure = Structure(input_file)

    # construct supercell from the structure
    supercell_positions = []

    for i in range(-1*(supercell_size-2),supercell_size-1):
        for j in range(-1*(supercell_size-2),supercell_size-1):
            for k in range(-1*(supercell_size-2),supercell_size-1):
                for position in structure.positions_abs:
                    supercell_positions.append(position+np.dot(structure.lattice.T, np.array([i, j, k])))

    # loop over atoms in original cell and determine their distribution functions
    for position in structure.positions_abs:
        # set up array to contain the bin occupancies
        atom_shell_occupations = np.zeros(len(shell_list), dtype=int)

        # loop over every atom in the supercell - exclude the atom itself
        for supercell_position in supercell_positions:
            distance_vector = supercell_position - position
            distance = np.sqrt(distance_vector.dot(distance_vector))

            if np.isclose(distance, 0, atol=shell_width/10000):
                # if shell_width is reasonably chosen, no actual atom should be this close to the reference atom
                # thus this distance is 0 and represents the distance between the atom and itself
                continue

            for i in range(len(shell_list)-1):
                if distance >= shell_list[i] and distance < shell_list[i+1]:
                    atom_shell_occupations[i] += 1
                    break

        shell_occupations.append(atom_shell_occupations.copy())

shell_occupations = np.array(shell_occupations)
# average along the columns to get the elements for each equivalent shell over all arrays
average_shell_occupations = np.mean(shell_occupations, axis=0)
# normalise by the shell volumes
normalised_shell_occupations = average_shell_occupations/shell_volumes

# write output data
datafile = open("rdf.dat", "w")
datafile.write("# shell middle r [{}]; g(r)\n".format(structure.length_units))
for i in range(len(shell_list)-1):
    datafile.write("{0: <.{2}f} {1: >.5f}\n".format(shell_list[i]+half_shell_width, normalised_shell_occupations[i], num_decimals))
datafile.close()

# write plotting commands
plotfile = open("rdf.gnu", "w")
plotfile.write("set terminal postscript eps colour font 'Helvectica,20'\n")
plotfile.write("set output '| epstopdf --filter --outfile=rdf.pdf'\n")
plotfile.write("set mxtics 2\n")
plotfile.write("set mytics 2\n")
plotfile.write("set boxwidth {}\n".format(shell_width))
plotfile.write("set style fill solid\n")
plotfile.write("set xlabel 'r [{}]'\n".format(structure.length_units))
plotfile.write("set ylabel 'g(r)'\n")
plotfile.write("set xrange [0:{}]\n".format(max_radius))
plotfile.write("set yrange [0:]\n")
plotfile.write("plot 'rdf.dat' u 1:2 w boxes lc rgb '#DC143C' notitle")
plotfile.close()
