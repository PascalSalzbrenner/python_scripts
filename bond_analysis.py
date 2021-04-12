# script to do various kinds of analysis and plots of the bonds in a crystal
# currently works using the bond analysis in a .castep output file of a CASTEP run
# I could write parsers for the structure input files of different DFT codes and then do the length / RDF type analyses for them

# Features:
# PDF: bins all bonds, divides by the number of atoms in the cell, and plots a histogram of frequency vs length
#      if a structure is given, essentially repeats the analysis of the .castep file
#      that is, PDF analysis is confined to unit cell
# RDF: proper RDF analysis - constructs 9x9x9 supercell from the structure (which can be read from the .castep or structure file)
#      and calculates proper RDF as defined on slide 11.10 of the atomistic modelling lecture notes
#      by iterating over all atoms in the central cell and averaging their individual distribution functions
#      (see equation 8 in https://scripts.iucr.org/cgi-bin/paper?S0021889800019993)
# bond length analysis: Takes the indices (bonds ranked by length from shortest to longest) of an arbitrary number of bonds and tracks their
#                       length(s) as a function of pressure
# bond population analysis: Likewise, but for the population. Only works from a .castep file

# if one of the latter two tasks is selected, the code reads all relevant files in the directory, and the pressure is determined from them

# written by Pascal Salzbrenner, pts28

import os
from operator import itemgetter
from exceptions import InputError
from utilities import lattice_basis_to_cartesian, get_structure_from_castep, get_structure_from_cell

# define two functions to determine the bonds, either from a structure, or from the CASTEP file, where they are already present and must
# only be read
def get_bonds_from_castep(filename):
    """Function to read the data about the bonds directly from the .castep file, where a Mulliken bond analysis is done
    :param str filename: the name of the .castep file

    :returns list bonds: a list in the format [("atom_name_number_1-atom_name_number_2", length, population)], sorted by bond length"""

    # set up output list
    bonds = []

    bonds_file = open("{}".format(filename), "r")

    for line in bonds_file:
        if "Bond" in line:
            # read past the next line; each line after that until a line full of "=" contains a bond
            bonds_file.readline()

            for bonds_line in bonds_file:
                if "===" in bonds_line:
                    break
                else:
                    bonds_data = bonds_line.split()
                    bonds.append(("{}{} - {}{}".format(bonds_data[0], bonds_data[1], bonds_data[3], bonds_data[4]), (float(bonds_data[-1]),
                                  float(bonds_data[-2]))))

    bonds_file.close()

    # sort by bond length
    bonds.sort(key=itemgetter(1))

    return bonds

# get necessary input
task = input("Which task would you like to run? (PDF, RDF, bond_length, bond_population - see the code header for descriptions) ").lower()

if task == "rdf" or task == "pdf":
    # different input files are possible
    input_file = input("What is the name of the input file? (currently supported: .castep, .cell files) ")
    input_file_type = input_file.split(".")[-1]
elif task == "bond_length":
    input_file_type = input("What is the type of input file? (currently supported: castep, cell) ").lower()
elif task == "bond_population":
    # only possile from a .castep file
    input_file_type = "castep"
else:
    # a non-implemented task has been requested
    raise InputError("task",
    "You have requested a task which is not implemented. Task must be one of the following: PDF, RDF, bond_length, bond_population.")


if input_file_type=="cell":
    get_structure_from_cell

if task == "pdf":

    if input_file_type=="castep":
        get_bonds_from_castep
    else:
        # any structure format will already have been parsed into the universal internal format
        get_bonds_from_structure

elif task == "rdf":

    if input_file_type=="castep":
        get_structure_from_castep

    # construct supercell from the structure, loop over atoms in original cell
elif task == "bond_length" or task == "bond_population":

    # for this and bond_population, loop over all files of the right type in the directory

    if task == "bond_length":
        data_index = 0 # the index which we will read from the dictionary values

        if input_file_type=="castep":
            get_bonds_from_castep
        else:
            get_bonds_from_structure
    else:
        data_index = 1

        if input_file_type=="castep":
            get_bonds_from_castep
        else:
            raise InputError("Bond population analysis",
            "You have requested a bond population plot, but this is only possible with a .castep file as input.")

elif task == "bond_population":

    get_bonds_from_castep

# in the relevant if blocks, do the analysis
