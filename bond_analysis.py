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
# bonds: simply prints out data about the bonds in the bonds.dat file

# if one of the latter two tasks is selected, the code reads all relevant files in the directory, and the pressure is determined from them

# written by Pascal Salzbrenner, pts28

import os
from operator import itemgetter
from structure import Structure
from exceptions import InputError

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

if task == "rdf" or task == "pdf" or task == "bonds":
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

if task == "pdf":

    if input_file_type=="castep":
        bonds = get_bonds_from_castep(input_file)
    else:
        structure = Structure(input_file)
        structure.get_bonds()
        bonds = structure.bonds

elif task == "rdf":

    structure = Structure(input_file)

    # construct supercell from the structure, loop over atoms in original cell

elif task == "bond_length" or task == "bond_population":

    # for these two tasks, loop over all files of the right type in the directory
    # determine all those files

    # set up container
    analysis_files = []
    for item in os.listdir():
        if "input_file_type" in item:
            analysis_files.append(item)

elif task == "bonds":

    if input_file_type == "castep":
        bonds = get_bonds_from_castep(input_file)
        length_units = "Angstrom"
    else:
        structure = Structure(input_file)
        structure.get_bonds()
        bonds = structure.bonds
        length_units = structure.length_units

    outfile = open("bonds.dat", "w")

    first_line_str = "# Bond; Length [{}]; ".format(length_units)

    if input_file_type == "castep":
        # write out bond population as well
        first_line_str += "Population; bond number\n"
    else:
        # only write the bond number
        first_line_str += "bond number\n"

    outfile.write(first_line_str)

    for i in range(len(bonds)):

        write_str = "{} {} ".format(bonds[i][0], bonds[i][1])

        if input_file_type == "castep":
            # write out bond population as well
            write_str += "{} {}\n".format(bonds[i][2], i+1)

        outfile.write(write_str)

    outfile.close()
