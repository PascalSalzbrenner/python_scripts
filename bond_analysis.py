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
from utilities import lattice_basis_to_cartesian, get_structure_from_castep, get_structure_from_cell
from exceptions import InputError

# define two functions to determine the bonds, either from a structure, or from the CASTEP file, where they are already present and must
# only be read
def get_bonds_from_castep

def get_bonds_from_structure
    # structure must already have been read and will be passed as input argument

# get necessary input
task = input("Which task would you like to run? (PDF, RDF, bond_length, bond_population - see the code header for descriptions) ").lower()

if task == "rdf" or task == "pdf" or task == "bond_length":
    # different input files are possible
    input_file = input("Which input file type should be used? (currently supported: CASTEP, cell) ").lower()
elif task == "bond_population":
    # only possile from a .castep file
    input_file = "castep"
else:
    # a non-implemented task has been requested
    raise InputError("task", "You have requested a task which is not implemented. Task must be one of the following: PDF, RDF, bond_length, bond_population.")

if input_file=="cell":
    get_structure_from_cell

if task == "pdf":

    if input_file=="castep":
        get_bonds_from_castep
    else:
        # any structure format will already have been parsed into the universal internal format
        get_bonds_from_structure

elif task == "rdf":

    if input_file=="castep":
        get_structure_from_castep

elif task == "bond_length":

    if input_file=="castep":
        get_bonds_from_castep
    else:
        get_bonds_from_structure

elif task == "bond_population":

    get_bonds_from_castep

# in the relevant if blocks, do the analysis
