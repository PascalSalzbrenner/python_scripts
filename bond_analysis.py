# script to do various kinds of analysis and plots of the bonds in a crystal
# currently works using the bond analysis in a .castep output file of a CASTEP run
# I could write parsers for the structure input files of different DFT codes and then do the length / RDF type analyses for them

# Features:
# PDF: bins all bonds, divides by the number of atoms in the cell, and plots a histogram of frequency vs length
#      if a structure is given, essentially repeats the analysis of the .castep file
#      that is, PDF analysis is confined to unit cell
# RDF: proper RDF analysis - constructs 9x9x9 supercell from the structure (which can be read from the .castep or structure file)
#      and calculates proper RDF as defined on slide 11.10 of the atomistic modelling lecture notes (note that file is actually Lecture 10)
#      by iterating over all atoms in the central cell and averaging their individual distribution functions
#      (see equation 8 in https://scripts.iucr.org/cgi-bin/paper?S0021889800019993)
# bond length analysis: Takes the indices (bonds ranked by length from shortest to longest) of an arbitrary number of bonds and tracks their
#                       length(s) as a function of pressure
# bond population analysis: Likewise, but for the population. Only works from a .castep file
# bonds: simply prints out data about the bonds - name, length, (if read from a .castep file) populations - in the bonds.dat file

# if bond_length or bond_population is selected, the code reads all relevant files in the directory, with the pressure determined from them

# written by Pascal Salzbrenner, pts28

import os
import re
import numpy as np
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
        if "External pressure/stress" in line:
            # this block is only present once; the first element of the next line gives the pressure
            pressure = float(bonds_file.readline().split()[0])
        elif "Total number of ions in cell" in line:
            # last element of this line gives the number of atoms
            num_atoms = float(line.split()[-1])
        elif "Bond" in line:
            # read past the next line; each line after that until a line full of "=" contains a bond
            bonds_file.readline()

            for bonds_line in bonds_file:
                if "===" in bonds_line:
                    break
                else:
                    bonds_data = bonds_line.split()
                    bonds.append(("{}{} - {}{}".format(bonds_data[0], bonds_data[1], bonds_data[3], bonds_data[4]), float(bonds_data[-1]),
                                  float(bonds_data[-2])))

    bonds_file.close()

    # sort by bond length
    bonds.sort(key=itemgetter(1))

    return bonds, pressure, num_atoms

# get necessary input
task = input("Which task would you like to run? (PDF, RDF, bond_length, bond_population, bonds - see the code header for descriptions) ").lower().strip()

if task == "rdf" or task == "pdf" or task == "bonds":
    # different input files are possible
    input_file = input("What is the name of the input file? (currently supported: .castep, .cell files) ").strip()
    input_file_type = input_file.split(".")[-1]
    input_file_root = "_" + ".".join(input_file.split(".")[:-1])
elif task == "bond_length":
    input_file_type = input("What is the type of input file? (currently supported: castep, cell) ").lower().strip()
    input_file_root = ""
elif task == "bond_population":
    # only possile from a .castep file
    input_file_type = "castep"
    input_file_root = ""
else:
    # a non-implemented task has been requested
    raise InputError("task",
    "You have requested a task which is not implemented. Task must be one of the following: PDF, RDF, bond_length, bond_population, bonds.")

if task != "bonds":
    # write Gnuplot file
    # these bits will be common to all tasks (except bonds), whereas others depend on the task and are written in the corresponding if-block
    plotfile = open("{}{}.gnu".format(task, input_file_root), "w")
    plotfile.write("set terminal postscript eps colour\n")
    plotfile.write("set output '| epstopdf --filter --outfile={}{}.pdf'\n".format(task, input_file_root))

if task == "pdf":

    # read the width of the bins in which the PDF data will be grouped
    bin_width = float(input("What bin width (length units should be those of your input file) would you like to use for the PDF? "))
    half_bin_width = bin_width/2
    num_decimals = len(str(half_bin_width)) - 2

    if input_file_type=="castep":
        bonds, pressure, num_atoms = get_bonds_from_castep(input_file)
        length_units = "[A]"
    else:
        structure = Structure(input_file)
        structure.get_bonds()
        bonds = structure.bonds
        length_units = structure.length_units
        num_atoms = structure.num_atoms

    # determine the length of the longest bond in the system
    max_length = bonds[-1][1]

    # set up a list for the lower bounds of the different bins
    bin_value = 0
    bin_list = []

    # first bin value
    bin_list.append(bin_value)

    # generate bin values
    while bin_value <= max_length:

        bin_value += bin_width

        bin_list.append(bin_value)

    # set up an array to contain the occupations of the different bins
    bin_occupations = np.zeros(len(bin_list), dtype=int)

    # the final value in bin_list will be the upper bound of the final bin (no length can be larger than it by definition)

    # sort bond lengths into bins
    for bond in bonds:
        bond_length = bond[1]

        for i in range(len(bin_list)-1):
            if bond_length >= bin_list[i] and bond_length < bin_list[i+1]:
                bin_occupations[i] += 1
                break

    bin_occupations = bin_occupations / num_atoms

    # write output data
    datafile = open("pdf{}.dat".format(input_file_root), "w")
    datafile.write("# bin middle {}; number of bonds in bin\n".format(length_units))
    for i in range(len(bin_list)-1):
        datafile.write("{0: <.{2}f} {1: >.5f}\n".format(bin_list[i]+half_bin_width, bin_occupations[i], num_decimals))
    datafile.close()

    # write plotting commands
    plotfile.write("set boxwidth {}\n".format(bin_width))
    plotfile.write("set style fill solid\n")
    plotfile.write("set xlabel 'Bond length {}'\n".format(length_units))
    plotfile.write("set ylabel 'Number per atom'\n")
    plotfile.write("set xrange [0:]\n")
    plotfile.write("set mxtics 2\n")
    plotfile.write("set yrange [0:]\n")
    plotfile.write("plot 'pdf{}.dat' u 1:2 w boxes lc rgb '#DC143C' notitle".format(input_file_root))

elif task == "rdf":

    # read the width of the radial shell
    shell_width = float(input("How wide (length units should be those of your input file) should the shell used to evaluate the RDF be? "))
    max_radius = float(input("What is the maximum radius (length units of the input file) at which to evalue the RDF? "))
    half_shell_width = shell_width/2
    num_decimals = len(str(half_shell_width)) - 2

    structure = Structure(input_file)

    # construct supercell from the structure
    supercell_positions = []

    for i in range(-1,2):
        for j in range(-1,2):
            for k in range(-1,2):
                for position in structure.positions_abs:
                    supercell_positions.append(position+np.dot(structure.lattice.T, np.array([i, j, k])))

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
    # normalise by the shell volumes and the density (= N/V, hence multiply by V/N)
    # note that as the supercell is 27 replicas of the primitive cell, the density is of course the same
    normalised_shell_occupations = (structure.volume/structure.num_atoms)*average_shell_occupations/shell_volumes

    # write output data
    datafile = open("rdf{}.dat".format(input_file_root), "w")
    datafile.write("# shell middle r [{}]; g(r)\n".format(structure.length_units))
    for i in range(len(shell_list)-1):
        datafile.write("{0: <.{2}f} {1: >.5f}\n".format(shell_list[i]+half_shell_width, normalised_shell_occupations[i], num_decimals))
    datafile.close()

    # write plotting commands
    plotfile.write("set boxwidth {}\n".format(shell_width))
    plotfile.write("set style fill solid\n")
    plotfile.write("set xlabel 'r [{}]'\n".format(structure.length_units))
    plotfile.write("set ylabel 'g(r)'\n")
    plotfile.write("set xrange [0:{}]\n".format(max_radius))
    plotfile.write("set mxtics 2\n")
    plotfile.write("set yrange [0:]\n")
    plotfile.write("plot 'rdf{}.dat' u 1:2 w boxes lc rgb '#DC143C' notitle, '' u 1:2 smooth mcsplines lw 2 lc rgb '#000080' notitle".format(input_file_root))

elif task == "bond_length" or task == "bond_population":

    indices = input("What bonds (indexed from the shortest=1, or named by <initial_atom><initial_atom_number>-<final_atom><final_atom_number>) do you want to take into account?\nPass a whitespace-separated list of an arbitrary number of bonds.\nIf using indices, a-b will lead to an average of all bonds in the range [a,b]. ").split()

    # check if the bonds are named or indexed
    if re.match(r"^[a-zA-Z]", indices[0]):
        # the first element of the string is a letter, hence we are dealing with a named bond
        named_bonds = True
    else:
        named_bonds = False

    # for these two tasks, loop over all files of the right type in the directory
    # determine all those files
    analysis_files = []
    for item in os.listdir():
        if input_file_type in item:
            analysis_files.append(item)

    if not analysis_files:
        # not matching files found
        raise InputError("Finding input files for {}".format(task),
        "You specified {} as the input file type for {}, but no files if this type were found.".format(format(input_file_type, task)))

    # determine which element of the data we need to read (1 for length, 2 for population)
    if task == "bond_length":
        data_index = 1
    else:
        data_index = 2

    # data container
    data_list=[]

    for file in analysis_files:
        if input_file_type == "castep":
            bonds, pressure, num_atoms = get_bonds_from_castep(file)
            pressure_units = "GPa"
            length_units = "A"

            if named_bonds:
                # set up dictionary of bonds for easy access of named bonds
                bonds_dict = {}
                for bond in bonds:
                    bonds_dict[bond[0]] = (bond[1], bond[2])
        else:
            structure = Structure(file)
            structure.get_bonds()
            bonds = structure.bonds
            bonds_dict = structure.bonds_dict
            pressure = structure.pressure
            pressure_units = structure.pressure_units
            length_units = structure.length_units

        # initialise string for the line
        write_str = "{: <10d}".format(int(pressure))

        # iterate over indices
        for index in indices:

            if named_bonds:
                # bonds are given by names
                split_name = index.split("-")
                write_str += " {: <10.5f}".format(bonds_dict["{} - {}".format(split_name[0], split_name[1])][data_index-1])
                # note about the previous line that the input bond names are not whitespace separated, but those in the dict are
                # also, the tuple misses the first element (the bond name), which is the library key
            else:
                # bonds are given by indices
                # check if we are averaging over a number of bonds
                if "-" in index:
                    start_end_indices = index.split("-")
                    start_index = int(start_end_indices[0])
                    end_index = int(start_end_indices[1])
                    num_points = end_index - start_index + 1
                    average = 0

                    for i in range(start_index-1, end_index):
                        average += bonds[i][data_index]

                    write_str += " {: <10.5f}".format(average/num_points)
                else:
                    # only a single point
                    write_str += " {: <10.5f}".format(bonds[index-1][data_index])

        data_list.append(write_str)

    # determine the y-units depending on the task
    if task == "bond_length":
        y_units = " [{}]".format(length_units)
    else:
        y_units = ""

    # sort data in ascending order of pressure
    data_list.sort()

    # format the task string so as to use it for the output
    task_strings = task.split("_")
    task_string_nice = "{} {}".format(task_strings[0].capitalize(), task_strings[1].capitalize())

    # write data
    datafile = open("{}.dat".format(task), "w")
    title_str = "# Pressure [{}];".format(pressure_units)
    for index in indices:
        title_str += " {} of bond(s) {}{};".format(task_strings[-1], index, y_units)
    datafile.write("{}\n".format(title_str.rstrip(";")))
    datafile.write("\n".join(data_list))
    datafile.close()

    # plotting stuff
    plotfile.write("set key top right\n")
    plotfile.write("set key box lt -1 lw 2 width 2 height 1.5 opaque\n")

    # redefine gnuplot linetypes with nice colours
    plotfile.write("set linetype 1 lc rgb '#3633F5'\n") # same dark blue as used for molecular bonds in VESTA
    plotfile.write("set linetype 2 lc rgb '#3633F5'\n")
    plotfile.write("set linetype 3 lc rgb '#FF0000'\n") # classic red - used for longer molecular bonds in VESTA
    plotfile.write("set linetype 4 lc rgb '#FF0000'\n")
    plotfile.write("set linetype 5 lc rgb '#16B0FF'\n") # same light blue as used for intermolecular bonds in VESTA
    plotfile.write("set linetype 6 lc rgb '#16B0FF'\n")
    plotfile.write("set linetype 7 lc rgb '#D95F02'\n") # orange
    plotfile.write("set linetype 8 lc rgb '#D95F02'\n")
    plotfile.write("set linetype 9 lc rgb '#E6AB02'\n") # yellow
    plotfile.write("set linetype 10 lc rgb '#E6AB02'\n")
    plotfile.write("set linetype 11 lc rgb '#66A61E'\n") # dark green
    plotfile.write("set linetype 12 lc rgb '#66A61E'\n")
    plotfile.write("set linetype 13 lc rgb '#8000C4'\n") # purple
    plotfile.write("set linetype 14 lc rgb '#8000C4'\n")
    plotfile.write("set linetype 15 lc rgb '#7570B3'\n") # silver-ish medium blue
    plotfile.write("set linetype 16 lc rgb '#7570B3'\n")
    plotfile.write("set linetype 17 lc rgb '#E7298A'\n") # pink
    plotfile.write("set linetype 18 lc rgb '#E7298A'\n")
    plotfile.write("set linetype 19 lc rgb '#1B9E77'\n") # medium green
    plotfile.write("set linetype 20 lc rgb '#1B9E77'\n")
    plotfile.write("set linetype 21 lc rgb '#B8860B'\n") # gold
    plotfile.write("set linetype 22 lc rgb '#B8860B'\n")
    plotfile.write("set linetype 23 lc rgb '#20C2C2'\n") # turqouise
    plotfile.write("set linetype 24 lc rgb '#20C2C2'\n")
    plotfile.write("set linetype cycle 24\n")

    plotfile.write("set xlabel 'Pressure [{}]'\n".format(pressure_units))

    plotfile.write("set ylabel '{}{}'\n".format(task_string_nice, y_units))

    # compile string for plotting the data
    plot_string = "plot"

    for i in range(len(indices)):
        plot_string += " '{0}.dat' u 1:{1} w points pt 7 ps 1.5 title 'Bond(s) {2}', '{0}.dat' u 1:{1} smooth bezier lw 2.5 notitle,".format(task, i+2, indices[i])

    plotfile.write(plot_string.rstrip(","))

elif task == "bonds":

    if input_file_type == "castep":
        bonds, pressure, num_atoms = get_bonds_from_castep(input_file)
        length_units = "Angstrom"
    else:
        structure = Structure(input_file)
        structure.get_bonds()
        bonds = structure.bonds
        length_units = structure.length_units

    outfile = open("bonds{}.dat".format(input_file_root), "w")

    first_line_str = "# Bond; Length [{}]; ".format(length_units)

    if input_file_type == "castep":
        # write out bond population as well
        first_line_str += "Population; bond number\n"
    else:
        # only write the bond number
        first_line_str += "bond number\n"

    outfile.write(first_line_str)

    for i in range(len(bonds)):

        write_str = "{:<15} {:<10.5f} ".format(bonds[i][0], bonds[i][1])

        if input_file_type == "castep":
            # write out bond population as well
            write_str += "{:< 10.2f} {:>5}\n".format(bonds[i][2], i+1)
        else:
            # only write out bond number
            write_str += "{:10}\n".format(i+1)

        outfile.write(write_str)

    outfile.close()

if task != "bonds":
    plotfile.close()
