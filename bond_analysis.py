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
# bonds: simply prints out data about the bonds - name, length, (if read from a .castep file) populations - in the bonds.dat file

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
        if "External pressure/stress" in line:
            # this block is only present once; the first element of the next line gives the pressure
            pressure = float(bonds_file.readline().split()[0])
        elif "Bond" in line:
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

    return (bonds, pressure)

# get necessary input
task = input("Which task would you like to run? (PDF, RDF, bond_length, bond_population, bonds - see the code header for descriptions) ").lower()

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

# write Gnuplot file
plotfile = open("{}.gnu".format(task), "w")
plotfile.write("set terminal postscript eps colour\n")
plotfile.write("set output '| epstopdf --filter --outfile={}.pdf'\n".format(task))
plotfile.write("set key top right\n")
plotfile.write("set key box lt -1 lw 2 width 2 height 1.5 opaque\n")
# these parts will be common to all tasks, whereas others epend in the task and are written in the corresponding if-block

if task == "pdf":

    if input_file_type=="castep":
        bonds, pressure = get_bonds_from_castep(input_file)
    else:
        structure = Structure(input_file)
        structure.get_bonds()
        bonds = structure.bonds

elif task == "rdf":

    structure = Structure(input_file)

    # construct supercell from the structure, loop over atoms in original cell

elif task == "bond_length" or task == "bond_population":

    indices = "What bonds (indexed from the shortest=1) do you want to take into account?\nPass a whitespace-separated list of an arbitrary number.\na-b will lead to an average of all bonds in the range [a,b]. ".split()

    # for these two tasks, loop over all files of the right type in the directory

    # determine all those files
    analysis_files = []
    for item in os.listdir():
        if "input_file_type" in item:
            analysis_files.append(item)

    # determine which element of the data we need to read (1 for length, 2 for population)
    if task == "bond_length":
        data_index = 1
    else:
        data_index = 2

    # data container
    data_list=[]

    for file in analysis_files:
        if input_file_type == "castep":
            bonds, pressure = get_bonds_from_castep(file)
            pressure_units = "GPa"
            length_units = " [A]"

        else:
            structure = Structure(file)
            structure.get_bonds()
            bonds = structure.bonds
            pressure = structure.pressure
            pressure_units = structure.pressure_units
            length_units = structure.length_units

        # initialise string for the line
        write_str = "{}".format(pressure)

        # iterate over indices
        for index in indices:

            # check if we are averaging over a number of bonds
            if "-" in index:
                start_end_indices = index.split("-")
                start_index = int(start_end_indices[0])
                end_index = int(start_end_indices[1])
                num_points = end_index - start_index + 1
                average = 0

                for i in range(start_index-1, end_index):
                    average += bonds[i][data_index]

                write_str += " {}".format(average/num_points)
            else:
                # only a single point
                write_str += " {}".format(bonds[index-1][data_index])

    # determine the y-units depending on the task
    if task == "bond_length":
        y_units = " [{}]".format()
    else:
        y_units = ""

    # sort data in ascending order of pressure
    data_list.sort()

    # write data
    datafile = open("{}.dat".format(task), "w")
    datafile.write("# Pressure [{}]; {}{}\n".format(pressure_units, task, y_units))
    datafile.write("\n".join(data_list))
    datafile.close()


    # plotting stuff
    # redefine gnuplot linetypes with nice colours
    plotfile.write("set linetype 1 lc rgb '#DC143C'\n")
    plotfile.write("set linetype 2 lc rgb '#DC143C'\n")
    plotfile.write("set linetype 3 lc rgb '#D95F02'\n")
    plotfile.write("set linetype 4 lc rgb '#D95F02'\n")
    plotfile.write("set linetype 5 lc rgb '#E6AB02'\n")
    plotfile.write("set linetype 6 lc rgb '#E6AB02'\n")
    plotfile.write("set linetype 7 lc rgb '#66A61E'\n")
    plotfile.write("set linetype 8 lc rgb '#66A61E'\n")
    plotfile.write("set linetype 9 lc rgb '#8000C4'\n")
    plotfile.write("set linetype 10 lc rgb '#8000C4'\n")
    plotfile.write("set linetype 11 lc rgb '#7570B3'\n")
    plotfile.write("set linetype 12 lc rgb '#7570B3'\n")
    plotfile.write("set linetype 13 lc rgb '#E7298A'\n")
    plotfile.write("set linetype 14 lc rgb '#E7298A'\n")
    plotfile.write("set linetype 15 lc rgb '#1E90FF'\n")
    plotfile.write("set linetype 16 lc rgb '#1E90FF'\n")
    plotfile.write("set linetype 17 lc rgb '#1B9E77'\n")
    plotfile.write("set linetype 18 lc rgb '#1B9E77'\n")
    plotfile.write("set linetype 19 lc rgb '#B8860B'\n")
    plotfile.write("set linetype 20 lc rgb '#B8860B'\n")
    plotfile.write("set linetype 21 lc rgb '#20C2C2'\n")
    plotfile.write("set linetype 22 lc rgb '#20C2C2'\n")
    plotfile.write("set linetype cycle 22\n")

    plotfile.write("set xlabel '{} [{}]'\n".format(pressure, pressure_units))

    # format the task string so as to use it for the y-axis label
    task_strings = task.split("_")
    task_string_nice = "{} {}".format(task_strings[0].capitalize(), task_strings[1].capitalize())
    plotfile.write("set ylabel '{}{}'\n".format(task_string_nice, y_units))

    # compile string for plotting the data
    plot_string = "plot"

    for i in range(len(indices)):
        plot_string.append(" {0}.dat u 1:{1} w points pt 7 ps 1.5 title 'Bond(s) {2}', {0}.dat u 1:{1} smooth mcsplines notitle".format(task, i+2, indices[i]))

    plotfile.write(plot_string)


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

plotfile.close()
