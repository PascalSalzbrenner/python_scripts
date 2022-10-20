# a script to analyse the output of an MD run in the form of a .xyz(e) trajectory
# basically a wrapper to some ASE functionalities

# written by Pascal Salzbrenner, pts28@cam.ac.uk

import os
import shutil
import ase, ase.io, ase.geometry.analysis, ase.md.analysis

import numpy as np
import matplotlib.pyplot as plt

from exceptions import InputError

######################################## define function to build a new xyz file containing only certain atoms from the original file ##################################################

def build_reduced_xyz_file(fileroot, first_atom, last_atom):
    """In the case of coexistence simulations (or maybe also other situations), it is desirable to only include certain atoms in a calculation.
    This function creates a new xyz file, where the unit cell is unchanged, but the only atoms included are those specified

    :param str fileroot: the root of the <fileroot>.xyz file
    :param int first_atom: the first atom to be included
    :param int final_atom: the final atom to be included

    :returns str reduced_filename: the name of the new .xyz file"""

    input_filename = "{}.xyz".format(fileroot)
    reduced_filename = "{}_atoms_{}_{}.xyz".format(fileroot, first_atom, last_atom)

    # calculate the new number of atoms
    num_atoms_output = last_atom - first_atom + 1 # eg if first_atom = 1, and last_atom = 2, we want num_atoms_output to be 2 = 2 - 1 + 1

    input_file = open(input_filename, "r")
    reduced_file = open(reduced_filename, "w")

    # read number of atoms
    num_atoms_input = int(input_file.readline().split()[0])

    # rewind file
    input_file.seek(0)

    # iterate over file
    for line in input_file:
        # first line of a new frame - number of atoms - replace this with the new number
        reduced_file.write("{}\n".format(num_atoms_output))

        # directly write the line relating to the lattice to the new file - it will now be wrong, but for the calculation of the diffusion coefficient at least, this should not matter
        reduced_file.write(input_file.readline())

        # reset counter
        atom_counter = 0

        while atom_counter < num_atoms_input:
            positions_line = input_file.readline()

            if atom_counter >= first_atom-1 and atom_counter <= last_atom-1:
                # the range we want to keep
                reduced_file.write(positions_line)

            atom_counter += 1

    input_file.close()
    reduced_file.close()

    return reduced_filename

# user_defined input
input_filename = input("What is the name of the trajectory file? ")

# check if the input file exists
ls = os.listdir()

if not input_filename in ls:
    raise InputError("Unable to read input file", "There is no file named {} in the current directory.".format(input_filename))


fileroot = input_filename.rstrip("e").replace(".xyz", "")

# xyz and xyze both work, but the xyze file suffix confuses ASE, so in this case we copy the file to one with the acceptable .xyz suffix
if input_filename.endswith(".xyze"):
    shutil.copyfile("{}.xyze".format(fileroot), "{}.xyz".format(fileroot))

# the calculation will need some time to converge - we ignore this in the calculation
ignore_images = int(input("At which image in the trajectory would you like to start the calculation? "))-1

coexistence = input("Is the run a coexistence simulation? [y/[n]] ").lower()

if coexistence.startswith("y"):
    coexistence = True
else:
    coexistence = False

# in a coexistence simulation, we usually want to do these calculations separately for the solid and liquid phases
# here, we let the user specify a subset of all atoms over which the calculation is carried out - numbered from 1 to natoms in the .xyz file
if coexistence:
    first_atom = int(input("What is the first atom you want to include in the calculation? "))
    last_atom = int(input("What is the last atom you want to include in the calculation? "))

task = input("What task would you like to carry out? [diffusion, rdf] ").lower()

if task.startswith("d"):
    # diffusion coefficient

    if coexistence:
        filename = build_reduced_xyz_file(fileroot, first_atom, last_atom)
    else:
        filename = input_filename

    time_step = float(input("What was the time step between images in {}? [fs] ".format(filename)))
    num_segments = int(input("Into how many segments (for the purpose of averaging) would you like to split the data? "))

    trajectory = ase.io.read(filename, ":")

    diffusion_coefficient_object = ase.md.analysis.DiffusionCoefficient(trajectory, timestep=time_step*ase.units.fs)

    diffusion_coefficient_object.calculate(ignore_n_images=ignore_images,number_of_segments=num_segments)

    # get number of elements
    num_elements = len(diffusion_coefficient_object.types_of_atoms)

    # the diffusion coefficients are output for each atom, "in alphabetical order"
    # this order is the same as that returned by diffusion_coefficient_object.types_of_atoms
    diffusion_coefficients, standard_deviations = diffusion_coefficient_object.get_diffusion_coefficients()

    outfile = open("diffusion_coefficient.txt", "w")
    outfile.write("# Element; Diffusion Coefficient [A**2/fs]; Diffusion Coefficient [cm**2/s]; Standard Deviation [A**2/fs]; Standard Deviation [cm**2/s]\n")

    for i in range(num_elements):
        diffusion_coefficient = diffusion_coefficients[i]*ase.units.fs
        standard_deviation = standard_deviations[i]*ase.units.fs

        outfile.write("{: <5} {: <10.6f} {: <10.6f} {: <10.6f} {: <10.6f}\n".format(diffusion_coefficient_object.types_of_atoms[i], diffusion_coefficient, diffusion_coefficient/10, standard_deviation, standard_deviation/10))

    outfile.close()

    # plot using ASE built-in plotting utility
    axes = plt.gca()
    diffusion_coefficient_object.plot(ax=axes)
    plt.savefig("diffusion_coefficient.pdf")

    # the data for the average diffusion coefficient can be read out from the axes.lines objects
    # specifically, there are num_segments * (num_elements + 1) entries, where for each segment the last entry corresponds to a straight line separating different segments

    for i in range(num_elements):
        # open a file for each element, where the x-y data for all segments are written, with different segments separated by a blank line

        element_name = diffusion_coefficient_object.types_of_atoms[i]

        outfile = open("{}_MSD_data.dat".format(element_name), "w")
        outfile.write("# Time [fs]; MSD [A**2]\n")


        for j in range(num_segments):
            x_data, y_data = axes.lines[i+j*(num_elements+1)].get_data()

            for k in range(len(x_data)):
                outfile.write("{: <10.6f} {: <10.6f}\n".format(x_data[k], y_data[k]))

            outfile.write("\n")

        outfile.close()

elif task.startswith("r"):
    # rdf

    # the calculated rdf is an average, but we also print out snapshots in case there is a phase transition

    # read trajectory and generate analysis class

    trajectory = ase.io.read(input_filename, index="{}:".format(ignore_images))
    geometry = ase.geometry.analysis.Analysis(trajectory)

    rmax = float(input("What is the maximum distance at which you want to compute the RDF? "))
    nbins = int(input("Into how many bins would you like to divide the RDF? "))

    # generate x range
    x_step = rmax / nbins
    x_min = x_step / 2

    x_values = []

    for i in range(nbins):
        x_values.append(x_min + x_step*i)

    final_image = input("Up to which image would you like to calculate the RDF? (leave blank for all images) ")

    if not final_image:
        final_image = geometry.nImages
    else:
        final_image = int(final_image)

    # because we are only reading from ignore_images, the image indices are renumbered. We subtract that number from the user-supplied final_image, which will reference the full trajectory
    final_image -= ignore_images

    image_slice = slice(0, final_image)

    snapshot_step = int(input("At which frequency would you like to print snapshots of the RDF? "))

    # it is possible to specify "elements" which are included in the RDF - these can be either a range of atoms, or a set of actual chemical elements, but not a combination of both
    element_list = []

    if not coexistence:
        # if it is a coexistence calculation, the user will have specified a range of atoms, so they can't also specify a set of chemical elements
        element = input("We will now compile a list of elements you want to include in the RDF calculation. Once you have supplied all desired elements, press 'enter' without supplying another element. If you want to include all elements, press 'enter' now: ")

        while element:
            element_list.append(element)
            element = input("Please enter the next element: ")

        if not element_list:
            # the element_list is empty, set it to all elements
            element_list = list(set(trajectory[0].get_chemical_symbols()))

    else:
        # in this case, we create a list including only the atoms in the range specified by first_atom, last_atom

        for i in range(first_atom-1, last_atom):
            element_list.append(i)

    # actually calculate the RDF - calculates one RDF per included image
    rdf = geometry.get_rdf(rmax,nbins,imageIdx = image_slice, elements=element_list)

    average_rdf = np.mean(rdf,axis=0)

    # write output
    datafile = open("average_rdf_{}_{}_images.dat".format(ignore_images, final_image + ignore_images), "w")
    
    datafile.write("# bin centre [A]; RDF\n")

    for i in range(len(average_rdf)):
        datafile.write("{} {}\n".format(x_values[i], average_rdf[i]))

    datafile.close()

    with open("average_rdf_{}_{}_images.gnu".format(ignore_images+1, final_image + ignore_images+1), "w") as plotfile:

        plotfile.write("set terminal postscript eps colour font 'Helvectica,20'\n")
        plotfile.write("set output '| epstopdf --filter --outfile=average_rdf_{}_{}_images.pdf'\n".format(ignore_images+1, final_image + ignore_images+1))
        plotfile.write("set mxtics 2\n")
        plotfile.write("set mytics 2\n")
        plotfile.write("set boxwidth {}\n".format(x_step))
        plotfile.write("set style fill solid\n")
        plotfile.write("set xlabel '|r| [A]'\n")
        plotfile.write("set ylabel 'g(r)'\n")
        plotfile.write("set xrange [0:{}]\n".format(rmax))
        plotfile.write("set yrange [0:]\n")
        plotfile.write("plot 'average_rdf_{}_{}_images.dat' u 1:2 w boxes lc rgb '#DC143C' notitle".format(ignore_images+1, final_image + ignore_images+1))

    # use matplotlib to generate some quick plots of the snapshots - not super nice-looking as they are intended for checking purposes only

    for i in range(0, len(rdf), snapshot_step):
        plt.plot(x_values, rdf[i])
        plt.xlabel("|r| [A]")
        plt.ylabel("g(r)")
        plt.savefig("rdf_snapshot_{}.pdf".format(ignore_images+i+1))
        plt.clf()
    
else:
    # not a valid option
    raise InputError("Task", "You have specified a non-implemented task. Currently supported: [d]iffusion, [r]df.")
