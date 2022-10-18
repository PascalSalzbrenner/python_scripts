# a script to analyse the output of an MD run in the form of a .xyz(e) trajectory
# basically a wrapper to some ASE functionalities

# written by Pascal Salzbrenner, pts28@cam.ac.uk

import os
import shutil
import ase, ase.io, ase.md.analysis

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

    # the diffusion coefficients are output for each atom, "in alphabetical order"
    # this order is the same as that returned by diffusion_coefficient_object.types_of_atoms
    diffusion_coefficients, standard_deviations = diffusion_coefficient_object.get_diffusion_coefficients()

    outfile = open("diffusion_coefficient.txt", "w")
    outfile.write("# Element; Diffusion Coefficient [A**2/fs]; Diffusion Coefficient [cm**2/s]; Standard Deviation [A**2/fs]; Standard Deviation [cm**2/s]\n")

    for i in range(len(diffusion_coefficient_object.types_of_atoms)):
        diffusion_coefficient = diffusion_coefficients[i]*ase.units.fs
        standard_deviation = standard_deviations[i]*ase.units.fs

        outfile.write("{} {} {} {} {}\n".format(diffusion_coefficient_object.types_of_atoms[i], diffusion_coefficient, diffusion_coefficient/10, standard_deviation, standard_deviation/10))

    outfile.close()



elif task.startswith("r"):
    # rdf
    print("To be implemented soon!")
else:
    # not a valid option
    raise InputError("Task", "You have specified a non-implemented task. Currently supported: [d]iffusion, [r]df.")
