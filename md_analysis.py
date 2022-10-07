# a script to analyse the output of an MD run in the form of a .xyz(e) trajectory
# basically a wrapper to some ASE functionalities

# written by Pascal Salzbrenner, pts28@cam.ac.uk

import os
import shutil
import ase, ase.io, ase.md.analysis

from exceptions import InputError

# user_defined input
input_filename = input("What is the name of the trajectory file? ")

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

    reduced_filename = "{}_atoms_{}_{}.xyz".format(fileroot, first_atom, last_atom)

    # calculate the new number of atoms
    num_atoms_output = last_atom - first_atom + 1 # eg if first_atom = 1, and last_atom = 2, I want num_atoms_output to be 2 = 2 - 1 + 1

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

        # read past the line giving the geometry - it will now be wrong and ASE can work out the actual geometry from the atom positions
        input_file.readline()

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

    filename = reduced_filename
else:
    filename = input_filename

task = input("What task would you like to carry out? [diffusion, rdf] ").lower()

if task.startswith("d"):
    # diffusion coefficient
    time_step = float(input("What was the time step of the calculation? [fs] "))
    ignore_images = int(input("At what image would you like to start the calculation? "))
    num_segments = int(input("Into how many segments (for the purpose of averaging) would you like to split the data? "))

    trajectory = ase.io.read(filename, ":")

    diffusion_coefficient = ase.md.analysis.DiffusionCoefficient(trajectory, timestep=time_step*ase.units.fs)

    diffusion_coefficient.calculate(ignore_n_images=ignore_images,number_of_segments=num_segments)[0]

elif task.startswith("r"):
    # rdf

else:
    # not a valid option
    raise InputError("Task", "You have specified a non-implemented task. Currently supported: [d]iffusion, [r]df.")
