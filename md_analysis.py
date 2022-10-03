# a script to analyse the output of an MD run in the form of a .xyz(e) trajectory
# basically a wrapper to some ASE functionalities

# written by Pascal Salzbrenner, pts28@cam.ac.uk

import os
import shutil
import ase, ase.io, ase.md.analysis

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

     input_file = open(input_filename, "r")
     reduced_file = open(reduced_filename, "w")

     input_file.close()
     reduced_file.close()

task = input("What task would you like to carry out? [diffusion, rdf] ").lower()
