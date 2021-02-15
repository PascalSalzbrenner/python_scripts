# short script to fit to the Birch-Murnaghan equation of state (DOI: 10.1103/PhysRev.71.809)
# input: energy.txt file as produced by my geometry_optimisation_cubic routine

# determines whether the structure is FCC, BCC, or simple cubic, and implements the right prefactor in the determination of the lattice
# parameter from the volume
# so that the lattice parameter a enters as follows into the descriptions of the different lattices:

# simple cubic:
# a 0 0
# 0 a 0
# 0 0 a

# bcc
# -a a a
# a -a a
# a a -a

# fcc
# 0 a a
# a 0 a
# a a 0

# and of course any permutation of the lattice vectors

# written by Pascal Salzbrenner - pts28@cam.ac.uk

import sys
import math
import numpy as np
import spglib as sl
from scipy.optimize import curve_fit

# define Birch-Murnaghan equation of state for E(V) - V is input, the other variables are
# parameters to be fit
def birch_murnaghan_eos(V, E_0, V_0, B_0, B_prime):
    """The Birch-Murnaghan equation of state as defined in DOI: 10.1103/PhysRev.71.809"""

    return E_0 + 9.0*V_0*B_0/(16.0) * ((((V_0/V)**(2))**(1.0/3.0)-1)**(3)*B_prime
           + (((V_0/V)**(2))**(1.0/3.0) - 1)**(2)*(6.0-4.0*((V_0/V)**(2))**(1.0/3.0)))

dft_code = sys.argv[1].lower()

# set up lattice, positions, and numbers lists needed by spglib
lattice = []
positions = []
numbers = []

# auxiliary list of atoms
atoms_list = []

if dft_code == "castep" or dft_code == "profess":
    fileroot = sys.argv[2]

    if dft_code == "castep":
        structure_file = open("{}.cell_for_fit".format(fileroot), "r")
    else:
        structure_file = open("{}.ion_for_fit".format(fileroot), "r")

    # structure parser - I need lattice, positions, numbers, as per https://spglib.github.io/spglib/python-spglib.html
    for line in structure_file:

        if not line.split():
            continue # the line is empty

        if line.lower().startswith("%b") or line.lower().startswith("% b"):
            # the beginning of an input block

            if "lattice" in line.lower():

                for lattice_vector in structure_file:

                    if lattice_vector.lower().startswith("%e") or lattice_vector.lower().startswith("% e"):
                        break # the end of the input block
                    elif lattice_vector.split():
                        lattice.append([float(el) for el in lattice_vector.split()[:3]])

            elif "position" in line.lower():

                for position in structure_file:

                    if position.lower().startswith("%e") or position.lower().startswith("% e"):
                        break # the end of the input block
                    elif position.split():
                        position_line_split = position.split()
                        # first element of the list is the name of the element
                        if position_line_split[0] not in atoms_list:
                            atoms_list.append(position_line_split[0])

                        numbers.append(atoms_list.index(position_line_split[0]))

                        positions.append([float(el) for el in position_line_split[1:4]])

            else:
                continue

        else:
            continue


elif dft_code == "vasp":
    structure_file = open("POSCAR_for_fit", "r")

    # structure parser

    structure_file.readline()

    # the next line contains a scaling factor applied to all lattice vectors
    scaling_factor = float(structure_file.readline().split()[0])

    # the next three lines contain the lattice_vectors
    for i in range(3):
        lattice_vector = structure_file.readline()
        lattice.append([float(el)*scaling_factor for el in lattice_vector.split()[:3]])

    # next line contains the atom names; read past it as they are not necessary
    structure_file.readline()

    # next line contains the how often each atom occurs in sequence
    atoms_list = [int(el) for el in structure_file.readline().split()] # this will fail if there's a comment in this line

    #  next line describes whether the coordinates are direct or Cartesian
    # direct is assumed as the entire geometry optimisation would in any case fail if it weren't
    structure_file.readline()

    # initialise number and atom_counter
    number = 0
    atom_counter = 0

    # read the atom positions
    for line in structure_file:

        atom_counter += 1

        numbers.append(number)
        positions.append([float(el) for el in line.split()[:3]])

        if atom_counter == atoms_list[number]:
            # all atoms of the type corresponding to number have been read
            number += 1
            atom_counter = 0

structure_file.close()

cell  = (lattice, positions, numbers)

spacegroup = sl.get_spacegroup(cell)

if spacegroup.startswith("P"):
    # simple cubic
    prefactor = 1
elif spacegroup.startswith("F"):
    # FCC
    prefactor = 0.5
elif spacegroup.startswith("I"):
    # BCC
    prefactor = 0.25

# define x and y data
volumes = []
energies = []

# read energy and volume out of the energy.txt file
with open("energy.txt", "r") as energyfile:

    for line in energyfile:
        volumes.append(float(line.split()[1])) # volume
        energies.append(float(line.split()[2])) # energy

parameters, covariance = curve_fit(birch_murnaghan_eos, volumes, energies,
bounds=((2*min(energies), volumes[0], 0, -np.inf), (max(energies), volumes[-1], np.inf, np.inf)),
maxfev=100000) # impose very general boundaries to ensure correct result
# if the sampling in geometry_optimisation_cubic is anywhere near reasonable,
# this reliably gives the right result

with open("birch_murnaghan_eos_parameters.txt", "w") as outfile:

    outfile.write("E_0 = {} eV\n".format(parameters[0]))
    outfile.write("V_0 = {} A^3\n".format(parameters[1]))
    outfile.write("lattice parameter a = {} A\n".format((prefactor*np.abs(parameters[1]))**(1.0/3.0)))
    outfile.write("B_0 = {} eV/A^3 = {} GPa\n".format(parameters[2], parameters[2]*1.60217662*10*10)) # conversion eV/A^3 to GPa
    outfile.write("B_prime = {} [dimensionless]".format(parameters[3]))

# add writing of gnuplot file to plot fit
with open("energy_volume.gnu", "w") as plotfile:

    plotfile.write("set terminal postscript eps colour\n")
    plotfile.write("set style data points\n")
    plotfile.write("set key box lt -1 lw 2 width 2 height 1.5 opaque\n")
    plotfile.write("set output '| epstopdf --filter --outfile=energy_volume.pdf'\n")
    plotfile.write("E_0 = {}\n".format(parameters[0]))
    plotfile.write("V_0 = {}\n".format(parameters[1]))
    plotfile.write("B_0 = {}\n".format(parameters[2]))
    plotfile.write("B_prime = {}\n".format(parameters[3]))
    plotfile.write("birch_murnaghan_eos(x) = E_0 + 9.0*V_0*B_0/(16.0)*(((V_0/x)**(2.0/3.0)-1)**(3)*B_prime + ((V_0/x)**(2.0/3.0)-1)**(2)*(6-4*(V_0/x)**(2.0/3.0)))\n")
    plotfile.write("set xrange [{}:{}]\n".format(math.floor(min(volumes)), math.ceil(max(volumes))))
    plotfile.write("set xlabel 'Volume [A^3]'\n")
    plotfile.write("set ylabel 'Energy [eV]'\n")
    plotfile.write("plot 'energy.txt' u 2:3 w points pt 7 ps 1.5 lc rgb '#DC143C' title 'Energy-Volume data from DFT calculation', birch_murnaghan_eos(x) lw 2 lc rgb '#DC143C' title 'Least-Squares fit to Birch-Murnaghan EOS'")
