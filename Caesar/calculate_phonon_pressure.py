# when calculating the contribution nuclear vibrations make to the energy, we must also take into account that the pressure is changed
# compared to the static lattice calculations by the vibrational contribution

# this is done by noting the P = -dE/dV
# so the pressure is found by fitting a polynomial to the energy-volume curve, where the energy includes the ZPE
# the derivatives at the different volumes where we have explicit calculations are then easily found

# implemented to work with my add_zpe_enthalpy.py script, but could probably be extended (somewhat) easily

import os
import re
import sys
import shutil
import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial import Polynomial

############################################################# input and setup #############################################################

# define pressure conversion factor from eV/A**3 to GPa
pressure_conversion = 160.21766208

# hardcode polynomial rank, but can be easily changed here
rank = 5

# compile a regex pattern to correctly recognise the pressure in a .res filename
pressure_pattern = re.compile(r"_\d*p0")

# delete and set up fresh enthalpy_plot_ZPE_correct_pressure
shutil.rmtree("../enthalpy_plot_ZPE_correct_pressure", ignore_errors=True)
os.mkdir("../enthalpy_plot_ZPE_correct_pressure")

# read the reference structure - this must agree with the structure_name
reference_structure = sys.argv[1]

# we will consider all AIRSS .res files in the directory, grouped by the different indices associated with the different structures
# we determine these indices here

# set up dictionary to contain all files associated with each index/structure
files_dict = {}

# set up energy and volume dictionaries, as well as one for the static lattice pressure
volumes = {}
energies = {}
static_pressures = {}

################################################### reading pressure-volume-energy data ###################################################

# determine all files in the directory
ls = os.listdir()
ls.sort()

for file in ls:
    if ".res" in file:
        structure_name = file.split("p0")[1].lstrip("-").rstrip(".res")

        if structure_name not in files_dict.keys():
            # this is the first structure with this index, so we have to generate the appropriate list in the dictionary
            files_dict[structure_name] = [file]
        else:
            # the appropriate list already exists
            files_dict[structure_name].append(file)

# iterate over and read all files for each structure
for structure, structure_files in files_dict.items():

    volumes[structure] = []
    energies[structure] = []
    static_pressures[structure] = []

    for file in structure_files:

        data_file = open(file, "r")

        # this is the AIRSS .res file format, so the pressure [GPa], volume [Ang**3], and energy [eV] are the third, fourth, and fifth elements of the first line respectively
        important_line = data_file.readline().split()
        static_pressures[structure].append(float(important_line[2]))
        volumes[structure].append(float(important_line[3]))
        # we need to fit to the energy, not the enthalpy, so we remove the PV contribution
        energies[structure].append(float(important_line[4])-(static_pressures[structure][-1]/pressure_conversion)*volumes[structure][-1])

        data_file.close()

    volumes[structure] = np.array(volumes[structure])
    energies[structure] = np.array(energies[structure])

########################################### fitting the polynomial and finding the real pressure ###########################################

    # polynomial of degree 5 seems to work fine, but the fit is plotted for visual, and the residuals calculated for quantitative confirmation
    # a more sophisticated way of doing this would be to iterate the order of the polynomial until the residuals are below some threshold
    fit, additional_info = Polynomial.fit(volumes[structure], energies[structure], rank, full=True)
    fit_volumes, fit_energies = fit.linspace(500)
    plt.plot(volumes[structure], energies[structure], 'o', fit_volumes, fit_energies, '-')
    plt.savefig("../enthalpy_plot_ZPE_correct_pressure/energy_volume_polynomial_fit_{}.pdf".format(structure))
    plt.close()

    residual_sum = additional_info[0][0]

    # find the first derivative - the pressure
    first_derivative = fit.deriv(m=1)

    # find the pressure at all the explicit volume points
    pressures = []

    for volume in volumes[structure]:
        pressures.append(-pressure_conversion*first_derivative(volume))

    pressures.sort()

####################################### writing generated data for plotting and other postprocessing #######################################

    # open file to write the old and new pressures to
    pressure_file = open("../enthalpy_plot_ZPE_correct_pressure/static_phonon_pressure_{}.dat".format(structure), "w")
    pressure_file.write("# rank of polynomial: {}\n".format(rank))
    pressure_file.write("# sum of energy-volume fit residuals: {}\n\n".format(residual_sum))
    pressure_file.write("# static-lattice pressure [GPa]; vibration-corrected pressure [GPa]\n")

    # copy over the files we have operated on, with the correct pressure replacing that of the static lattice
    for i in range(len(structure_files)):
        # the structure_files and pressures lists will be of the same length, and the elements in the same place will correspond to one another

        correct_pressure = pressures[i]
        relevant_file = structure_files[i]

        pressure_file.write("{} {}\n".format(static_pressures[structure][i], correct_pressure))

        read_file = open("{}".format(relevant_file), "r")
        write_file = open("../enthalpy_plot_ZPE_correct_pressure/{}".format(pressure_pattern.sub("_{}p0".format(int(correct_pressure)),
                                                                                                                relevant_file)), "w")

        # the first line is the only one that must be changed
        first_line = read_file.readline().split()

        # the second element contains the pressure
        first_line[1] = pressure_pattern.sub("_{}p0".format(int(correct_pressure)), first_line[1])

        # third element is the pressure
        first_line[2] = str(correct_pressure)

        write_file.write(" ".join(first_line))
        write_file.write("\n")

        # rest of the lines are copied over as they are
        for line in read_file:
            write_file.write(line)

        read_file.close()
        write_file.close()

    pressure_file.close()

# TODO implement plotting in Gnuplot
# important features:
# data files for each structure - pressure in column 1, energy per atom in column 2
# calculate the pressure on a dense enough spacing - perhaps 0.01 GPa (but can play around with)
# I can use a polynomial fit for this
# the points are not as smoothly distributed as at the static lattice, so this will erase some information, but I feel like it's still the best option
# use the same range (the range of the reference structure, which will become an input) for all structures to enable subtraction
# if there are structures with a narrower pressure range, all values not in its range will be set to 0 - when the
