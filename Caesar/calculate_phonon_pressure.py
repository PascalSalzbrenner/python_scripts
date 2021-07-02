# when calculating the contribution nuclear vibrations make to the energy, we must also take into account that the pressure is changed
# compared to the static lattice calculations by the vibrational contribution

# this is done by noting the P = -dE/dV
# so the pressure is found by fitting a polynomial to the energy-volume curve, where the energy includes the ZPE
# the derivatives at the different volumes where we have explicit calculations are then easily found

# implemented to work with my add_zpe_enthalpy.py script, but could probably be extended (somewhat) easily

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial import Polynomial
from scipy.interpolate import interp1d

############################################################# input and setup #############################################################

# define pressure conversion factor from eV/A**3 to GPa
pressure_conversion = 160.21766208

# hardcode polynomial rank, but can be easily changed here
rank = 5

# hardcode step size for writing out the pressure-volume data, can again be changed easily here
pressure_increment = 0.01

# read the reference structure - this must agree with the structure_name
reference_structure = sys.argv[1]

# we will consider all AIRSS .res files in the directory, grouped by the different indices associated with the different structures
# we determine these indices here

# set up dictionary to contain all files associated with each index/structure
files_dict = {}

# set up energy and volume dictionaries, as well as one for the static lattice pressure
volumes = {}
energies = {}
full_energies = {}
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

# reference structure is handled separately

volumes[reference_structure] = []
energies[reference_structure] = []
full_energies[reference_structure] = []
static_pressures[reference_structure] = []

for file in files_dict[reference_structure]:

    data_file = open(file, "r")

    # this is the AIRSS .res file format
    # the pressure [GPa], volume [Ang**3], and energy [eV] are the third, fourth, and fifth elements of the first line respectively
    important_line = data_file.readline().split()
    static_pressures[reference_structure].append(float(important_line[2]))
    volumes[reference_structure].append(float(important_line[3]))
    # we need to fit to the energy, not the enthalpy, so we remove the PV contribution
    energies[reference_structure].append(float(important_line[4])-
                                        (static_pressures[reference_structure][-1]/pressure_conversion)*volumes[reference_structure][-1])
    # however, we need the full energy when plotting the pressure-energy curve at the end
    full_energies[reference_structure].append(float(important_line[4]))

    # we also require the number of atoms in the unit cell
    # this will be the same for every file, so it is not a problem that it's read every time
    reference_structure_natoms = int(important_line[7])

    data_file.close()

volumes[reference_structure] = np.array(volumes[reference_structure])
energies[reference_structure] = np.array(energies[reference_structure])
full_energies[reference_structure] = np.array(full_energies[reference_structure])/reference_structure_natoms

# polynomial of degree 5 seems to work fine, but the fit is plotted for visual, and the residuals calculated for quantitative confirmation
# a more sophisticated way of doing this would be to iterate the order of the polynomial until the residuals are below some threshold
fit, additional_info = Polynomial.fit(volumes[reference_structure], energies[reference_structure], rank, full=True)
fit_volumes, fit_energies = fit.linspace(500)
plt.plot(volumes[reference_structure], energies[reference_structure], 'o', fit_volumes, fit_energies, '-')
plt.savefig("energy_volume_polynomial_fit_{}.pdf".format(reference_structure))
plt.close()

residual_sum = additional_info[0][0]

# find the first derivative - the pressure
first_derivative = fit.deriv(m=1)

# find the pressure at all the explicit volume points
reference_pressures = []

for volume in volumes[reference_structure]:
    reference_pressures.append(-pressure_conversion*first_derivative(volume))

reference_pressures.sort()

# cubically interpolate the pressure-energy data
reference_pressure_energy_interpolation = interp1d(reference_pressures, full_energies[reference_structure], kind="cubic")

pressure_fit_pressures = np.arange(reference_pressures[0],reference_pressures[-1],0.01)
pressure_fit_energies = reference_pressure_energy_interpolation(pressure_fit_pressures)
plt.plot(reference_pressures, full_energies[reference_structure], 'o', pressure_fit_pressures, pressure_fit_energies, '-')
plt.savefig("pressure_energy_interpolation_{}.pdf".format(reference_structure))
plt.close()

for structure, structure_files in files_dict.items():

    # the data for the reference structure has already been read and processed
    # only carry out these operations if we are dealing with another structure
    if structure != reference_structure:

        volumes[structure] = []
        energies[structure] = []
        full_energies[structure] = []
        static_pressures[structure] = []

        for file in structure_files:

            data_file = open(file, "r")

            # this is the AIRSS .res file format
            # the pressure [GPa], volume [Ang**3], and energy [eV] are the third, fourth, and fifth elements of the first line respectively
            important_line = data_file.readline().split()
            static_pressures[structure].append(float(important_line[2]))
            volumes[structure].append(float(important_line[3]))
            # we need to fit to the energy, not the enthalpy, so we remove the PV contribution
            energies[structure].append(float(important_line[4])-(static_pressures[structure][-1]/pressure_conversion)*volumes[structure][-1])
            full_energies[structure].append(float(important_line[4]))

            natoms = int(important_line[7])

            data_file.close()

        volumes[structure] = np.array(volumes[structure])
        energies[structure] = np.array(energies[structure])
        full_energies[structure] = np.array(full_energies[structure])/natoms

########################################### fitting the polynomial and finding the real pressure ###########################################

        # polynomial of degree 5 seems to work fine, but the fit is plotted for visual, and the residuals calculated for quantitative confirmation
        # a more sophisticated way of doing this would be to iterate the order of the polynomial until the residuals are below some threshold
        fit, additional_info = Polynomial.fit(volumes[structure], energies[structure], rank, full=True)
        fit_volumes, fit_energies = fit.linspace(500)
        plt.plot(volumes[structure], energies[structure], 'o', fit_volumes, fit_energies, '-')
        plt.savefig("energy_volume_polynomial_fit_{}.pdf".format(structure))
        plt.close()

        residual_sum = additional_info[0][0]

        # find the first derivative - the pressure
        first_derivative = fit.deriv(m=1)

        # find the pressure at all the explicit volume points
        pressures = []

        for volume in volumes[structure]:
            pressures.append(-pressure_conversion*first_derivative(volume))

        pressures.sort()

        # cubically interpolate the pressure-energy data here
        pressure_energy_interpolation = interp1d(pressures, full_energies[structure], kind="cubic")

        pressure_fit_pressures = np.arange(pressures[0],pressures[-1],0.01)
        pressure_fit_energies = pressure_energy_interpolation(pressure_fit_pressures)
        plt.plot(pressures, full_energies[structure], 'o', pressure_fit_pressures, pressure_fit_energies, '-')
        plt.savefig("pressure_energy_interpolation_{}.pdf".format(structure))
        plt.close()

########################################## write out pressure-volume data for subsequent plotting ##########################################

        # determine starting pressure - the higher of the lowest pressure for the reference structure and that of the current one
        if pressures[0] >= reference_pressures[0]:
            initial_pressure = pressures[0]
        else:
            initial_pressure = reference_pressures[0]

        # determine starting pressure - the lower of the highest pressure for the reference structure and that of the current one
        if pressures[-1] <= reference_pressures[-1]:
            final_pressure = pressures[-1]
        else:
            final_pressure = reference_pressures[-1]

        pressure_energy_file = open("phonon_pressure_energy_{}.dat".format(structure), "w")

        pressure_energy_file.write("# Pressure [GPa]; Gibbs Free Energy [meV/atom], relative to the Gibbs Free Energy of {}\n".format(reference_structure))

        while initial_pressure < final_pressure:
            pressure_energy_file.write("{} {}\n".format(initial_pressure,
            (1000*(pressure_energy_interpolation(initial_pressure))-1000*(reference_pressure_energy_interpolation(initial_pressure)))))

            initial_pressure += pressure_increment

        pressure_energy_file.close()

######################################### write out the static and corresponding phonon pressures #########################################

    else:
        # this is the reference structure - we need to set the pressures list here for the following part to work
        pressures = reference_pressures[:]

    # open file to write the old and new pressures to
    pressure_file = open("static_phonon_pressure_{}.dat".format(structure), "w")
    pressure_file.write("# rank of polynomial: {}\n".format(rank))
    pressure_file.write("# sum of energy-volume fit residuals: {}\n\n".format(residual_sum))
    pressure_file.write("# static-lattice pressure [GPa]; vibration-corrected pressure [GPa]\n")

    # copy over the files we have operated on, with the correct pressure replacing that of the static lattice
    for i in range(len(structure_files)):
        # the structure_files and pressures lists will be of the same length, and the elements in the same place will correspond to one another

        correct_pressure = pressures[i]
        relevant_file = structure_files[i]

        pressure_file.write("{} {}\n".format(static_pressures[structure][i], correct_pressure))

    pressure_file.close()

######################################################### generate Gnuplot script #########################################################

plotfile = open("phonon_pressure_energy.gnu", "w")

plotfile.write("set terminal postscript eps colour font 'Helvetica,20'\n")
plotfile.write("set style data lines\n")
plotfile.write("set output '| epstopdf --filter --outfile=phonon_pressure_energy.pdf'\n")

plotfile.write("set key top right\n")
plotfile.write("set key box lt -1 lw 2 width 2 height 1.5 opaque font 'Helvectica,15'\n")

# redefine gnuplot linetypes with nice colours
plotfile.write("set linetype 1 lc rgb '#DC143C'\n")
plotfile.write("set linetype 2 lc rgb '#D95F02'\n")
plotfile.write("set linetype 3 lc rgb '#E6AB02'\n")
plotfile.write("set linetype 4 lc rgb '#66A61E'\n")
plotfile.write("set linetype 5 lc rgb '#8000C4'\n")
plotfile.write("set linetype 6 lc rgb '#7570B3'\n")
plotfile.write("set linetype 7 lc rgb '#E7298A'\n")
plotfile.write("set linetype 8 lc rgb '#1E90FF'\n")
plotfile.write("set linetype 9 lc rgb '#1B9E77'\n")
plotfile.write("set linetype 10 lc rgb '#B8860B'\n")
plotfile.write("set linetype 11 lc rgb '#20C2C2'\n")
plotfile.write("set linetype cycle 11\n")

plotfile.write("set xlabel 'Pressure [GPa]'\n")
plotfile.write("set ylabel 'Enthalpy [meV/atom]'\n")

# pick somewhat sensible defaults for the xrange using the fact that the reference structure must span the whole range at the static level
plotfile.write("set xrange [{}:{}]\n".format(int(static_pressures[reference_structure][0]),
                                             int(static_pressures[reference_structure][-1])+100))
plotfile.write("set xtics 50\n")
plotfile.write("set mxtics 2\n")

# generate plot command
plot_string = "plot 0 w lines lw 5 title '{}'".format(reference_structure)

for structure in files_dict.keys():

    # skip reference structure
    if structure != reference_structure:
        plot_string += ", 'phonon_pressure_energy_{0}.dat' u 1:2 w lines lw 5 title '{0}'".format(structure)

plotfile.write("{}\n".format(plot_string))

plotfile.close()
