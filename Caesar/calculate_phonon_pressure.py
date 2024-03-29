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
from scipy import interpolate
from natsort import natsorted
from numpy.polynomial import Polynomial

############################################################# input and setup #############################################################

# define pressure conversion factor from eV/A**3 to GPa
pressure_conversion = 160.21766208

# hardcode polynomial rank, but can be easily changed here
rank = 5

# pressure increment set to None, will be set below
pressure_increment = None

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
ls = natsorted(ls)

for file in ls:
    if ".res" in file:
        structure_name = file.split("p0")[1].lstrip("-").replace(".res", "")

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

    # we also require the number of atoms in the unit cell
    # this will be the same for every file, so it is not a problem that it's read every time
    reference_structure_natoms = int(important_line[7])

    data_file.close()

# flip the lists as we require the volumes to be sorted from lowest (corresponding to the highest pressure) to highest (vice versa)
volumes[reference_structure] = np.flip(np.array(volumes[reference_structure]))
energies[reference_structure] = np.flip(np.array(energies[reference_structure]))

splines = interpolate.splrep(volumes[reference_structure], energies[reference_structure], s=0)
fit_volumes = np.arange(volumes[reference_structure][0], volumes[reference_structure][-1],0.001)
fit_energies = interpolate.splev(fit_volumes, splines, der=0)

plt.plot(volumes[reference_structure], energies[reference_structure], 'o', fit_volumes, fit_energies, '-')
plt.savefig("energy_volume_polynomial_fit_{}.pdf".format(reference_structure))
plt.close()

reference_pressures = interpolate.splev(volumes[reference_structure], splines, der=1)
reference_pressures = -pressure_conversion*np.array(reference_pressures)

# find the Gibbs energies per atom using the correct pressure - add the PV term back in
# pressures are currently sorted from highest to lowest - this corresponds to the order of the energies and volumes
full_energies[reference_structure] = np.array(energies[reference_structure]
                                     + reference_pressures/pressure_conversion * volumes[reference_structure])/reference_structure_natoms

# fit a rank 5 polynomial to the pressure-energy data for interpolation
reference_pressure_energy_fit = Polynomial.fit(reference_pressures, full_energies[reference_structure], rank)

pressure_fit_pressures = np.arange(reference_pressures[-1],reference_pressures[0],0.01)
pressure_fit_energies = reference_pressure_energy_fit(pressure_fit_pressures)
plt.plot(reference_pressures, full_energies[reference_structure], 'o', pressure_fit_pressures, pressure_fit_energies, '-')
plt.savefig("pressure_energy_fit_{}.pdf".format(reference_structure))
plt.close()

# we also create a new directory with the .res files, where the energy is substituted to be correct for the full pressure
if not ("correct_pressure_energy_resfiles" in os.listdir()):
    os.mkdir("correct_pressure_energy_resfiles")

for file in files_dict[reference_structure]:

    infile = open(file, "r")
    outfile = open("correct_pressure_energy_resfiles/{}".format(file), "w")

    # in the first line, we read the pressure and then subsititute the energy from the fit at that pressure
    data = infile.readline().split()
    pressure = float(data[2])
    data[4] = str(reference_pressure_energy_fit(pressure)*reference_structure_natoms) # the fit is per atom, so we multiply here
    outfile.write(" ".join(data)+"\n")

    # copy over the rest of the lines
    for line in infile:
        outfile.write(line)

    infile.close()
    outfile.close()

# at this point we can flip the reference_pressures to the order lowest -> highest for writing out
reference_pressures = np.flip(reference_pressures)

for structure, structure_files in files_dict.items():

    # the data for the reference structure has already been read and processed
    # only carry out these operations if we are dealing with another structure
    if structure != reference_structure:

        volumes[structure] = []
        energies[structure] = []
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

            natoms = int(important_line[7])

            data_file.close()

        # flip the lists as we require the volumes to be sorted from lowest (corresponding to the highest pressure) to highest (vice versa)
        volumes[structure] = np.flip(np.array(volumes[structure]))
        energies[structure] = np.flip(np.array(energies[structure]))

########################################### fitting the polynomial and finding the real pressure ###########################################

        # use cubic splines with no smoothing to fit to the energy-volume data - gives better precision than polynomial interpolation
        splines = interpolate.splrep(volumes[structure], energies[structure], s=0)
        fit_volumes = np.arange(volumes[structure][0], volumes[structure][-1],0.001)
        fit_energies = interpolate.splev(fit_volumes, splines, der=0)

        plt.plot(volumes[structure], energies[structure], 'o', fit_volumes, fit_energies, '-')
        plt.savefig("energy_volume_polynomial_fit_{}.pdf".format(structure))
        plt.close()

        # calculate the first derivative of the interpolation to get the pressure
        pressures = interpolate.splev(volumes[structure], splines, der=1)
        pressures = -pressure_conversion*np.array(pressures)


        # find the Gibbs energies per atom using the correct pressure - add the PV term back in
        # pressures are currently sorted from highest to lowest - this corresponds to the order of the energies and volumes
        full_energies[structure] = np.array(energies[structure]
                                           + pressures/pressure_conversion * volumes[structure])/natoms

        # fit a rank 5 polynomial to the pressure-energy data for interpolation
        pressure_energy_fit = Polynomial.fit(pressures, full_energies[structure], rank)

        pressure_fit_pressures = np.arange(pressures[-1],pressures[0],0.01)
        pressure_fit_energies = pressure_energy_fit(pressure_fit_pressures)
        plt.plot(pressures, full_energies[structure], 'o', pressure_fit_pressures, pressure_fit_energies, '-')
        plt.savefig("pressure_energy_fit_{}.pdf".format(structure))
        plt.close()

        # at this point we can flip the pressure to the order lowest -> highest for writing out
        pressures = np.flip(pressures)

        for file in structure_files:

            infile = open(file, "r")
            outfile = open("correct_pressure_energy_resfiles/{}".format(file), "w")

            # in the first line, we read the pressure and then subsititute the energy from the fit at that pressure
            data = infile.readline().split()
            pressure = float(data[2])
            data[4] = str(pressure_energy_fit(pressure)*natoms) # the fit is per atom, so we multiply here

            outfile.write(" ".join(data)+"\n")

            # copy over the rest of the lines
            for line in infile:
                outfile.write(line)

            infile.close()
            outfile.close()

########################################## write out pressure-volume data for subsequent plotting ##########################################

        # determine starting pressure - the higher of the lowest pressure for the reference structure and that of the current one
        if pressures[0] >= reference_pressures[0]:
            initial_pressure = pressures[0]
        else:
            initial_pressure = reference_pressures[0]

        # determine final pressure - the lower of the highest pressure for the reference structure and that of the current one
        if pressures[-1] <= reference_pressures[-1]:
            final_pressure = pressures[-1]
        else:
            final_pressure = reference_pressures[-1]

        if not pressure_increment:
            # set the pressure increment such that 10 000 steps are used to interpolate
            pressure_increment = 0.0001 * (final_pressure - initial_pressure)

        pressure_energy_file = open("phonon_pressure_energy_{}.dat".format(structure), "w")

        pressure_energy_file.write("# Pressure [GPa]; Gibbs Free Energy [meV/atom], relative to the Gibbs Free Energy of {}\n".format(reference_structure))

        while initial_pressure < final_pressure:
            pressure_energy_file.write("{} {}\n".format(initial_pressure,
            (1000*(pressure_energy_fit(initial_pressure))-1000*(reference_pressure_energy_fit(initial_pressure)))))

            initial_pressure += pressure_increment

        pressure_energy_file.close()

######################################### write out the static and corresponding phonon pressures #########################################

    else:
        # this is the reference structure - we need to set the pressures list here for the following part to work
        pressures = reference_pressures[:]

    # open file to write the old and new pressures to
    pressure_file = open("static_phonon_pressure_{}.dat".format(structure), "w")
    pressure_file.write("# rank of polynomial: {}\n".format(rank))
    pressure_file.write("# static-lattice pressure [GPa]; vibration-corrected pressure [GPa]\n")

    # copy over the files we have operated on, with the correct pressure replacing that of the static lattice
    for i in range(len(structure_files)):
        # the structure_files and pressures lists will be of the same length, and the elements in the same place will correspond to one another
        pressure_file.write("{} {}\n".format(static_pressures[structure][i], pressures[i]))

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
plotfile.write("set ylabel 'Gibbs Free Energy [meV/atom]'\n")

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
