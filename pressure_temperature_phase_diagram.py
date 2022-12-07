# script to plot the pressure-temperature phase diagram, taking into account that phonons shift the pressure

# when calculating the contribution nuclear vibrations make to the energy, we must also take into account that the pressure is changed
# compared to the static lattice calculations by the vibrational contribution

# this is done by noting the P = -dE/dV
# so the pressure is found by fitting a polynomial to the energy-volume curve, where the energy includes the ZPE
# the derivatives at the different volumes where we have explicit calculations are then easily found

# implemented to work with my calculate_free_energy_wobble script

# written by Pascal Salzbrenner, pts28@cam.ac.uk

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from natsort import natsorted
from numpy.polynomial import Polynomial

############################################################ helper functions ############################################################

def connect_boundary_list(boundary_list, t_step):
    """Function that takes a list of points indicating the boundary between two phases and reshapes it so as to draw a connected boundary
       :param list boundary_list: a list of [pressure, temperature] points
       :param float t_step: the difference between temperatures in the input data

       :returns list connected_boundary_list: a list of [pressure, temperature] points, sorted to create a connected boundary
    """

    # sort by termperature
    boundary_list.sort(key = lambda x: x[1])

    # move the first element from the input to the output list
    connected_boundary_list = [boundary_list.pop(0)]

    # the way this will work is we always find the element closest to the already existing element by iterating over the entire list
    # not particularly elegant but it should work - I can't think of a situation where this is not the connectivity we want
    for i in range(len(boundary_list)):

        current_point = np.array(connected_boundary_list[-1])

        # set up the shortest distance to a large number so that the first attempt is guaranteed a hit
        shortest_distance = 100000

        for j in range(len(boundary_list)):

            point = np.array(boundary_list[j])

            # naively, or generally, we would like to find the absolute closest point, but this doesn't work when temperature is on a
            # much sparser scale than pressure - so instead we artifically restrict it to points neighbouring in temperature and take the one
            # with the closest pressure

            # also, to be fair, we do want to connect to a point adjacent in temperature

            if np.abs(point[1]-current_point[1]) < (t_step+1):

                distance = np.abs(point[0]-current_point[0])

                if distance < shortest_distance:
                    nearest_point_index = j
                    shortest_distance = distance

        connected_boundary_list.append(boundary_list.pop(nearest_point_index))

    return connected_boundary_list

def round_to_nearest_larger_five(number):
    """Function to round to the largest number divisible by 5*10^(num_digits-2) - eg 17.3 would be rounded to 20, 1167 to 1500, etc
       # will break if the number is smaller than 0, but at least for the application here this shouldn't occur
       :param float number: number to be rounded

       :returns int rounded_number: rounded number"""

    # convert number to string for a variety of useful operations
    str_number = str(number)

    # find the parts we do and don't round
    rounded_part = float(str_number[1:])
    str_rounded_part = str(rounded_part)
    not_rounded_remainder = int(number - rounded_part)

    # find the number of digits before the comma for the rounded part
    # also handle it if the rounded_part is smaller than 1

    if rounded_part >= 1:
        num_digits = len(str_rounded_part.split(".")[0])
    elif not_rounded_remainder > 0:
        # if not_rounded_remainder is an integer > 0, we always want to round to tenths
        num_digits = 0
    else:
        # the number of digits behind the comma corresponding to the first non-zero value is equal to the length before the comma of 1/rounded_part
        num_digits = 1 - len(str(int(1/rounded_part)))
    # find the multiple of 5 we round to
    multiple_of_five = 5*(10**(num_digits-1))

    # if this is smaller than the part we want to round, then we multiply by 2
    # ie we always want to round up, and the next-largest number divisible by the multiple of 5 is twice it
    # eg, if we have 1167, we want to round to 1500, but if we have 1667, we want to round to 2000
    if multiple_of_five < rounded_part:
        multiple_of_five *= 2

    # put the two parts together
    rounded_number = not_rounded_remainder + multiple_of_five

    return rounded_number

############################################################# input and setup #############################################################
# define pressure conversion factor from eV/A**3 to GPa
pressure_conversion = 160.21766208

# hardcode polynomial rank, but can be easily changed here
rank = 5

# hardcode step size for writing out the pressure-volume data, can again be changed easily here
pressure_increment = 0.01

# dictionary to contain the temperature - pressure - energy data for each structure
struc_temp_pressure_energies = {}

# list of temperatures
temp_list = []

################################################### reading pressure-volume-energy data ###################################################

ls_top = os.listdir()
ls_top = natsorted(ls_top)
temp_dirs = [temp_dir for temp_dir in ls_top if "temp_" in temp_dir]

# work out temperature step for phase boundary list sorting
t_step = float(temp_dirs[1].split("_")[1]) - float(temp_dirs[0].split("_")[1])

for temp_dir in temp_dirs:

    # we will consider all AIRSS .res files in the directory, grouped by the different indices associated with the different structures
    # we determine these indices here

    # set up dictionary to contain all files associated with each index/structure
    files_dict = {}

    # set up energy and volume dictionaries, as well as one for the static lattice pressure
    volumes = {}
    energies = {}
    full_energies = {}
    static_pressures = {}

    temperature = temp_dir.split("_")[1]
    temp_list.append(temperature)

    # determine all files in the directory
    ls = os.listdir(temp_dir)
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

    for structure, structure_files in files_dict.items():

        volumes[structure] = []
        energies[structure] = []
        static_pressures[structure] = []

        for file in structure_files:

            data_file = open("{}/{}".format(temp_dir, file), "r")

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
        full_energies[structure] = np.array(energies[structure] + pressures/pressure_conversion * volumes[structure])/natoms

        # fit a rank 5 polynomial to the pressure-energy data for interpolation
        pressure_energy_fit = Polynomial.fit(pressures, full_energies[structure], rank)

        pressure_fit_pressures = np.arange(pressures[-1],pressures[0],0.01)
        pressure_fit_energies = pressure_energy_fit(pressure_fit_pressures)
        plt.plot(pressures, full_energies[structure], 'o', pressure_fit_pressures, pressure_fit_energies, '-')
        plt.savefig("pressure_energy_fit_{}.pdf".format(structure))
        plt.close()

        # at this point we can flip the pressure to the order lowest -> highest for writing out
        pressures = np.flip(pressures)

########################################## write out pressure-volume data for subsequent plotting ##########################################
    
        # use the static pressures to delimit the pressure range

        initial_pressure = static_pressures[structure][0]
        final_pressure = static_pressures[structure][-1]

        pressure_energy_file = open("phonon_pressure_energy_{}.dat".format(structure), "w")

        pressure_energy_file.write("# Pressure [GPa]; Gibbs Free Energy [meV/atom]\n")

        while initial_pressure <= final_pressure:
            pressure_energy_file.write("{} {}\n".format(initial_pressure, 1000*pressure_energy_fit(initial_pressure)))

            initial_pressure += pressure_increment

        pressure_energy_file.close()

        # we reset the initial pressure here, so we have access to it when constructing the phase diagram
        # these should be the same for every structure anyways
        initial_pressure = static_pressures[structure][0]

######################################### write out the static and corresponding phonon pressures #########################################

        # open file to write the old and new pressures to
        pressure_file = open("static_phonon_pressure_{}.dat".format(structure), "w")
        pressure_file.write("# rank of polynomial: {}\n".format(rank))
        pressure_file.write("# static-lattice pressure [GPa]; vibration-corrected pressure [GPa]\n")

        # copy over the files we have operated on, with the correct pressure replacing that of the static lattice
        for i in range(len(structure_files)):
            # the structure_files and pressures lists will be of the same length, and the elements in the same place will correspond to one another
            pressure_file.write("{} {}\n".format(static_pressures[structure][i], pressures[i]))

        pressure_file.close()

        # store pressure_energy_fit
        if structure not in struc_temp_pressure_energies.keys():
            struc_temp_pressure_energies[structure] = {temperature: pressure_energy_fit}
        else:
            # not the first temperature for this structure
            struc_temp_pressure_energies[structure][temperature] = pressure_energy_fit

####################################################### generation of phase diagram #######################################################

# find lowest-energy structure for each temperature-pressure pair, as well as boundary between them
phase_diagram_data = []

structure_list = list(struc_temp_pressure_energies.keys())

# write out list of structures if it does not exist yet - if it does exist, we use it as an input for naming the structures

if not "structures_indices.dat" in os.listdir():
    with open("structures_indices.dat", "w") as index_file:
        for i in range(len(structure_list)):
            index_file.write("{}: {}\n".format(i, structure_list[i]))
else:
    with open("structures_indices.dat", "r") as index_file:
        for line in index_file:
            contents = line.split()
            index = int(contents[0].rstrip(":"))
            structure_name = contents[1]
            structure_list[index] = structure_name

# output file
pd_data_file = open("phase_diagram_data.dat", "w")
pd_data_file.write("# Pressure [GPa]; Temperature [K]; Ground state structure\n")

for temperature in temp_list:

    initial_pressure_copy = initial_pressure

    while initial_pressure_copy <= final_pressure:
        
        energies = []

        for structure in struc_temp_pressure_energies.keys():
            energies.append(struc_temp_pressure_energies[structure][temperature](initial_pressure_copy))

        energies = np.array(energies)
        minimum_energy_index = np.argmin(energies)

        phase_diagram_data.append([initial_pressure_copy,temperature,minimum_energy_index,structure_list[minimum_energy_index]])
        pd_data_file.write("{} {} {}\n".format(initial_pressure_copy, temperature, structure_list[minimum_energy_index]))

        initial_pressure_copy += pressure_increment

pd_data_file.close()

# determine phase transition (pressure, temperature) points
# read out the data corresponding to the first point to start us off, then iterate over the rest
# we write to file as well as storing these values

# make this a dictionary, where the keys are the indices of the structure
# as this is purely for plotting, it does not matter, for example, whether this is 2-3 or 3-2, so we group those together

phase_transition_points = {}
pt_points_file = open("phase_transition_points.dat", "w")
pt_points_file.write("# Points where a phase transition takes place\n")
pt_points_file.write("# Structure_1-Structure_2; Pressure [GPa]; Temperature [K]\n")

previous_pressure = phase_diagram_data[0][0]
previous_temperature = phase_diagram_data[0][1]
previous_mindex = phase_diagram_data[0][2]

lowest_pressure_prev_temp = previous_pressure
lowest_pressure_prev_temp_mindex = previous_mindex

for point in phase_diagram_data[1:]:

    pressure = point[0]
    temperature = point[1]
    mindex = point[2]

    if temperature != previous_temperature:
        # phase_diagram_data is structured such that all the points for a given temperature are grouped together - this detects that we have shifted to the next block
        # thus, we need to compare not to the previous point, but the the lowest pressure point for the previous temperature

        if mindex != lowest_pressure_prev_temp_mindex:
            # temperature has caused a phase transition here - we place it halfway between the two temperatures
            index_str="{}-{}".format(structure_list[min(lowest_pressure_prev_temp_mindex, mindex)], structure_list[max(lowest_pressure_prev_temp_mindex,mindex)])

            if index_str not in phase_transition_points.keys():
                phase_transition_points[index_str] = [pressure, (float(temperature)+float(previous_temperature))/2]
            else:
                phase_transition_points[index_str].append([pressure, (float(temperature)+float(previous_temperature))/2])

            pt_points_file.write("{} {} {}\n".format(index_str, pressure, (float(temperature)+float(previous_temperature))/2))

        previous_pressure = pressure
        previous_temperature = temperature
        previous_mindex = mindex

        lowest_pressure_prev_temp = previous_pressure
        lowest_pressure_prev_temp_mindex = previous_mindex

        continue

    if mindex != previous_mindex:
        # phase transition point
        index_str="{}-{}".format(structure_list[min(previous_mindex, mindex)], structure_list[max(previous_mindex,mindex)])

        if index_str not in phase_transition_points.keys():
            phase_transition_points[index_str] = [[(pressure+previous_pressure)/2, float(temperature)]]
        else:
            phase_transition_points[index_str].append([(pressure+previous_pressure)/2, float(temperature)])

        pt_points_file.write("{} {} {}\n".format(index_str, (pressure+previous_pressure)/2, float(temperature)))

    previous_pressure = pressure
    previous_temperature = temperature
    previous_mindex = mindex

pt_points_file.close()

########################################################### plotting #####################################################################

# plot with just the lines

plt.xlabel("Pressure [GPa]")
plt.ylabel("Temperature [K]")

# determine maximum temperature and pressure to set nice upper boundaries for the plot
# set low starting values
max_press = 0
max_temp = 0

for index_str, pt_line in phase_transition_points.items():

    # check if we have found new maximum values
    current_max = np.amax(pt_line, axis=0)
    current_max_press = current_max[0]
    current_max_temp = current_max[1]

    if current_max_press > max_press:
        max_press = current_max_press
    if current_max_temp > max_temp:
        max_temp = current_max_temp

    x, y = zip(*connect_boundary_list(pt_line, t_step))

    plt.plot(x, y, "#000080")
    plt.text(x[0], y[0], index_str)

# set plot parameters
plt.xlim(0, round_to_nearest_larger_five(max_press))
plt.ylim(0, round_to_nearest_larger_five(max_temp))

plt.savefig("phase_diagram_boundaries_only.pdf")

# plot with regions coloured according to different structures
for index_str, pt_line in phase_transition_points.items():

    x, y = zip(*connect_boundary_list(pt_line, t_step))
    plt.plot(x, y, "black")



plt.savefig("phase_diagram.pdf")



