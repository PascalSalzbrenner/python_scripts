# script to plot the pressure-temperature phase diagram, taking into account that phonons shift the pressure

# when calculating the contribution nuclear vibrations make to the energy, we must also take into account that the pressure is changed
# compared to the static lattice calculations by the vibrational contribution

# this is done by noting that P = -dE/dV
# so the pressure is found by fitting a polynomial to the energy-volume curve, where the energy includes the ZPE
# the derivatives at the different volumes where we have explicit calculations are then easily found

# implemented to work with my calculate_free_energy_wobble script

# written by Pascal Salzbrenner, pts28@cam.ac.uk
# with thanks to Se Hun Joo for introducing me to pcolormesh

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

from copy import deepcopy
from scipy import interpolate
from natsort import natsorted
from numpy.polynomial import Polynomial
from matplotlib.colors import ListedColormap

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

            if np.abs(point[1]-current_point[1]) < (2*t_step-t_step/10):

                distance = np.abs(point[0]-current_point[0])

                if distance < shortest_distance:
                    nearest_point_index = j
                    shortest_distance = distance

        try:
            connected_boundary_list.append(boundary_list.pop(nearest_point_index))
        except IndexError:
            return connected_boundary_list            

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

# hardcode step size for writing out the pressure-volume data as well as temperatures, can again be changed easily here
pressure_increment = 0.01
t_step = 5

# dictionaries to contain the temperature - pressure / volume - energy data for each structure, input as well as output
struc_temp_input_data = {}
struc_temp_pressure_energy_fits = {}

################################################### reading pressure-volume-energy data ###################################################

print("reading input data")

ls_top = os.listdir()
ls_top = natsorted(ls_top)
temp_dirs = [temp_dir for temp_dir in ls_top if "temp_" in temp_dir]

for temp_dir in temp_dirs:

    # we will consider all AIRSS .res files in the directory, grouped by the different indices associated with the different structures
    # we determine these indices here

    # set up dictionary to contain all files associated with each index/structure
    files_dict = {}

    temperature = temp_dir.split("_")[1]

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

        # check if the structure is already in the input dictionary; add it if not
        if structure not in struc_temp_input_data.keys():
            struc_temp_input_data[structure] = {}

        volumes = []
        energies = []
        static_pressures = []

        for file in structure_files:

            data_file = open("{}/{}".format(temp_dir, file), "r")

            # this is the AIRSS .res file format
            # the pressure [GPa], volume [Ang**3], and energy [eV] are the third, fourth, and fifth elements of the first line respectively
            important_line = data_file.readline().split()
            natoms = int(important_line[7])
            static_pressures.append(float(important_line[2]))
            # volumes and energies are both calculated per atom to enable comparison between structures
            volumes.append(float(important_line[3])/natoms)
            # we need to fit to the energy, not the enthalpy, so we remove the PV contribution
            energies.append((float(important_line[4])/natoms-(static_pressures[-1]/pressure_conversion)*volumes[-1]))

            data_file.close()

        # flip the lists as we require the volumes to be sorted from lowest (corresponding to the highest pressure) to highest (vice versa)
        volumes = np.flip(np.array(volumes))
        energies = np.flip(np.array(energies))
        static_pressures = np.flip(np.array(static_pressures))

        struc_temp_input_data[structure][temperature] = [volumes[:], energies[:], static_pressures[:]]

############################ interpolate the temperature-free energy curves for every volume with a polynomial ############################

print("Interpolating temperature-energy curves")

for structure in struc_temp_input_data.keys():

    temperatures = []
    energies = []
    fits = []

    for temperature in struc_temp_input_data[structure].keys():
        temperatures.append(float(temperature))

        # note that energies is itself a list, so for every temperature we have several energy values corresponding to the different volumes
        # we do one fit for each entry, the order of the fits corresponds to that of the volumes
        energies.append(struc_temp_input_data[structure][temperature][1])

        # the volumes and static pressures will be the same for every temperature - here we read them out for later use
        volumes = struc_temp_input_data[structure][temperature][0]
        static_pressures = struc_temp_input_data[structure][temperature][2]

    energies = np.array(energies)

    # do the fits
    for i in range(len(energies[0])):
        temperature_energy_fit = Polynomial.fit(temperatures, energies[:,i], rank)
        fits.append(temperature_energy_fit)

        # plot fit to assess its quality
        temperature_fit_temperatures = np.arange(temperatures[0], temperatures[-1], t_step)
        temperature_fit_energies = temperature_energy_fit(temperature_fit_temperatures)
        plt.plot(temperatures, energies[:,i], 'o', temperature_fit_temperatures, temperature_fit_energies, '-')
        plt.savefig("temperature_energy_fit_{}_volume_{}.pdf".format(structure, volumes[i]))
        plt.close()

    # write fits to input dict
    initial_temperature = temperatures[0]
    final_temperature = temperatures[-1]

    while initial_temperature < final_temperature:

        # check if we already have an entry for this temperature

        if not np.isclose(initial_temperature, temperatures, atol=t_step/10).any():

            fitted_energies = []

            for fit in fits:
                fitted_energies.append(fit(initial_temperature))

            struc_temp_input_data[structure][str(initial_temperature)] = [volumes[:], fitted_energies[:], static_pressures[:]]

        initial_temperature += t_step

########################################### fitting the polynomial and finding the real pressure ###########################################

print("Calculating phonon-modified pressure and interpolating it")

# set up dictionary to contain output data
full_energies = {}

# set up list of temperatures, both as strings and as numbers
temp_list = []
temp_list_numbers = []

for structure in struc_temp_input_data.keys():

    for temperature in struc_temp_input_data[structure].keys():

        volumes = struc_temp_input_data[structure][temperature][0]
        energies = struc_temp_input_data[structure][temperature][1]
        static_pressures = struc_temp_input_data[structure][temperature][2]

        # use cubic splines with no smoothing to fit to the energy-volume data - gives better precision than polynomial interpolation
        splines = interpolate.splrep(volumes, energies, s=0)
        fit_volumes = np.arange(volumes[0], volumes[-1], 0.001)
        fit_energies = interpolate.splev(fit_volumes, splines, der=0)

        # output plots every 100 K
        if np.isclose(float(temperature) % 100, 0):

            plt.plot(volumes, energies, 'o', fit_volumes, fit_energies, '-')
            plt.savefig("energy_volume_polynomial_fit_{}_{}_K.pdf".format(structure, temperature))
            plt.close()

        # calculate the first derivative of the interpolation to get the pressure
        pressures = interpolate.splev(volumes, splines, der=1)
        pressures = -pressure_conversion*np.array(pressures)

        # find the Gibbs energies per atom using the correct pressure - add the PV term back in
        # pressures are currently sorted from highest to lowest - this corresponds to the order of the energies and volumes
        full_energies[structure] = np.array(energies + pressures/pressure_conversion * volumes)

        # fit a rank 5 polynomial to the pressure-energy data for interpolation
        pressure_energy_fit = Polynomial.fit(pressures, full_energies[structure], rank)

        pressure_fit_pressures = np.arange(pressures[-1],pressures[0],0.01)
        pressure_fit_energies = pressure_energy_fit(pressure_fit_pressures)

        # output plots every 100 K
        if np.isclose(float(temperature) % 100, 0):
            
            plt.plot(pressures, full_energies[structure], 'o', pressure_fit_pressures, pressure_fit_energies, '-')
            plt.savefig("pressure_energy_fit_{}_{}_K.pdf".format(structure, temperature))
            plt.close()

########################################## write out pressure-volume data for subsequent plotting ##########################################

        print("writing out pressure-volume data")

        # use the static pressures to delimit the pressure range

        initial_pressure = static_pressures[-1]
        final_pressure = static_pressures[0]

        # out data every 100 K - should be sufficient to interpolate if I ever need it in a different script
        if np.isclose(float(temperature) % 100, 0):

            pressure_energy_file = open("phonon_pressure_energy_{}_{}_K.dat".format(structure, temperature), "w")

            pressure_energy_file.write("# Pressure [GPa]; Gibbs Free Energy [meV/atom]\n")

            while initial_pressure <= final_pressure:
                pressure_energy_file.write("{} {}\n".format(initial_pressure, 1000*pressure_energy_fit(initial_pressure)))

                initial_pressure += pressure_increment

            pressure_energy_file.close()

            # we reset the initial pressure here, so we have access to it when constructing the phase diagram
            # these should be the same for every structure anyways
            initial_pressure = static_pressures[-1]

######################################### write out the static and corresponding phonon pressures #########################################

        print("writing out static and phonon-modified pressures")
        
        # output presssure every 100 K to diagnose the trend
        if np.isclose(float(temperature) % 100, 0):
        
            # open file to write the old and new pressures to
            pressure_file = open("static_phonon_pressure_{}_{}_K.dat".format(structure, temperature), "w")
            pressure_file.write("# rank of polynomial: {}\n".format(rank))
            pressure_file.write("# static-lattice pressure [GPa]; vibration-corrected pressure [GPa]\n")

            # copy over the files we have operated on, with the correct pressure replacing that of the static lattice
            for i in range(len(files_dict[structure])):
                # the structure_files and pressures lists will be of the same length, and the elements in the same place will correspond to one another
                pressure_file.write("{} {}\n".format(static_pressures[i], pressures[i]))

            pressure_file.close()

        # store pressure_energy_fit
        if structure not in struc_temp_pressure_energy_fits.keys():
            struc_temp_pressure_energy_fits[structure] = {temperature: pressure_energy_fit}
        else:
            # not the first temperature for this structure
            struc_temp_pressure_energy_fits[structure][temperature] = pressure_energy_fit

        if temperature not in temp_list:
            # add the temperature to the list of temperatures - needed for reading the data out correctly
            temp_list.append(temperature)
            temp_list_numbers.append(float(temperature))

# sort temp lists
temp_list = natsorted(temp_list)
temp_list_numbers.sort()

####################################################### generation of phase diagram #######################################################

print("generating phase diagram")

# find lowest-energy structure for each temperature-pressure pair, as well as boundary between them
phase_diagram_data = []

structure_list = list(struc_temp_pressure_energy_fits.keys())

# write out list of structures if it does not exist yet - if it does exist, we use it as an input for naming the structures
# this must NOT contain the data for the liquid, it will mess things up if it does!

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

# create list to hold pressures
pressure_list = []

for temperature in temp_list:

    initial_pressure_copy = initial_pressure

    while initial_pressure_copy < final_pressure + 0.1 * pressure_increment:
        # sometimes a bit of error accumulates when adding the pressure increment, which can lead to the last step being slightly larger than final_pressure
        # that's why we compare to a number slightly larger than final_pressure, but strictly smaller than final_pressure + pressure_increment

        if len(pressure_list) < round(((final_pressure-initial_pressure)/pressure_increment)+1):
            pressure_list.append(initial_pressure_copy)
        
        energies = []

        for structure in struc_temp_pressure_energy_fits.keys():
            energies.append(struc_temp_pressure_energy_fits[structure][temperature](initial_pressure_copy))

        energies = np.array(energies)
        minimum_energy_index = np.argmin(energies)

        phase_diagram_data.append([initial_pressure_copy,temperature,minimum_energy_index,structure_list[minimum_energy_index]])

        initial_pressure_copy += pressure_increment

################################ check if melting data has been given; if so, adapt the phase diagram data ################################

if "melt_curve.dat" in ls_top:

    print("Reading in melt curve")

    # read pressure-temperature data for the melt curve and fit a polynomial for interpolation
    melt_pressures = []
    melt_temperatures = []

    # assign index to melt which is the next unused one
    melt_index = len(structure_list)

    # put "liquid" into structure_list
    structure_list.append("liquid")

    with open("melt_curve.dat", "r") as melt_file:

        # skip first line
        melt_file.readline()

        # read data from other lines
        for line in melt_file:

            data=line.split()

            melt_pressures.append(float(data[0]))
            melt_temperatures.append(float(data[1]))

    # fit melt curve
    melt_curve = Polynomial.fit(melt_pressures, melt_temperatures, rank)

    melt_curve_fit_pressures = np.arange(melt_pressures[0],melt_pressures[-1],0.01)
    melt_curve_fit_temperatures = melt_curve(melt_curve_fit_pressures)

    # plot melt curve to make sure everything has worked
    plt.plot(melt_pressures, melt_temperatures, 'o', melt_curve_fit_pressures, melt_curve_fit_temperatures, '-')
    plt.savefig("melt_curve.pdf")
    plt.close()

    # iterate over phase diagram data
    # at every point, for a given pressure, where temperature > melt temperature, we change the index to that for the liquid
    for i in range(len(phase_diagram_data)):

        if float(phase_diagram_data[i][1]) > melt_curve(phase_diagram_data[i][0]):
            phase_diagram_data[i][2] = melt_index
            phase_diagram_data[i][3] = "liquid"

########################################################### Put data into plottable form #################################################

print("Putting data into plottable form")

# create list to hold indices
minimum_index_list = []

# write phase diagram data to output file - format such that reading can be done with pandas
pd_data_file = open("phase_diagram_data.dat", "w")
pd_data_file.write("Pressure [GPa], Temperature [K], Ground state structure\n")

for point in phase_diagram_data:
    minimum_index_list.append(point[2])
    pd_data_file.write("{}, {}, {}\n".format(point[0], point[1], point[3]))

pd_data_file.close()

# reshape the minimum_index_list to pass it to pandas
minimum_index_array=np.array(minimum_index_list).reshape([len(temp_list), len(pressure_list)])

##################################################### determine the phase boundaries #####################################################

print("Determining the phase boundaries")

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
                phase_transition_points[index_str] = [[pressure, (float(temperature)+float(previous_temperature))/2]]
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

# treatment for extremely obscure problem I faced with lead
# very near the phase transition for a tiny bit it's first Fm-3m-liquid, then P63/mmc-liquid, then Fm-3m-liquid
# this leads to regions where visually, there are two melt boundaries - to avoid this, we split lists where successive points
# are separated by more than t_step (as the outer loop is over temperature, this is what implies a discontinuous line)

split_phase_transition_points = {}

for index_str, pt_line in phase_transition_points.items():

    # sort in order of ascending pressure
    pt_line.sort(key = lambda x: x[0])

    split_pt_lines = []
    split_pt_line = [pt_line[0]]

    for pt_point in pt_line[1:]:
        if np.isclose(abs(pt_point[1]-split_pt_line[-1][1]), t_step):
            # continuous line
            split_pt_line.append(pt_point)
        else:
            # there is a split in the line here
            # add the existing list to the list of split lines and start a new list
            split_pt_lines.append(deepcopy(split_pt_line))
            split_pt_line = [pt_point]

    # add the final split_pt_line to the superlist
    split_pt_lines.append(deepcopy(split_pt_line))

    # add all the lines to the dictionary of split lines
    for i in range(len(split_pt_lines)):
        split_phase_transition_points["{}_{}".format(index_str, i)] = split_pt_lines[i]

########################################################### plotting #####################################################################

print("Plotting phase diagram")

plt.xlabel("Pressure [GPa]")
plt.ylabel("Temperature [K]")

# define plotting colourmap
cmap = ListedColormap(["#E6AB02", "#66A61E", "#8000C4", "#7570B3", "#E7298A", "#1E90FF", "#1B9E77", "#20C2C2", "#D95F02", "#DC143C"])

# plot colour mesh of phases
plt.pcolormesh(pressure_list,temp_list_numbers,minimum_index_array, cmap=cmap, vmin=0, vmax=len(structure_list)-1, alpha=0.75)

# plot phase boundaries

# determine maximum temperature and pressure to set nice upper boundaries for the plot
# set low starting values
max_press = 0
max_temp = 0

for index_str, pt_line in split_phase_transition_points.items():

    # check if we have found new maximum values
    current_max = np.amax(pt_line, axis=0)
    current_max_press = current_max[0]
    current_max_temp = current_max[1]

    if current_max_press > max_press:
        max_press = current_max_press
    if current_max_temp > max_temp:
        max_temp = current_max_temp

    x, y = zip(*pt_line)

    plt.plot(x, y, "#000080")

# set plot parameters

# when we have given melt data, set the xrange to the limits of the melt data
if "melt_curve.dat" in ls_top:
    x_limits = [int(min(melt_pressures)), int(max(melt_pressures))]
else:
    x_limits = [0, round_to_nearest_larger_five(max_press)]
y_limits = [0, round_to_nearest_larger_five(max_temp)]

plt.xlim(x_limits[0], x_limits[1])
plt.ylim(y_limits[0], y_limits[1])

plt.savefig("phase_diagram.png", dpi=300)
plt.savefig("phase_diagram.pdf")
plt.close()


