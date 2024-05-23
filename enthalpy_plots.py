# script to plot the enthalpies in *enthalpy*.agr files (generated by cryan and renamed)
# need to ensure all files have the same structures, and are referenced to the same structure
# Structures can be named "Structure : Quality", where the script considers structures for which "Structure" is the same identical
# Quality is a differentiating characteristic, for instance "No SOC" vs "SOC", or different DFT XC functionals

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

from copy import deepcopy
from scipy import interpolate

# note that the default labels in the _enthalpy.agr files are chemical_formula-structure_name

# we convert to int, there is no good reason really to start at fractional pressures, and int makes things cleaner
low_press = int(sys.argv[1])
high_press = int(sys.argv[2])
reference_structure = sys.argv[3]

press_range = high_press - low_press

# check if user-supplied labels are present
# labels.txt format:
# structure_name_in_file label_for plot
# ...

ls = os.listdir()

labels = {}
if "labels.txt" in ls:
	with open("labels.txt", "r") as labels_file:
		for line in labels_file:
			label = line.split()
			labels[label[0]] = label[1]
else:
	# for unified access, we fill the labels dictionary anyways, even though the label is just identical to the name
	labels[reference_structure] = reference_structure

# read input data
enthalpies = {reference_structure: [np.linspace(low_press, high_press, high_press-low_press+1), np.zeros(high_press-low_press+1)]}
structure_list = [reference_structure]
quality_list = []

agr_files = [file for file in ls if "enthalpy" in file and file.endswith("agr")]
agr_files.sort()
#agr_files.sort(reverse=True)

for file in agr_files:

	infile = open(file, "r")

	for line in infile:

		# the data we want is stored in the blocks at the end of the file; the beginning of these blocks is indicated by a few header lines setting the style. This is what we detect.
		if "legend" in line and "s" in line:
			# extra check for "s" to avoid detecting the lines including "legend" in the general file preamble setting global style
			
			if reference_structure in line:
				# we have already set the data for the reference structure to 0; ignore it
				continue
			else:
				# this we read
				# extract the structure name
				structure_name_parts = ' '.join(line.split()[3:]).strip('"').split(" : ")
				structure_name = structure_name_parts[0]
				structure_quality = structure_name_parts[-1]

				if not structure_name in structure_list:
					structure_list.append(structure_name)

				if not structure_name in labels.keys():
					labels[structure_name] = structure_name

				if not structure_quality in quality_list:
					quality_list.append(structure_quality)

				# skip the next 6 lines which are headers
				for i in range(6):
					infile.readline()

				# now we get to the data
				pressure = []
				enthalpy = []

				for data_line in infile:

					if data_line.startswith("&"):
						# "&" marks the end of a data block
						if not structure_quality in enthalpies.keys():
							enthalpies[structure_quality] = {structure_name: [np.array(deepcopy(pressure)), np.array(deepcopy(enthalpy))]}
						else:
							enthalpies[structure_quality][structure_name] = [np.array(deepcopy(pressure)), np.array(deepcopy(enthalpy))]
						break

					else:
						data = data_line.split()

						current_pressure = float(data[0])

						if (current_pressure > low_press or np.isclose(current_pressure, low_press)) and (current_pressure < high_press or np.isclose(current_pressure, high_press)):

							pressure.append(current_pressure)
							enthalpy.append(1000*float(data[1])) # convert from eV/atom to meV/atom

		else:
			continue

	infile.close()

# plotting
# define colour list
colours = ["#E6AB02", "#66A61E", "#8000C4", "#E7298A", "#1E90FF", "#1B9E77", "#20C2C2", "#D95F02", "#DC143C"]

# define plotting dash types
linestyles = ["-", "--", "-.", ":"]

plt.xlabel("Pressure [GPa]")
plt.ylabel("Enthalpy [meV/formula unit]")
plt.xlim(low_press, high_press)

# determine sensible tick spacing on the x-axis
max = round(press_range/50)*5

for i in range(max, 0, -1):

	if press_range%i == 0:
		tick_step = i
		break

if max <= 1:
	tick_step = 1

plt.xticks(np.arange(low_press, high_press+1, tick_step))

# to set a sensible enthalpy window, determine the largest and lowest enthalpies
max_enthalpy = 0
min_enthalpy = 0

# plot DFT data
for structure in enthalpies.keys():

	if structure == reference_structure:
		plt.plot(enthalpies[structure][0], enthalpies[structure][1], color=colours[structure_list.index(structure)],
			     linestyle="solid", label=labels[structure])
		continue

	# if this isn't true, "structure" is actually a differentiating characteristic
	for actual_structure, pressure_enthalpy in enthalpies[structure].items():
		plt.plot(pressure_enthalpy[0], pressure_enthalpy[1], color=colours[structure_list.index(actual_structure)],
			     linestyle=linestyles[quality_list.index(structure)], label="{} - {}".format(labels[actual_structure], structure))

		if np.amax(pressure_enthalpy[1]) > max_enthalpy:
			max_enthalpy = np.amax(pressure_enthalpy[1])
		if np.amin(pressure_enthalpy[1]) < min_enthalpy:
			min_enthalpy = np.amin(pressure_enthalpy[1])

# round maximum and minimum enthalpies to the next 10 to get y range
y_max = round(max_enthalpy, -1)
y_min = round(min_enthalpy, -1)

# make sure the entire range is included
if y_max < max_enthalpy:
	y_max += 10
if y_min > min_enthalpy:
	y_min -= 10

plt.ylim(y_min, y_max)

plt.legend()

plt.savefig("enthalpies.png", dpi=300)
#plt.savefig("dft_eddp_dev_enthalpies.pdf")
plt.close()



































