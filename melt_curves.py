# script to plot the melt curves in enthalpy*.dat files, including their standard deviations
# format: P Tm std_dev
# first line contains "# label"

import os
import numpy as np
import matplotlib.pyplot as plt

ls = os.listdir()
melt_files = [file for file in ls if "melt_curve" in file and file.endswith("dat")]
melt_files.sort()
#agr_files.sort(reverse=True)

label_list = []
melting_points = {}

for file in melt_files:

	infile = open(file, "r")

	label = " ".join(infile.readline().split()[1:])

	label_list.append(label)
	melting_points[label] = [[], [], []]

	for line in infile:

		if line.startswith("#"):
			continue

		data=line.split()

		# three lists, one each for pressure, melting point, standard deviation
		melting_points[label][0].append(float(data[0]))
		melting_points[label][1].append(float(data[1]))
		melting_points[label][2].append(float(data[2]))

	infile.close()

# read experimental data
if "experimental_data.dat" in ls:
	exp_data = open("experimental_data.dat", "r")

	exp_pressures = []
	exp_melting_points = []
	exp_errors = []

	for line in exp_data:
	    
	        if line.startswith("#"):
	            continue
	    
	        data = line.split()
	    
	        exp_pressures.append(float(data[0]))
	        exp_melting_points.append(float(data[1]))
	        exp_errors.append(float(data[2]))

# plotting
# define colour list
colours = ["#E6AB02", "#66A61E", "#8000C4", "#E7298A", "#1E90FF", "#1B9E77", "#20C2C2", "#D95F02", "#DC143C"]

plt.xlabel("Pressure [GPa]")
plt.ylabel(r'$\mathrm{T_m}$ [K]')

# just use the last dataset we iterated over
low_press=melting_points[label][0][0]
high_press=melting_points[label][0][-1]
press_range=high_press-low_press
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
plt.minorticks_on()

for label, data in melting_points.items():
	plt.plot(data[0], data[1], color=colours[label_list.index(label)], linestyle="solid", label=label)
	plt.fill_between(data[0], np.array(data[1])+np.array(data[2]), np.array(data[1])-np.array(data[2]), color=colours[label_list.index(label)], alpha=0.5)

if "experimental_data.dat" in ls:
	plt.errorbar(exp_pressures, exp_melting_points, yerr=exp_errors, fmt="x", color="#000080", label="Experimental data")

plt.legend()

plt.savefig("melt_curves.png", dpi=300)
plt.close()



































