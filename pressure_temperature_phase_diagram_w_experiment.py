# script to generate a pressure-temperature phase diagram, including points for experimental data
# requires four input files:
# phase_diagram_data.dat, which has the same format as the output of / will likely have been generated by my pressure_temperature_phase_diagram.py script
# phase_transition_points.dat, likewise generated by the pressure_temperature_phase_diagram.py script
# structure_indices.dat, like the two preceding scripts generated by the pressure_temperature_phase_diagram.py script
# phase_diagram_data_experiment.dat, which has the same format as phase_diagram_data.dat, but will contain much fewer points and be usually generated by hand

# written by Pascal Salzbrenner, pts28@cam.ac.uk
# with thanks to Se Hun Joo for the code using Pandas to read the phase diagram data and for introducing me to pcolormesh

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from copy import deepcopy
from matplotlib.colors import ListedColormap

################################################################ hardcode some input variables ################################################################

# temperature increment of interpolated input data
t_step = 5

# define colour for phase boundaries and experimental symbols
pb_colour="#000080"

# define plotting colourmap
cmap = ListedColormap(["#E6AB02", "#66A61E", "#8000C4", "#7570B3", "#E7298A", "#1E90FF", "#1B9E77", "#20C2C2", "#D95F02", "#DC143C"])

# define symbols for experimental data
symbols = ["x", "o", "v", "^", "<", ">", "s", "+", "H", "D"]


################################################################### reading of input files ###################################################################

# read structures_indices.dat file
structure_list = []

with open("structures_indices.dat", "r") as index_file:
    for line in index_file:
        structure_list.append(line.split()[1])

# read phase_diagram_data.dat file and put it into plottable shape
df_sim = pd.read_csv("phase_diagram_data.dat", sep=", ", engine="python")

colx = "Pressure [GPa]"
coly = "Temperature [K]"
colz = "Ground state structure"

df_sim = df_sim.pivot(columns=colx, index=coly, values=colz)

X = df_sim.columns.values
Y = df_sim.index.values

# convert structure names into indices
minimum_index_list = []

for row in df_sim.values:
	
    # the data is already organised into a grid, so we go row by row
    single_row = []

    for structure in row:

    	if structure.lower() == "liquid":
    		single_row.append(len(structure_list))
    	else:
    		single_row.append(structure_list.index(structure))

    minimum_index_list.append(deepcopy(single_row))

# set up mesh of pressure, temperature data
pressure_list,temp_list=np.meshgrid(X,Y)

# read phase_transition_points.dat file
pt_file = open("phase_transition_points.dat", "r")

# define dictionary in which to keep unprocessed input data
pt_lines = {}

# read past first two lines, which are comments
for i in range(2):
	pt_file.readline()

# read data
for line in pt_file:
	data = line.split()
	pt_index = data[0]
	pt_pressure = float(data[1])
	pt_temp = float(data[2])

	if pt_index in pt_lines.keys():
		pt_lines[pt_index].append([pt_pressure, pt_temp])
	else:
		pt_lines[pt_index] = [[pt_pressure, pt_temp]]

pt_file.close()

# treatment for extremely obscure problem I faced with lead
# very near the phase transition for a tiny bit it's first Fm-3m-liquid, then P63/mmc-liquid, then Fm-3m-liquid
# this leads to regions where visually, there are two melt boundaries - to avoid this, we split lists where successive points
# are separated by more than t_step (as the outer loop is over temperature, this is what implies a discontinuous line)

split_phase_transition_points = []

for pt_line in pt_lines.values():

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
        split_phase_transition_points.append(split_pt_lines[i])

# read experimental_phase_diagram_data.dat and process data
df_exp = pd.read_csv("phase_diagram_data_experiment.dat", sep=", ", engine="python")

data_labels = list(df_exp[colz].unique())

experimental_data = {}

for i in range(len(data_labels)):

	data_type = data_labels[i]

	pressure_data = df_exp.loc[df_exp[colz] == data_type, colx].to_list()
	temperature_data = 	df_exp.loc[df_exp[colz] == data_type, coly].to_list()

	experimental_data[data_type] = [pressure_data, temperature_data, symbols[i]]

########################################################################## plotting ##########################################################################

plt.xlabel("Pressure [GPa]")
plt.ylabel("Temperature [K]")

# plot colour mesh of phases
plt.pcolormesh(pressure_list,temp_list,minimum_index_list, cmap=cmap, vmin=0, vmax=len(structure_list)-1, alpha=0.75)

# plot phase boundaries
# determine maximum temperature and pressure to set nice upper boundaries for the plot
# set low starting values
max_press = 0
max_temp = 0

for pt_line in split_phase_transition_points:

    # check if we have found new maximum values
    current_max = np.amax(pt_line, axis=0)
    current_max_press = current_max[0]
    current_max_temp = current_max[1]

    if current_max_press > max_press:
        max_press = current_max_press
    if current_max_temp > max_temp:
        max_temp = current_max_temp

    x, y = zip(*pt_line)

    plt.plot(x, y, pb_colour)

# plot experimental data
for data_type, data in experimental_data:
	plt.scatter(data[0], data[1], c=pb_colour, marker=data[2], label=data_type)

plt.legend()

# set plot parameters
# when we have given melt data, set the xrange to the limits of the melt data
if "melt_curve.dat" in ls_top:
    x_limits = [int(min(melt_pressures)), int(max(melt_pressures))]
else:
    x_limits = [0, round_to_nearest_larger_five(max_press)]
y_limits = [0, round_to_nearest_larger_five(max_temp)]

plt.xlim(x_limits[0], x_limits[1])
plt.ylim(y_limits[0], y_limits[1])

plt.savefig("phase_diagram_w_experiment.png", dpi=300)
#plt.savefig("phase_diagram_w_experiment.pdf")
plt.close()
