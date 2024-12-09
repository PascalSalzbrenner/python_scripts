# script to plot band structures from CASTEP output
# will plot all seed-*.bands structures contained in the directory where the script is run
# they must all be calculated along the same path, the density will automatically be scaled
# the high-symmetry path must be supplied in the hsp.dat file - format is:
# kx ky kz name
# can ideally just be copied from a .cell file used to initialise one of the calculations

# written by Pascal Salzbrenner, pts28@cam.ac.uk

import os
import re
import sys
import numpy as np
import matplotlib.pyplot as plt

### define a few constants ###

# CASTEP .bands files use Hartree for the energy - multiply by constant below to convert to eV
hartree_to_eV = 27.211407953
# this is divided by the number of k-points between high-symmetry points to determine the distance between points
adaptive_distance = 1

# plot y-axis minima and maxima
y_min = float(sys.argv[1])
y_max = float(sys.argv[2])

### function definitions ###

def read_hsp():
    """Function to read hsp.dat
    :returns list[list[str]] hsp: the list of high-symmetry points as a string, as they are mostly used for parsing"""

    hsp = []

    with open('hsp.dat') as f:
        for line in f:
            hsp.append(line.split())

    return hsp

def locate_hsps(filename, hsp_list):
    """Function to locate high-symmetry points in a .bands file using the information from hsp.dat
    :param str filename: name of the file in which we search the hsp data
    :param list[list[str]] hsp_list: list of high-symmetry points produced by read_hsp

    :returns list[int] k_points: the numbers of the k-points"""

    k_points = []

    # read contents of file
    with open(filename) as f:
        file_text = f.read()

    for k_point in hsp_list:
        pattern = re.compile(r"K-point\s+[0-9]+\s+{}\.*0+\s+{}\.*0+\s+{}\.*0+\s+".format(k_point[0], k_point[1], k_point[2]))

        for match in pattern.findall(file_text):
            index = int(match.split()[1])
            if index not in k_points:
                k_points.append(index)

    # sort because if the path contains the same point several times, we will detect all the matches at once
    return sorted(k_points)

### main code ###

# initialise dictionary to hold x and y data for each .bands file
# format {file_identifier: [x_data], [y_data]}
data_dict = {}

# get high-symmetry points
hsp = read_hsp()

bands_files = sorted([file for file in os.listdir() if file.endswith(".bands")])

for file in bands_files:

    # initialisation
    file_kpoints = locate_hsps(file, hsp)
    bands_file = open(file, "r")

    # lists for x and y data
    x_data = []
    y_data = []

    # create and fill a directory for the x-position corresponding to a given index
    # this is necessary because the k-points in CASTEP .bands files are not necessarily ordered
    x_positions = [0]

    for i in range(file_kpoints[-1]-1):
        # update displacement when i reaches a high-symmetry point

        if i+1 in file_kpoints:
            # next high-symmetry point reached
            hsp_index = file_kpoints.index(i+1)

            x_disp = (adaptive_distance * np.linalg.norm(np.array(hsp[hsp_index + 1][:3], dtype=float) -
                                                         np.array(hsp[hsp_index][:3], dtype=float)) /
                      (file_kpoints[hsp_index + 1] - file_kpoints[hsp_index]))

        # incorrect position update for last k-point but that is immaterial, as we will never reach here again
        x_positions.append(x_positions[-1] + x_disp)

    # fill y_data with empty lists to contain the data for the correct points
    for i in range(file_kpoints[-1]):
        y_data.append([])

    # skip past header
    for i in range(9):
        header_line = bands_file.readline()

        # read Fermi energy
        if header_line.startswith("Fermi"):
            fermi_energy = float(header_line.split()[-1])

    for line in bands_file:

        if line.startswith("K-point"):

            k_point_index = int(line.split()[1])

            x_data.append(x_positions[k_point_index-1])

        elif line.startswith("Spin"):
            continue

        else:
            # actual k-point
            y_data[k_point_index-1].append((float(line.split()[0])-fermi_energy) * hartree_to_eV)

    bands_file.close()

    file_identifier = "-".join(file.split("-")[1:]).rstrip(".bands")

    # y-data is automatically sorted due to the way we've set it up
    # x-data must be sorted to correspond to this
    # transpose the y_data to put it in a format easy to plot with pyplot
    data_dict[file_identifier] = sorted(x_data).copy(), np.array(y_data[:]).transpose()

### plotting ###

# define colour palette
colour_palette = ["#E6AB02", "#66A61E", "#8000C4", "#7570B3", "#E7298A", "#1E90FF", "#1B9E77", "#20C2C2", "#D95F02", "#DC143C"]

for key, value in data_dict.items():

    for i in range(len(value[1])):
        if i == 0:
            label=key
            colour = colour_palette.pop(0)
        else:
            label = None

        plt.plot(value[0], value[1][i], label=label, linestyle='-', color=colour)

plt.ylabel("Energy [eV]")

# draw in the Fermi energy
plt.axhline(y=0, color='black', linestyle='--')

# high symmetry points
xtick_pos = []
xtick_labels = []

for point in file_kpoints:

    label = hsp[file_kpoints.index(point)][3]

    if label == "G":
        label = "Γ"

    xtick_labels.append(label)

    x_pos = x_positions[point-1]
    xtick_pos.append(x_pos)

    plt.axvline(x=x_pos, color='black', linestyle='-')

plt.xticks(xtick_pos, xtick_labels)
plt.tick_params(axis='x', direction='in')

plt.xlim(x_positions[0], x_positions[-1])

plt.legend()
plt.ylim(y_min, y_max)
plt.savefig("band_structure.png", dpi=300)
