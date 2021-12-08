# script to calculate the corrections in the hopping elements arising due to the displacement of atoms by Caesar compared to the static lattice
# ie for each hopping element, we subtract the static lattice value
# run in directory containing the configurations directory, or in the configurations directory

# input files
# static wannier90_hr.dat
# Caesar configuration wannier90_hr.dat

# output files
# correction_hr.dat - the corrections at each hopping term, written in the same directory as the original wannier90_hr.dat

# written by Pascal Salzbrenner, pts28@cam.ac.uk

import os
import sys
import numpy as np

########################################################## define helper function ##########################################################

def to_hr_list(hopping_dict):
    """Function to map a hopping dict to a list, which can then be written to an _hr.dat file via "\n".join(hopping_list)
    :param dict hopping_dict: dict of the form {vector_tuple : hopping_matrix (size num_wann by num_wann), ...}

    :param returns hopping_list: list of strings giving one line each in an _hr.dat file"""

    # set up hopping_list
    hopping_list = []

    # iterate over elements of hopping_dict
    for hop_vector, hop_matrix in hopping_dict.items():
        # access all elements of the hopping matrix
        for i in range(len(hop_matrix)):
            for j in range(len(hop_matrix)):
                hopping_list.append("{}    {}    {}    {}    {}    {}    {}".format(hop_vector[0], hop_vector[1], hop_vector[2], j+1, i+1,
                                                                                    hop_matrix[j][i].real, hop_matrix[j][i].imag))
                # entries in matrix numbered 0-len(hop_matrix)-1, but numbered 1-len(hop_matrix) in _hr.dat file

    return hopping_list

############################################################### parse input ###############################################################

# point under consideration passed as first input argument
point = sys.argv[1]

# this list will contain all the lines which will be written to correction_hr.dat
correction_lines = []
correction_lines.append("Electron-phonon coupling correction terms as the difference between the static lattice and the configuration")

# check if we are in the configurations directory by checking if there is a configurations directory in the directory we are in
# if there is no configurations directory, it is assumed we are in it
if not "configurations" in os.listdir():
    path = "./"
else:
    path  = "configurations/"

############################################### read data from static wannier90_hr.dat file ###############################################

# open static input file
static_hr = open(path+"static/wannier90_files/wannier90_hr.dat", "r")

# set up dictionary to contain static hopping elements
static_hopping_dict = {}

# skip comment line
static_hr.readline()

# next line contains the number of wannier orbitals per unit cell, the one after that the number of hopping vectors
num_wann = int(static_hr.readline())
correction_lines.append(str(num_wann))
num_vec = int(static_hr.readline())
correction_lines.append(str(num_vec))

# based on the number of vectors, determine how many lines to read past, using the fact that the file has 15 points per line
read_past = num_vec//15
remaining_elements = num_vec%15

# write the weight lines for the output file, setting all weights to 1
for i in range(read_past):
    correction_lines.append("    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1")
    static_hr.readline()

if remaining_elements:
    # remaining_elements != 0 => there is a line that is not completely filled - for writing: some more weights must be written
    additional_line = "    1"
    for i in range(remaining_elements-1):
        additional_line += "    1"

    correction_lines.append(additional_line)

    # for reading: one additional line with < 15 elements must be read past if there is a not fully filled line
    static_hr.readline()

for line in static_hr:

    # read the hopping_vector, the from_orbital, the to_orbital, and the real and imaginary components of the strength
    hopping_data = line.split()
    hopping_vector = tuple(hopping_data[0:3])
    from_orbital = int(hopping_data[3])-1 # renumber to conform to internal Python representation
    to_orbital = int(hopping_data[4])-1
    strength = complex(float(hopping_data[5]), float(hopping_data[6]))

    # set up hopping matrix for this vector if the entry doesn't exist yet
    if hopping_vector not in static_hopping_dict.keys():
        static_hopping_dict[hopping_vector] = np.zeros((num_wann, num_wann), dtype=np.cdouble)

    static_hopping_dict[hopping_vector][from_orbital][to_orbital] = strength

static_hr.close()

############################################ read data from configuration wannier90_hr.dat file ############################################

# open configuration input file
configuration_hr = open(path+"point.{}/".format(point)+"wannier90_files/wannier90_hr.dat", "r")

# set up dictionary to contain configuration hopping elements
configuration_hopping_dict = {}

# skip comment line
configuration_hr.readline()

# next line contains the number of wannier orbitals per unit cell, the one after that the number of hopping vectors
num_wann_configuration = int(configuration_hr.readline())
num_vec_configuration = int(configuration_hr.readline())

# if the static calculation and the configuration are part of the same set of Caesar calculations, these values should agree
if num_wann_configuration != num_wann:
    raise ValueError("The static tight-binding model is constructed on a basis of {} WFs, whereas point {}'s model is based on {} WFs. These will be very different models which should not be treated using this script".format(
    num_wann, point, num_wann_configuration))
if num_vec_configuration != num_vec:
    raise ValueError("The static tight-binding model contains {} different lattice vectors, whereas point {}'s model contains {}. This is unexpected. Are you sure you are comparing two compatible models?".format(
    num_vec, point, num_vec_configuration))

for i in range(read_past):
    configuration_hr.readline()

for line in configuration_hr:

    # read the hopping_vector, the from_orbital, the to_orbital, and the real and imaginary components of the strength
    hopping_data = line.split()
    hopping_vector = tuple(hopping_data[0:3])
    from_orbital = int(hopping_data[3])-1 # renumber to conform to internal Python representation
    to_orbital = int(hopping_data[4])-1
    strength = complex(float(hopping_data[5]), float(hopping_data[6]))

    # set up hopping matrix for this vector if the entry doesn't exist yet
    if hopping_vector not in configuration_hopping_dict.keys():
        configuration_hopping_dict[hopping_vector] = np.zeros((num_wann, num_wann), dtype=np.cdouble)

    configuration_hopping_dict[hopping_vector][from_orbital][to_orbital] = strength

configuration_hr.close()

##################################### calculate corrections and write them to corrections_hr.dat file #####################################

# open output file
correction_hr = open(path+"point.{}/".format(point)+"wannier90_files/correction_hr.dat", "w")

# open warning file for potential warnings relating to the absence of certain hopping vectors in the specific configuration
warning_file = open(path+"point.{}/".format(point)+"/wannier90_files/correction_calculation_warning.txt", "w")

# write out preample
correction_hr.write("\n".join(correction_lines))
correction_hr.write("\n")

for hop_vector, hop_matrix in static_hopping_dict.items():
    if hop_vector in configuration_hopping_dict.keys():
        correction_lines = to_hr_list({hop_vector: configuration_hopping_dict[hop_vector]-hop_matrix})
        correction_hr.write("\n".join(correction_lines))
        correction_hr.write("\n")
    else:
        warning_file.write("Warning: static hopping vector ({}, {}, {}) not in this configuration's hoppings. Skipping.".format(hop_vector[0], hop_vector[1], hop_vector[2]))
        continue

correction_hr.close()
warning_file.close()
