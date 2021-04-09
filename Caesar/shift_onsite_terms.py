# script to rigidly shift all onsite terms of a wannier90_hr.dat tight-binding model
# this shifts every eigenvalue (and hence the Fermi level) by the amount of the shift, but leaves the dispersion unaffected
# subtract the energy of R = (0,0,0), from_orbital=to_orbital=1, from every onsite term

# input file
# Caesar configuration wannier90_hr.dat

# output file
# wannier90_shifted_hr.dat - the corrections at each hopping term, written in the same directory as the original wannier90_hr.dat

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
                hopping_list.append("{:>5}{:>5}{:>5}{:>5}{:>5}{:>12.6f}{:>12.6f}".format(hop_vector[0], hop_vector[1], hop_vector[2], j+1,
                                                                                         i+1, hop_matrix[j][i].real, hop_matrix[j][i].imag))
                # entries in matrix numbered 0-len(hop_matrix)-1, but numbered 1-len(hop_matrix) in _hr.dat file

    return hopping_list

############################################################### parse input ###############################################################

# point under consideration passed as first input argument
point = sys.argv[1]

# this list will contain all the lines which will be written to wannier90_shifted_hr.dat
shifted_lines = []

# check if we are in the configurations directory by checking if there is a configurations directory in the directory we are in
# if there is no configurations directory, it is assumed we are in it
if not "configurations" in os.listdir():
    path = "./"
else:
    path  = "configurations/"

######################### read data from configuration wannier90_hr.dat file and write to wannier90_shifted_hr.dat #########################

# open configuration input file
configuration_hr = open(path+"point.{}/".format(point)+"wannier90_files/wannier90_hr.dat", "r")

# set up dictionary to contain configuration hopping elements
configuration_hopping_dict = {}

# skip comment line
configuration_hr.readline()

# read the number of Wannier functions and number of vectors
num_wann = int(configuration_hr.readline())
num_vec = int(configuration_hr.readline())
shifted_lines.append(str(num_wann))
shifted_lines.append(str(num_vec))

# based on the number of vectors, determine how many lines to read past, using the fact that the file has 15 points per line
read_past = num_vec//15
remaining_elements = num_vec%15

# write the weight lines for the output file, setting all weights to 1
for i in range(read_past):
    shifted_lines.append(configuration_hr.readline().rstrip("\n"))

if remaining_elements:
    # remaining_elements != 0 => there is a line that is not completely filled - for writing: some more weights must be read and written

    shifted_lines.append(configuration_hr.readline().rstrip("\n"))

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

####################################### to wannier90_shifted_hr.dat file, with shifted onsite terms #######################################

# open output file
shifted_hr = open(path+"point.{}/".format(point)+"wannier90_files/wannier90_shifted_hr.dat", "w")

# add comment line to beginning list and write preamble to output files
shifted_lines = ["wannier90 TB model with onsite terms shifted by {} eV".format(configuration_hopping_dict[("0", "0", "0")][0][0])] + \
                                                                         shifted_lines
shifted_hr.write("\n".join(shifted_lines))
shifted_hr.write("\n")

for hop_vector, hop_matrix in configuration_hopping_dict.items():
    if hop_vector == ("0", "0", "0"):
        # shift the onsite terms
        shifted_lines = to_hr_list({hop_vector: hop_matrix-hop_matrix[0][0]*np.identity(num_wann)})
    else:
        shifted_lines = to_hr_list({hop_vector: hop_matrix})

    shifted_hr.write("\n".join(shifted_lines))
    shifted_hr.write("\n")

shifted_hr.close()
