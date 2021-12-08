# script to add the averaged corrections to the calculation for the static lattice in the primitive unit cell
# run in directory containing the two input files

# input files
# primitive static wannier90_hr.dat
# averaged configuration_correction_average_hr.dat

# output files
# wannier90_renormalised_hr.dat - the temperature-renormalised hopping parameters

# written by Pascal Salzbrenner, pts28@cam.ac.uk

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

########################################################### set up output lines ###########################################################

# this list will contain all the lines which will be written to wannier90_corrected_hr.dat
renorm_lines = []
renorm_lines.append("Primitive static lattice, with renormalisation terms due to electron-phonon coupling added")

############################################### read data from static wannier90_hr.dat file ###############################################

# open static input file
static_hr = open("wannier90_hr.dat", "r")

# set up dictionary to contain static hopping elements
static_hopping_dict = {}

# skip comment line
static_hr.readline()

# next line contains the number of wannier orbitals per unit cell, the one after that the number of hopping vectors
num_wann = int(static_hr.readline())
renorm_lines.append(str(num_wann))
num_vec = int(static_hr.readline())

# based on the number of vectors, determine how many lines to read past, using the fact that the file has 15 points per line
read_past = num_vec//15
remaining_elements = num_vec%15

if remaining_elements:
    #one additional line with < 15 elements must be read past if there is a not fully filled line
    read_past+=1

for i in range(read_past):
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

####################################### read data from configuration_correction_average_hr.dat file #######################################

# open configuration input file
correction_hr = open("configuration_correction_average_hr.dat", "r")

# set up dictionary to contain configuration hopping elements
correction_hopping_dict = {}

# skip comment line
correction_hr.readline()

# next line contains the number of wannier orbitals per unit cell, the one after that the number of hopping vectors
num_wann_correction = int(correction_hr.readline())
num_vec_correction = int(correction_hr.readline())

# if the static calculation and the configuration are part of the same set of Caesar calculations, these values should agree
if num_wann_correction != num_wann:
    raise ValueError("The primitive tight-binding model is constructed on a basis of {} WFs, whereas the corrections are based on {} WFs. This will not work".format(
    num_wann, point, num_wann_configuration))

# based on the number of vectors, determine how many lines to read past, using the fact that the file has 15 points per line
read_past = num_vec_correction//15
remaining_elements = num_vec_correction%15

if remaining_elements:
    #one additional line with < 15 elements must be read past if there is a not fully filled line
    read_past+=1

for i in range(read_past):
    correction_hr.readline()

for line in correction_hr:

    # read the hopping_vector, the from_orbital, the to_orbital, and the real and imaginary components of the strength
    hopping_data = line.split()
    hopping_vector = tuple(hopping_data[0:3])
    from_orbital = int(hopping_data[3])-1 # renumber to conform to internal Python representation
    to_orbital = int(hopping_data[4])-1
    strength = complex(float(hopping_data[5]), float(hopping_data[6]))

    # set up hopping matrix for this vector if the entry doesn't exist yet
    if hopping_vector not in correction_hopping_dict.keys():
        correction_hopping_dict[hopping_vector] = np.zeros((num_wann, num_wann), dtype=np.cdouble)

    correction_hopping_dict[hopping_vector][from_orbital][to_orbital] = strength

correction_hr.close()

######################## calculate renormalised hopping terms and write them to wannier90_renormalised_hr.dat file ########################

# iterate over elements in the corrections dict
# if an element isn't present in the static calculation, the correction is taken as the entire value, and added as an entry to static dict
for hop_vector, hop_matrix in correction_hopping_dict.items():
    if hop_vector in static_hopping_dict.keys():
        static_hopping_dict[hop_vector] += hop_matrix
    else:
        static_hopping_dict[hop_vector] = hop_matrix

num_vec_out = len(static_hopping_dict)
renorm_lines.append(str(num_vec_out))

# determine how many weight lines to write
full_lines = num_vec_out//15
remaining_elements = num_vec_out%15

# write the weight lines for the output file, setting all weights to 1
for i in range(full_lines):
    renorm_lines.append("    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1")

if remaining_elements:
    # remaining_elements != 0 => there is a line that is not completely filled - some more weights must be written
    additional_line = "    1"
    for i in range(remaining_elements-1):
        additional_line += "    1"

    renorm_lines.append(additional_line)

# open output file
renorm_hr = open("wannier90_renormalised_hr.dat", "w")

# write out preample
renorm_hr.write("\n".join(renorm_lines))
renorm_hr.write("\n")

for hop_vector, hop_matrix in static_hopping_dict.items():
    renorm_lines = to_hr_list({hop_vector: hop_matrix})
    renorm_hr.write("\n".join(renorm_lines))
    renorm_hr.write("\n")

renorm_hr.close()
