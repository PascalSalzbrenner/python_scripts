# script to average the wannier90 hopping parameter corrections due to electron-phonon coupling calculated for different on average
# equivalent positions of a supercell generated using the in-house Caesar code
# run in directory containing the configurations directory and the supercell.dat file

# input files
# correction_hr.dat for the configuration in question; wannier90.win for static
# supercell.dat file from Caesar
# orbitals_per_atom.dat file in the form:
# atom_type number_of_orbitals
# ... ...
# ... ...
# theoretically, this could be read out from the begin_projections block of the wannier90.win file, but there are a lot of cases to handle
# can be implemented in the future

# output
# primitive_cell_corrections_hr.dat: hopping elements averaged and mapped to the primitive orbitals
# primitive_cell.win: contains the same information as the static wannier90.win, mapped to the primitive cell
# primitive_cell_centres.xyz: contains the WCCs in the primitive cell, enforced to the positions of the atoms

# written by Pascal Salzbrenner, pts28

import sys
import numpy as np
from copy import deepcopy

########################## define functions to convert vectors quickly between the primitive and supercell bases ##########################

def get_expansion_coefficients(old_basis, new_basis):
    """Function to determine the coefficients need to expand the old_basis in the new_basis
    :param np.ndarray old_basis: an N by N set of vectors, to be expanded in new_basis
    :param np.ndarray new_basis: an N by N set of vectors in which we expand old_basis - must be invertible

    :returns np.ndarray expansions_coefficients: an N by N set of expansion coefficients of old_basis in new_basis.
                                                 The first row gives the expansion coefficients of first old_basis
                                                 vector in new_basis, that is, of the first row of old_basis. Etc
                                                 Note that for supercell construction, this is the inverse of the
                                                 supercell matrix"""

    return np.dot(old_basis, np.linalg.inv(new_basis))

def change_basis(vector, expansion_coefficients):
    """Function to convert the representation of the vector from an old basis to a new basis, where the two bases are
    related by the expansion_coefficient matrix
    :param np.ndarray vector: length N array which is to be converted
    :param np.ndarray expansion_coefficients: NxN array giving expansion coefficients of the old in the new basis

    :returns np.ndarray vector_in_new_basis: vector represented in the new basis"""

    return np.dot(expansion_coefficients.T, vector)

def cartesian_to_lattice_basis(point, lattice):
    """Function to convert a point given in Cartesian coordinates into coordinates in the lattice
    :param np.ndarray k_point: numpy array of length dimensions indicating a point in Cartesian coordinates
    :param np.ndarray lattice: numpy array of size dimensions x dimensions, containing the lattice vectors as rows

    :returns np.ndarray point_lattice_basis: numpy array of length dimensions giving the point in the lattice basis
    """

    return np.dot(np.linalg.inv(lattice.T), point)

######################################################### read out supercell_size #########################################################

# point under consideration passed as first input argument
point = sys.argv[1]

supercell_size = []

# works this way for diagonal supercells
with open("supercell.dat", "r") as supercell:

    for i in range(3):
        line = supercell.readline()
        supercell_size.append(int(line.split()[i]))

# determine number of copies of the primitive cell as the product of the 3 sizes in supercell_size
num_copies = np.product(supercell_size)

# each element of this list will be one line to write to the primitive_cell_corrections_hr.dat file
hr_lines = []

# comment line
hr_lines.append("constructed by averaging _hr.dat corrections for a {}x{}x{} supercell".format(supercell_size[0],supercell_size[1],
                                                                                                       supercell_size[2]))

##################################################### read number of orbitals per atom #####################################################

orbitals_per_atom = {}

with open("orbitals_per_atom.dat", "r") as orbital_data_file:

    for line in orbital_data_file:
        orbitals_per_atom[line.split()[0].strip()] = int(line.split()[1])

######################################################### read and write .win data #########################################################

# NB: kpoints from wannier90.win will not be converted in any way. This would lead to less dense sampling, but this block is not used in the
# tools described in the introduction anyways. It is included here for consistency's sake

win_in = open("configurations/static/wannier90_files/wannier90.win", "r")

# define variable to check if the model uses spinors (in which case half the orbitals on an atom are spin-up, and half spin-down)
# False by default
spinors = False

# each element of this list will be one line to write to the primitive_cell.win file
win_lines = []

for line in win_in:

    if line.lstrip().lower().startswith("num_bands") or line.lower().startswith("num_wann"):

        # determine the separator used - can be ":", "=", or " "
        if "=" in line:
            separator = "="
            elements = line.split("=")
        elif ":" in line:
            separator = ":"
            elements = line.split(":")
        else:
            # no exception handling as the user would not have gotten this far with a faulty wannier90.win file - " " assumed as separator
            separator = " "
            elements = line.split()

        # store num_wann for later use
        if elements[0].lower().startswith("num_wann"):
            num_wann = int(elements[1])

        # divide by the number of copies to get the number for the primitive cell
        win_lines.append("{}{}{}".format(elements[0], separator, int(elements[1])//num_copies))

    elif line.lstrip().lower().startswith("spinors"):

        # determine the separator used - can be ":", "=", or " "
        if "=" in line:
            separator = "="
            elements = line.split("=")
        elif ":" in line:
            separator = ":"
            elements = line.split(":")
        else:
            # no exception handling as the user would not have gotten this far with a faulty wannier90.win file - " " assumed as separator
            separator = " "
            elements = line.split()

        # detect if spinors is true
        if "t" in elements[1].lower():
            spinors = True

        win_lines.append(line.rstrip("\n"))

    elif line.lstrip().lower().startswith("begin unit_cell_cart"):
        win_lines.append(line.rstrip("\n"))

        # define basis matrices to obtain expansion_coefficients
        supercell_basis = []
        primitive_basis = []

        # check if units are given explicitly
        supercell_vector_line = win_in.readline()

        if supercell_vector_line.lstrip().lower().startswith("ang"):
            win_lines.append("Angstrom")
            supercell_vector_line = win_in.readline()
        elif supercell_vector_line.lstrip().lower().startswith("bohr"):
            win_lines.append("Bohr")
            supercell_vector_line = win_in.readline()
        # the line at this point is guaranteed to be the line corresponding to the first vector in the supercell

        # convert to the primitive lattice basis
        for i in range(3):
            supercell_vector = np.array(supercell_vector_line.split(), dtype=np.double)
            supercell_basis.append(supercell_vector)

            primitive_vector = supercell_vector/supercell_size[i]
            primitive_basis.append(primitive_vector)

            win_lines.append("{}  {}  {}".format(primitive_vector[0], primitive_vector[1], primitive_vector[2]))

            supercell_vector_line = win_in.readline()

        supercell_basis = np.array(supercell_basis)
        primitive_basis = np.array(primitive_basis)

        expansion_coefficients = get_expansion_coefficients(supercell_basis, primitive_basis)

        # current line is that ending unit_cell_cart
        win_lines.append(supercell_vector_line.rstrip("\n"))

    elif line.lstrip().lower().startswith("begin atoms_cart"):
        win_lines.append(line.rstrip("\n"))

        # define containers for atom positions in Cartesian coordinates, positions in the primitive_cell basis
        atom_pos_cartesian = {}
        atom_pos_primitive_cell = {}

        # define container for hopping vectors which are between cells in the primitive cell and within cells in the supercell
        hopping_vectors = [np.array([0,0,0])]
        # also define a container for the corresponding tuples for easier checking
        hopping_tuples = []

        # check if units are given explicitly
        atom_line = win_in.readline()

        if atom_line.lstrip().lower().startswith("ang"):
            win_lines.append("Angstrom")
        elif atom_line.lstrip().lower().startswith("bohr"):
            win_lines.append("Bohr")
        else:
            # units not given explicitly - the line contains the first atom
            atom_data = atom_line.split()
            atom_pos_primitive_cell[atom_data[0]] = []
            atom_pos_cartesian[atom_data[0]] = [np.array(atom_data[1:], dtype=np.double)]

        # read atom lines until we reach "end atoms_cart"
        for atom_line in win_in:
            if atom_line.lstrip().lower().startswith("end atoms_cart"):
                break
            else:
                atom_data = atom_line.split()
                if not atom_data[0] in atom_pos_primitive_cell.keys():
                    atom_pos_primitive_cell[atom_data[0]] = []
                    atom_pos_cartesian[atom_data[0]] = [np.array(atom_data[1:], dtype=np.double)]
                else:
                    atom_pos_cartesian[atom_data[0]].append(np.array(atom_data[1:], dtype=np.double))

        # determine which of the positions are in the same primitive cell, and which in different primitive cells
        # from the latter, we can obtain the additional hopping_vectors

        for atom, positions in atom_pos_cartesian.items():
            for position in positions:
                temp_position = np.copy(position)
                primitive_position = cartesian_to_lattice_basis(temp_position, primitive_basis)

                if (primitive_position < 1).all():
                    # dealing with an atom in the primitive uc
                    # in Cartesian coordinates, the position is the same
                    atom_pos_primitive_cell[atom].append(temp_position)
                else:
                    # dealing with an atom in a different primitive uc
                    # if we floor this position, we get a hopping vector to another primitive uc, which must be treated explicitly when
                    # going from supercell to primitive cell
                    # these hopping vectors are in the same order as the WCCs in wannier90_centres.xyz and hence as the orbitals
                    hopping_vector = np.floor(primitive_position)

                    # append each vector only once
                    if not tuple(hopping_vector) in hopping_tuples:
                        hopping_vectors.append(hopping_vector)
                        hopping_tuples.append(tuple(hopping_vector))

        # add atoms to win_lines for writing
        for atom, positions in atom_pos_primitive_cell.items():
            for position in positions:
                win_lines.append("{}  {}  {}  {}".format(atom, position[0], position[1], position[2]))

        # write line ending atoms_cart
        win_lines.append(atom_line.rstrip("\n"))

    else:
        win_lines.append(line.rstrip("\n"))

win_in.close()

with open("configurations/point.{}/wannier90_files/primitive_cell.win".format(point), "w") as win_out:
     win_out.write("\n".join(win_lines))

################################################## write primitive_cell_centres.xyz file ##################################################

# each element of these lists will be one line to write to the respective _centres.xyz file
xyz_atomic_lines = []

# the first line of the _centres.xyz file has the number position lines, which is orbitals_per_uc + num_atoms
xyz_atomic_lines.append(str(int(num_wann/num_copies + len(atom_pos_primitive_cell))))
xyz_atomic_lines.append("Wannier centres for the primitive cell, enforced to the atomic positions")

# if spinors are used, we have to iterate over the atomic positions twice
if spinors:
    outer_range = 2
else:
    outer_range = 1

# iterate over the different atomic positions and add that atomic position orbitals_per_atom times to xyz_atomic_lines
for i in range(outer_range):
    for atom, positions in atom_pos_primitive_cell.items():
        for position in positions:
            for j in range(int(orbitals_per_atom[atom]/outer_range)):
                xyz_atomic_lines.append("X    {}  {}  {}".format(position[0], position[1], position[2]))

# iterate over the different atomic positions and add them
for atom, positions in atom_pos_primitive_cell.items():
    for position in positions:
        xyz_atomic_lines.append("{}    {}  {}  {}".format(atom, position[0], position[1], position[2]))

with open("configurations/point.{}/wannier90_files/primitive_cell_centres.xyz".format(point), "w") as xyz_out:
    xyz_out.write("\n".join(xyz_atomic_lines))

# generate list of tuples giving the start of an atom block (defined as same type, same spin), and the number of primitive orbitals in it
atom_block_list = []
# define start variable to be iteratively increased
start = 0

for i in range(outer_range):
    # iterate over all the atoms
    for atom, num_orbitals in orbitals_per_atom.items():
        atom_block_list.append((start, int(len(atom_pos_primitive_cell[atom])*num_orbitals/outer_range)))
        # the number of supercell entries corresponding to this atom will be num_atoms*orbitals_per_atom(/2, if spinors == True)
        start += int(len(atom_pos_primitive_cell[atom])*num_copies*num_orbitals/outer_range)

################################################## read, process, and write _hr.dat data ##################################################

# read data from corrections_hr.dat file
hr_in = open("configurations/point.{}/wannier90_files/correction_hr.dat".format(point), "r")

# skip comment line
hr_in.readline()

# next line contains the number of wannier orbitals per unit cell, the one after that the number of hopping vectors
num_wann = int(int(hr_in.readline())/num_copies)
hr_lines.append(str(num_wann))

# define dictionary to contain hopping elements: {vector: [num_wann*num_wann nested list]}
# each element of the nested list is itself a list of all the hopping strengths with that index - to be averaged
hopping_dict = {}

# set up that nested list
hopping_list = []
# make inner list
inner_list = []

for i in range(num_wann):
    inner_list.append(0)

for i in range(num_wann):
    hopping_list.append(inner_list[:])

num_vec = int(hr_in.readline())

# based on the number of vectors, determine how many lines to read past, using the fact that the file has 15 points per line
read_past = num_vec//15

if num_vec%15 != 0:
    read_past+=1
    # one additional line with < 15 elements must be read past if there is a not fully filled line

for i in range(read_past):
    hr_in.readline()

for line in hr_in:

    # read the hopping_vector, the from_orbital, the to_orbital, and the real and imaginary components of the strength
    hopping_data = line.split()
    hopping_vector_conventional = np.array(hopping_data[0:3], dtype=int)
    from_orbital = int(hopping_data[3])
    to_orbital = int(hopping_data[4])
    strength = complex(float(hopping_data[5]), float(hopping_data[6]))

    # determine hopping vector in the primitive basis
    hopping_vector_primitive = change_basis(hopping_vector_conventional, expansion_coefficients)

    # determine which atom set we are dealing with
    # set from_index and to_index to None to allow us to detect when they have been set
    from_index = None
    to_index = None

    # the correct index is that where the orbital is smaller than or equal to the start of the next block
    # if this doesn't happen by the last block, we know the orbital is in that block
    for i in range(len(atom_block_list)):
        if i == len(atom_block_list)-1:
            # checking for type rather than simply not from/to_index, as that returns True for 0 as well
            if not isinstance(from_index, int):
                from_index = i
            if not isinstance(to_index, int):
                to_index = i
        else:
            if from_orbital <= atom_block_list[i+1][0] and not isinstance(from_index, int):
                from_index = i
            if to_orbital <= atom_block_list[i+1][0] and not isinstance(to_index, int):
                to_index = i

        # terminate loop when both from_index and to_index have been determined
        if isinstance(from_index, int) and isinstance(to_index, int):
            break

    # determine vector connecting the primitive cell of the from_orbital to the primitive cell of the to_orbital
    offset_vector = hopping_vectors[(to_orbital-atom_block_list[to_index][0]-1)//atom_block_list[to_index][1]] -\
                    hopping_vectors[(from_orbital-atom_block_list[from_index][0]-1)//atom_block_list[from_index][1]]

    # add offset vector to hopping vector
    hopping_vector_primitive += offset_vector

    # map orbitals to primitive and add offset depending on which set of primitive orbitals the hopping points to
    # orbitals in wannier90_hr.dat are numbered 1-N, but labels in Python are 0-N-1 - this conversion includes that
    from_orbital = (from_orbital-atom_block_list[from_index][0]-1)%atom_block_list[from_index][1]
    to_orbital = (to_orbital-atom_block_list[to_index][0]-1)%atom_block_list[to_index][1]

    # add an offset of the numbers of primitive orbitals associated with all other atom blocks
    for i in range(from_index):
        from_orbital += atom_block_list[i][1]

    for i in range(to_index):
        to_orbital += atom_block_list[i][1]

    int_hopping_vector_primitive = hopping_vector_primitive.astype(int)
    # add hopping element to correct element of hopping_dict
    if not tuple(int_hopping_vector_primitive) in hopping_dict.keys():
        # set up that key-value pair
        hopping_dict[tuple(int_hopping_vector_primitive)] = deepcopy(hopping_list)
    if hopping_dict[tuple(int_hopping_vector_primitive)][from_orbital][to_orbital] == 0:
        hopping_dict[tuple(int_hopping_vector_primitive)][from_orbital][to_orbital] = [strength]
    else:
        hopping_dict[tuple(int_hopping_vector_primitive)][from_orbital][to_orbital].append(strength)

hr_in.close()

# iterate over hopping_dict and average all the elements
for hop_matrix in hopping_dict.values():
    for i in range(num_wann):
        for j in range(num_wann):
            hop_matrix[i][j] = np.mean(hop_matrix[i][j])

# write out all the data

# line indicating the number of different hopping vectors
hr_lines.append(str(len(hopping_dict)))

# set all the weights to 1
# iterate over complete lines (15 elements)
for i in range(len(hopping_dict)//15):
    hr_lines.append("    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1")

remaining_elements = len(hopping_dict)%15

if remaining_elements:
    # remaining_elements != 0 => there is a line that is not completely filled
    additional_line = "    1"
    for i in range(remaining_elements-1):
        additional_line += "    1"

    hr_lines.append(additional_line)

with open("configurations/point.{}/wannier90_files/primitive_cell_correction_hr.dat".format(point), "w") as hr_out:
    # write all the preamble lines
    hr_out.write("\n".join(hr_lines))
    hr_out.write("\n")

    # clear hr_lines for new elements to be added to it
    hr_lines = []

    # iterate over all hopping elements and write them to primitive_cell_corrections_hr.dat file
    for hop_vector, hop_matrix in hopping_dict.items():
        for i in range(num_wann):
            for j in range(num_wann):
                hr_lines.append("   {}   {}   {}   {}   {}   {}   {}".format(hop_vector[0], hop_vector[1], hop_vector[2], j+1, i+1,
                                hop_matrix[j][i].real, hop_matrix[j][i].imag))
        # add 1 to j and i to go from 0 to N-1 back to 1 to N

        # write after every hopping vector to avoid potential memory problems
        hr_out.write("\n".join(hr_lines))
        hr_out.write("\n")

        # clear hr_lines for new elements to be added to it
        hr_lines = []
