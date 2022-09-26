# structure class to be used in any script that requires structures

# written by Pascal Salzbrenner, pts28@cam.ac.uk

import numpy as np
from copy import deepcopy
from operator import itemgetter
from exceptions import InputError
from utilities import lattice_basis_to_cartesian, cartesian_to_lattice_basis

degree_to_radian = np.pi/180.0 # the angles in lattice_abc are given in degree, but numpy functions require them in radians

class Structure:
    """Class to contain all the data about a structure one might find useful"""

    def __init__(self, filename, frame=None):
        self.filename = filename.strip()
        self.filetype = filename.split(".")[-1]
        self.length_units = "Angstrom"
        self.pressure = 0
        self.pressure_units = "GPa"
        self.atoms = [] # contains the element names in the same sequence as the positions are in their list
        self.atom_numbers = [] # follows the numbering convention of CASTEP (all elements numbered from 1 to N_element)
                               # but is implemented for all structures
        self.frame=frame
        # MD trajectory files, such as .xyz(e) files contain the positions at different time steps
        # Here, we create a Structure instance for a specific frame
        # if no frame is given, this defaults to the final frame

        # structure parsers must set the following attributes: self.lattice, self.positions_abs, self.positions_frac, self.atoms,
        # self.atom_numbers
        # additionally, they can alter some of the properties set above

        if self.filetype == "castep":
            self.get_structure_from_castep()
        elif self.filetype == "cell" or self.filetype == "ion":
            self.get_structure_from_cell()
        elif self.filetype == "res":
            self.get_structure_from_res()
        elif self.filetype == "xyze":
            self.get_structure_from_xyze()
        elif self.filename == "structure.dat":
            # Caesar structure format
            self.get_structure_from_caesar_structure()
        else:
            # filetype that is not (currently) implemented
            raise InputError("Reading structure",
            "You have passed a structure for which no parsers is implemented. Currently supported: .castep, .cell, .res, .xyze files and the Caesar structure.dat")

        self.num_atoms = len(self.atoms)
        self.volume = np.linalg.det(np.dstack(self.lattice))[0]
        self.density = self.num_atoms / self.volume

    def __str__(self):
        """This is what will be printed when str(Structure) is called"""

        return self.filename

    def get_structure_from_cell(self):
        """Function to read the lattice vectors, atomic positions, and pressure (if present) from a CASTEP .cell file"""

        structure_file = open("{}".format(self.filename), "r")

        for line in structure_file:
            # interesting data are: lattice vectors, atom coordinates, and pressure

            if "lattice" in line.lower():
                # the next lines indicate the lattice, either in Cartesian coordinates or in ABC alpha beta gamma format

                # set up container for all lines relating to the lattice, to handle the optional units indication
                lattice_lines = []

                for lattice_line in structure_file:

                    if "end" in lattice_line.lower():
                        # we have reached the end of the lattice block
                        break
                    else:
                        lattice_lines.append(lattice_line.rstrip("\n"))

                # differentiate between cart and ABC, with and without optional units
                if "cart" in line.lower():

                    if len(lattice_lines) == 3:
                        # only the three lattice vectors
                        index_shift = 0
                    else:
                        # with three lattice vectors, this means one (the first) gives coordinates
                        self.length_units = lattice_lines[0]
                        index_shift = 1

                    self.lattice = np.array([lattice_lines[index_shift].split()[:3], lattice_lines[1+index_shift].split()[:3],
                    lattice_lines[2+index_shift].split()[:3]], dtype=float)
                else:
                    # abc format

                    if len(lattice_lines) == 2:
                        # only the two lines specifying the lattice are given
                        index_shift = 0
                    else:
                        # lattice_abc has one line giving the lengths of the vectors, and the second line the angles between them
                        # if there are three lines, the third element gives the units
                        self.length_units = lattice_lines[0]
                        index_shift = 1

                    vector_lengths = lattice_lines[index_shift].split()
                    angles = lattice_lines[index_shift+1].split()

                    self.construct_lattice_from_abc(vector_lengths, angles)

            elif "positions" in line.lower():
                # the next lines indicate the atomic positions
                # it may contain unit specifications (in the _abs case)
                # this script does not support units different from those in which the lattice vectors are given

                # initialise positions container
                positions = []

                # initialise container to count the numbers of the different atoms
                atoms_numbers = {}

                # initialise atom units
                atom_units = "Angstrom"

                for atoms_line in structure_file:

                    if "end" in atoms_line.lower():
                        # we have reached the end of the atoms block
                        break
                    else:
                        atom_data = atoms_line.split()

                        if len(atom_data) < 4:
                            # not perfect, because in principle I think CASTEP does not forbid comments, so the unit line could still be longer
                            # but this will be how the unit line is detected

                            atom_units = atom_data[0]
                        else:
                            self.atoms.append(atom_data[0])
                            positions.append(np.array(atom_data[1:4], dtype=float))

                            if atom_data[0] not in atoms_numbers.keys():
                                # first atom of this type
                                atoms_numbers[atom_data[0]] = 1
                            else:
                                atoms_numbers[atom_data[0]] += 1

                            self.atom_numbers.append(atoms_numbers[atom_data[0]])

                if "abs" in line.lower():
                    # positions directly in Cartesian coordinates
                    self.positions_abs = deepcopy(positions)

                    # set up fractional coordinates container
                    # it is not guaranteed that the lattice is passed before the positions, so we can't do the conversion here
                    self.positions_frac = []
                elif "frac" in line.lower():
                    # positions in fractional coordinates - do the inverse of the above
                    self.positions_frac = deepcopy(positions)

                    # set up Cartesian coordinates container
                    self.positions_abs = []

            elif "pressure" in line.lower():
                # the next lines indicate the upper half of the pressure matrix
                # again there is the possibility of optional units

                # set up lines container
                pressure_lines = []

                for pressure_line in structure_file:

                    if "end" in pressure_line.lower():
                        # we have reached the end of the pressure block
                        break
                    else:
                        pressure_lines.append(pressure_line.rstrip("\n"))

                if len(pressure_lines) == 3:
                    # no units given
                    self.pressure = float(pressure_lines[2].split()[0])
                else:
                    # units given
                    self.pressure_units = pressure_lines[0].split()[0]
                    self.pressure = float(pressure_lines[3].split()[0])

        # end of file reading
        structure_file.close()

        # check whether the positions were given in absolute or fractional coordinates
        if not self.positions_frac:
            # the positions were given in absolute coordinates
            # fill the fractional coordinates

            for atom in self.positions_abs:
                self.positions_frac.append(cartesian_to_lattice_basis(atom, self.lattice))

            # in this case, the positions have units - check for consistency with lattice units

            if atom_units.lower() not in self.length_units.lower() and self.length_units.lower() not in atom_units.lower():
                # note that this should be able to handle any situation
                # if the units are not angstrom, they should be given in both blocks identically
                # or if they are angstrom but given in both, they should also be identical
                # if they are angstrom and not given in either, then we'll be comparing the identical internal defaults
                # if they are given in one, they will be either "ang" or "angstrom", so one of the two conditions will always be true
                # thus, if this is false, there were inconsistent units, or they were denoted inconsistently
                # and why would you do that

                raise InputError("Length units",
                "Length units were given inconsistently between atoms and lattice. CASTEP may allow this (but why would you do it?), but this script does not.")

        else:
            # the positions were given in fractional coordinates - do the inverse of the above

            for atom in self.positions_frac:
                self.positions_abs.append(lattice_basis_to_cartesian(atom, self.lattice))

    def get_structure_from_castep(self):
        """Function to read the lattice vectors, atomic positions, and pressure from a CASTEP run .castep file
        Note that this assumes the task was a geometry optimisation - there is no reason to read from the .castep rather than the .cell file
        if a singlepoint calculation is carried out, and for MD postprocessing dedicated scripts (mdtep) already exist"""

        structure_file = open("{}".format(self.filename), "r")

        for line in structure_file:

            if "External pressure/stress" in line:
                # this block is only present once; the first element of the next line gives the pressure
                self.pressure = float(structure_file.readline().split()[0])
            elif "Final Configuration" in line:
                # we are interested in the final configuration in a geometry optimisation run

                # read past 6 lines
                for i in range(6):
                    structure_file.readline()

                # the first three elements of each of the next three lines contain the lattice vectors
                # initialise lattice vector list
                lattice_vectors = []

                for i in range(3):
                    lattice_vectors.append(structure_file.readline().split()[:3])

                self.lattice = np.array(lattice_vectors, dtype=float)

                # read past the next 18 lines to get to the atomic positions
                for i in range(18):
                    structure_file.readline()

                # initialise positions_frac - the .castep file always contains them like this
                self.positions_frac = []

                for atoms_line in structure_file:

                    if "xxx" in atoms_line:
                        # demarcates the end of the atoms block
                        break
                    else:
                        atoms_data = atoms_line.split()
                        self.atoms.append(atoms_data[1])
                        self.atom_numbers.append(atoms_data[2])
                        self.positions_frac.append(np.array(atoms_data[3:6], dtype=float))

        # end of file reading
        structure_file.close()

        # fill in positions_abs
        self.positions_abs = []
        for atom in self.positions_frac:
            self.positions_abs.append(lattice_basis_to_cartesian(atom, self.lattice))

    def get_structure_from_caesar_structure(self):
        """Function to read the required structure information from a Caesar structure.dat file"""

        self.length_units = "Bohr" # Caesar uses atomic units

        structure_file = open("{}".format(self.filename), "r")

        # the structure.dat format is very simple - the first line contains "lattice", and the next three lines contain a lattice vector each
        lattice_vectors = []

        structure_file.readline()

        for line in structure_file:
            if "atoms" in line.lower():
                # end of the lattice block
                break
            else:
                lattice_vectors.append(np.array(line.split()[:3]))

        self.lattice = np.array(lattice_vectors, dtype=float)

        # read the atomic names (converted into atomic numbers) and positions
        self.positions_abs = [] # this is how they are always given in structure.dat

        # initialise container to count the numbers of the different atoms
        atoms_numbers = {}

        for line in structure_file:

            if "symmetry" in line.lower():
                # end of the atoms block
                break
            else:
                atom_data = line.split()

                self.atoms.append(atom_data[0])
                self.positions_abs.append(np.array(atom_data[2:5], dtype=float))

                # determine correct atom number
                if atom_data[0] not in atoms_numbers.keys():
                    # first atom of this type
                    atoms_numbers[atom_data[0]] = 1
                else:
                    atoms_numbers[atom_data[0]] += 1

                self.atom_numbers.append(atoms_numbers[atom_data[0]])

        structure_file.close()

        # fill in positions_frac
        self.positions_frac = []
        for atom in self.positions_abs:
            self.positions_frac.append(cartesian_to_lattice_basis(atom, self.lattice))

    def get_structure_from_res(self):
        """Function to read the structure data from an AIRSS .res output file"""

        structure_file = open("{}".format(self.filename), "r")

        # the first line's third element is the pressure
        pressure_line = structure_file.readline()
        self.pressure = float(pressure_line.split()[2])

        # the next line is the one starting with "CELL", which gives the lattice in a b c alpha beta gamma format
        for line in structure_file:

            if line.startswith("CELL"):
                lattice_data = line.split()
                vector_lengths = [float(length) for length in lattice_data[2:5]]
                angles = [float(angle) for angle in lattice_data[5:8]]
                self.construct_lattice_from_abc(vector_lengths, angles)

            elif line.startswith("SFAC"):
                # the following lines, until the line "END", give the atoms and their positions (in Fractional_coordinates)

                # set up fractional coordinates container
                self.positions_frac = []

                # initialise container to count the numbers of the different atoms
                atoms_numbers = {}

                for atoms_line in structure_file:

                    if atoms_line.startswith("END"):
                        break
                    else:
                        atom_data = atoms_line.split()
                        # first element is the name of the atom, third to fifth element give its coordinates
                        self.atoms.append(atom_data[0])
                        self.positions_frac.append(np.array(atom_data[2:5], dtype=float))

                        if atom_data[0] not in atoms_numbers.keys():
                            # first atom of this type
                            atoms_numbers[atom_data[0]] = 1
                        else:
                            atoms_numbers[atom_data[0]] += 1

                        self.atom_numbers.append(atoms_numbers[atom_data[0]])

        structure_file.close()

        # convert to absolute coordinates
        self.positions_abs = []

        for atom in self.positions_frac:
                self.positions_abs.append(lattice_basis_to_cartesian(atom, self.lattice))

    def get_structure_from_xyze(self):
        """Function to read and average structure data from an MD trajectory .xyze file. Specifically, this has been designed to work with
           the ramble MD code, although the format is fairly standard and should translate fairly easily to other codes.
        """

        structure_file = open("{}".format(self.filename), "r")

        # initialise counter variable to keep track of frames
        frame_counter = 1

        # read number of atoms
        num_atoms = int(structure_file.readline().split()[0])

        # rewind file
        structure_file.seek(0)

        if not self.frame:
            # no frame has been given - determine the last frame of the file
            # the number of frames is the number of lines divided by the num_atoms+2 - for the atom number line and the lattice line
            num_lines = 0

            # count lines
            for line in structure_file:
                num_lines += 1

            self.frame = int(num_lines/(num_atoms+2))

        # rewind file
        structure_file.seek(0)

        # initialise container to count the numbers of the different atoms
        atoms_numbers = {}

        # find and read the right frame
        # read lattice vectors and atom positions - always given in Cartesian coordinates
        for line in structure_file:

            if frame_counter < self.frame:
                # we have not reached the right frame yet - read the next num_atoms + 2 lines
                for i in range(num_atoms+1):
                    structure_file.readline()

                # increment frame counter
                frame_counter += 1
            else:
                # we are at the frame we want to read - note we are already at the line with the atom number
                # next line contains the info about the lattice vectors
                lattice_vector_data_raw = structure_file.readline().split()

                # the last element of the third lattice vector will have a trailing " - remove it
                lattice_vector_data = [el.rstrip('"') for el in lattice_vector_data_raw]

                # create set of lattice vectors
                self.lattice = np.array(lattice_vector_data[1:10], dtype=np.double).reshape((3,3))

                # the next num_atom lines will contain the atomic positions
                # set up lists to contain the positions of all the atoms at the frame

                # read over all atoms at this frame
                for i in range(num_atoms):
                    atom_data = structure_file.readline().split()
                    self.atoms.append(atom_data[0])

                    if atom_data[0] not in atoms_numbers.keys():
                        # first atom of this type
                        atoms_numbers[atom_data[0]] = 1
                    else:
                        atoms_numbers[atom_data[0]] += 1

                    self.atom_numbers.append(atoms_numbers[atom_data[0]])

                    self.positions_abs.append(np.array(atom_data[1:4], dtype=np.double))

                break

        # convert to fractional coordinates
        self.positions_frac = []

        for atom in self.positions_abs:
            self.positions_frac.append(cartesian_to_lattice_basis(atom, self.lattice))

    def construct_lattice_from_abc(self, vector_lengths, angles):
        """Given lattice data in the a b c alpha beta gamma format, this function constructs the lattice vectors
           see the explanation for determining lattice vectors from this format at
           https://en.wikipedia.org/wiki/Fractional_coordinates#In_crystallography

           :param list vector_lengths: a list containing [a, b, c]
           :param list angles: a list containing [alpha, beta, gamma] in degrees
        """

        v_1 = [float(vector_lengths[0]), 0, 0]
        v_2 = [float(vector_lengths[1])*np.cos(float(angles[2])*degree_to_radian),
        float(vector_lengths[1])*np.sin(float(angles[2])*degree_to_radian), 0]
        v_3_x = float(vector_lengths[2])*np.cos(float(angles[1])*degree_to_radian)
        v_3_y = float(vector_lengths[2])*(np.cos(float(angles[0])*degree_to_radian)-np.cos(float(angles[2])*degree_to_radian)
                *np.cos(float(angles[1])*degree_to_radian))/(np.sin(float(angles[2])*degree_to_radian))
        v_3_z = np.sqrt(float(vector_lengths[2])**2-v_3_x**2-v_3_y**2)

        self.lattice = np.array([v_1, v_2, [v_3_x, v_3_y, v_3_z]])

    def get_bonds(self):
        """Function to do a bond length analysis of the parsed structure
        Creates bonds list containing tuples ("atom_name_number_1-atom_name_number_2", length).
        The bonds list generated from a .castep file in the get_bonds_from_castep function in bond_analysis.py generates tuples with the
        population as a third element. As calculating them requires the wavefunctions, it can't be done simply from the structure, but in
        tasks that only require the length, either the output of this function or of get_bonds_from_castep can be easily and equivalently
        used, with the bond length accessed by bonds[bond_index][1]."""

        # set up the bonds list
        self.bonds = []

        # iterate over every atom with a higher index (to avoid double counting) for each atom
        for i in range(len(self.atoms)):
            for j in range(i+1, len(self.atoms)):

                # note that the "real" distance between two atoms is that which is the shortest
                # this does not correspond necessarily to the distance within the unit cell, but might be between i and one of the periodic
                # images of j
                # trivial implementation (slow, but stable): iterate over all possible (sensible) periodic images of atom j

                min_vector = self.positions_abs[j]-self.positions_abs[i]
                min_vector_len = np.sqrt(min_vector.dot(min_vector))

                # determine in which directions it is sensible to check - those where point i is closer to the next unit cell
                additional_variables = []
                for coordinate in self.positions_frac[i]:
                    if coordinate <= 0.5:
                        additional_variables.append(-1)
                    else:
                        additional_variables.append(1)

                for l in [0, additional_variables[0]]:
                    for m in [0, additional_variables[1]]:
                        for n in [0, additional_variables[2]]:
                            cell_shift = np.array([l, m, n])
                            difference_vector = min_vector + np.dot(self.lattice.T, cell_shift)
                            difference_vector_len = np.sqrt(difference_vector.dot(difference_vector))

                            if difference_vector_len < min_vector_len:
                                # we have found a shorter image
                                min_vector_len = difference_vector_len

                self.bonds.append(("{}{} - {}{}".format(self.atoms[i], self.atom_numbers[i], self.atoms[j],
                self.atom_numbers[j]), min_vector_len))

        # sort by bond length
        self.bonds.sort(key=itemgetter(1))

        # create dictionary containing bond names as keys and a length 1 tuple (for consistent access with the bonds read directly from the
        # .castep file) containing the bond length as the value
        self.bonds_dict = {}

        for bond in self.bonds:
            self.bonds_dict[bond[0]] = (bond[1])
