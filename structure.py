# structure class to be used in any script that requires structures

# written by Pascal Salzbrenner, pts28@cam.ac.uk

import numpy as np
from copy import deepcopy
from exceptions import InputError
from utilities import lattice_basis_to_cartesian, cartesian_to_lattice_basis

class Structure:
    """Class to contain all the data about a structure one might find useful"""

    def __init__(self, filename):
        self.pressure = 0
        self.filename = filename
        self.filetype = filename.split(".")[-1]
        self.length_units = "Angstrom"
        self.pressure_units = "GPa"
        self.atoms = [] # sequence will be the same as the list for the positions

        if self.filetype == "castep":
            self.get_structure_from_castep()
        elif self.filetype == "cell":
            self.get_structure_from_cell()
        else:
            # filetype that is not (currently) implemented
            raise InputError("Reading structure",
            "You have passed a structure for which no parsers is implemented. Currently supported: .castep, .cell files")

    def get_structure_from_cell():
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

                    self.lattice = np.array([lattice_lines[index_shift].split()[:3], lattice_lines[1+index_shift].split()[:3]],
                    lattice_lines[2+index_shift].split()[:3], dtype=float)
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

                    # see the explanation for determining lattice vectors from this format at
                    # https://en.wikipedia.org/wiki/Fractional_coordinates#In_crystallography
                    v_1 = [float(vector_lengths[0]), 0, 0]
                    v_2 = [float(vector_lengths[1])*np.cos(float(angles[2])), float(vector_lengths[1])*np.sin(float(angles[2])), 0]
                    v_3_x = float(vector_lengths[2])*np.cos(float(angles[1]))
                    v_3_y = float(vector_lengths[2])*(np.cos(float(angles[0]))-np.cos(float(angles[2]))*np.cos(float(angles[1])))/(np.sin(float(angles[2])))
                    v_3_z = np.sqrt(float(vector_lengths[2])**2-v_3_x**2-v_3_y**2)

                    self.lattice = np.array([v_1, v_2, [v_3_x, v_3_y, v_3_z]])

            elif "positions" in line.lower():
                # the next lines indicate the atomic positions
                # it may contain unit specifications (in the _abs case)
                # this script does not support units different from those in which the lattice vectors are given

                # initialise positions container
                positions = []

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

                if "abs" in line.lower():
                    # positions directly in Cartesian coordinates
                    self.positions_abs = deepcopy(positions)

                    # set up fractional coordinates container
                    # it is not guaranteed that the lattice is passed before the positions, so we can't do the conversion here
                    self.positions_frac = []
                elif "frac" line.lower():
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
                    self.pressure = float(pressure_lines[2]) # single element due to the CASTEP pressure convention
                else:
                    # units given
                    self.pressure_units = pressure_lines[0].split()[0]
                    self.pressure = float(pressure_lines[3])

        # end of file reading

        # check whether the positions were given in absolute or fractional coordinates
        if not self.positions_frac:
            # the positions were given in absolute coordinates
            # fill the fractional coordinates

            for atom in self.positions_abs:
                self.positions_frac.append(cartesian_to_lattice_basis(atom, self.lattice))

            # in this case, the positions have units - check for consistency with lattice units

            if atom_units.lower() not in self.length_units.lower() and self.length_units.lower not in atom_units.lower():
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
            # the positions were given in absolute coordinates - do the inverse of the above

            for atom in self.positions_frac:
                self.positions_abs.append(lattice_basis_to_cartesian(atom, self.lattice))

    def get_structure_from_castep():
        """Function to read the lattice vectors, atomic positions, and pressure from a CASTEP run .castep file"""
