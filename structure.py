# structure class to be used in any script that requires structures

# written by Pascal Salzbrenner, pts28@cam.ac.uk

import numpy as np
from exceptions import InputError
from utilities import lattice_basis_to_cartesian, cartesian_to_lattice_basis

class Structure:
    """Class to contain all the data about a structure one might find useful"""

    def __init__(self, filename):
        self.pressure = 0
        self.filename = filename
        self.filetype = filename.split(".")[-1]

        if self.filetype == "castep":
            self.get_structure_from_castep()
        elif self.filetype == "cell":
            self.get_structure_from_cell()
        else:
            # filetype that is not (currently implemented)
            raise InputError("Reading structure",
            "You have requested a structure for which no parsers is implemented. Currently supported: .castep, .cell files")

    def get_structure_from_cell():
        """Function to read the lattice vectors, atomic positions, and pressure (if present) from a CASTEP .cell file"""

        structure_file = open("{}".format(self.filename), "r")

        for line in structure_file:
            # interesting data are: lattice vectors, atom coordinates, and pressure

            if "lattice" in line.lower():
                # the next lines indicate the lattice, either in Cartesian coordinates or in ABC alpha beta gamma format

                # set up container for all lines relating to the lattice, to handle the optional coordinates indication
                lattice_lines = []

                for lattice_line in structure_file:

                    if "end" in lattice_line.lower():
                        # we have reached the end of the lattice block
                        break
                    else:
                        lattice_lines.append(line.rstrip("\n"))

                # differentiate between cart and ABC, with and without optional coordinates




def get_structure_from_castep
