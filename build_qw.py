# Python script to take a file roughly in the shape of a VASP POSCAR giving the single conventional unit cell and convert it into a
# POSCAR for a quantum well
# NB: works for cubic materials
# written by Pascal Salzbrenner, pts28
# Input file format:
"""
Barrier_Material-Well_Material-Barrier_Material
Conventional lattice parameter
lattice_vector_x
lattice_vector_y
lattice_vector_z
VASP POSCAR-like block for first material (atom names, numbers, positions)
VASP POSCAR-like block for second material
"""
# the order of the well and barrier materials in these blocks is unimportant

import numpy as np

filename = input("Enter the name of your position file: ")

element_order = input("Enter the order of the elements in your POTCAR (A-B-C...): ")
element_list = element_order.split("-")

qw_growth_direction = int(input("Enter the direction in which the QW is constructed (1, 2, 3): "))

barrier_width = int(input("Enter the number of layers of the barrier material simulated explicitly on either side of the well: "))
well_width = int(input("Enter the number of layers of the well material: "))

infile = open("{}".format(filename), "r")
outfile = open("POSCAR", "w")

positions_dict = {} # dictionary to hold all the final positions
for element in element_list:
    # set up empty entries for each element
    positions_dict[element] = []

material_line = infile.readline()
outfile.write(material_line)

material_order = material_line.split()[0].split("-")
barrier_material = material_order[0]
well_material = material_order[1]

lattice_parameter = float(infile.readline().split()[0])
outfile.write("1.000000\n")

lattice_vectors = {}

for i in range(1, 4):

    lattice_vector_line = infile.readline()

    # store individual cell lattice vectors for potential later use
    lattice_vector = np.array([float(coordinate) for coordinate in lattice_vector_line.split()])
    lattice_vector = lattice_parameter * lattice_vector
    lattice_vectors["lattice_vector_{}".format(i)] = lattice_vector

    if i == qw_growth_direction:
        lattice_vector = (2*barrier_width + well_width)*lattice_vector

    outfile.write("{:.6f}  {:.6f}  {:.6f}\n".format(lattice_vector[0], lattice_vector[1], lattice_vector[2]))

barrier_dict = {}
well_dict = {}

for i in range(2):
    # there will always be exactly two blocks because a QW is made of two different materials

    atoms = infile.readline().split()
    atom_numbers =  infile.readline().split()
    # read the two lines giving element name and number, such that atom_numbers[i] gives the number of the atom in atoms[i]

    # distinguish between the blocks for well and barrier
    is_well_material = True

    for atom in atoms:
        if atom.lower() in well_material.lower():
            continue # we know it's the well material if every atom enters into this clause
        else:
            is_well_material = False
            break  # if one of the atoms isn't in well_material, the material must be the barrier material

    coordinates = infile.readline()

    if coordinates.lower().startswith("c") or coordinates.lower().startswith("k"):
        is_direct = False
    else:
        is_direct = True
        # NB: direct coordinates assumed as default - ie anything that isn't "Cartesian" or "Kartesisch" is taken as direct.

    if is_well_material:
        for atom, number in zip(atoms, atom_numbers):
            well_dict[atom] = []
            for j in range(int(number)):
                # number number of lines correspond to atom atom and are read as arrays into its well_dict
                well_dict[atom].append(np.array([float(coordinate) for coordinate in infile.readline().split()]))
    else:
        for atom, number in zip(atoms, atom_numbers):
            barrier_dict[atom] = []
            for j in range(int(number)):
                # number number of lines correspond to atom atom and are read as arrays into its well_dict
                barrier_dict[atom].append(np.array([float(coordinate) for coordinate in infile.readline().split()]))

# write the chemical symbols of the QW elements to the POSCAR in the right order
element_string = ""
for element in element_list:
    element_string += "{} ".format(element)

element_string = element_string.rstrip()
element_string += "\n"

outfile.write(element_string)

lattice_vector_growth = lattice_vectors["lattice_vector_{}".format(qw_growth_direction)]

# iterate over atoms in the well
for atom in well_dict.keys():
    for i in range(well_width):
        # we need to put all atoms once for each layer
        for position in well_dict[atom]:
            if is_direct:
                # convert to Cartesian
                append_position = np.zeros(3)
                for j in range(len(position)):
                    append_position += position[j]*lattice_vectors["lattice_vector_{}".format(j+1)]
            else:
                # Cartesian already
                append_position = position[:]

            positions_dict[atom].append(append_position+(barrier_width+i)*lattice_vector_growth)

# iterate over atoms in the barrier
for atom in barrier_dict.keys():
    for i in range(barrier_width):
        # we need to put all atoms once for each layer
        for position in barrier_dict[atom]:
            if is_direct:
                # convert to Cartesian
                append_position = np.zeros(3)
                for j in range(len(position)):
                    append_position += position[j]*lattice_vectors["lattice_vector_{}".format(j+1)]
            else:
                # Cartesian already
                append_position = position[:]

            positions_dict[atom].append(append_position+i*lattice_vector_growth)
            positions_dict[atom].append(append_position+(barrier_width+well_width+i)*lattice_vector_growth)
            # append the position in the barrier on both sides of the well layer

num_string = ""

# write the correct number of elements for the correct element
for element in element_list:
    num_string += "{} ".format(len(positions_dict[element]))

num_string = num_string.rstrip()
num_string += "\n"

outfile.write(num_string)
outfile.write("Cartesian\n")

# write the atomic positions to the POSCAR
for element in element_list:
    for position in positions_dict[element]:
        pos_string = "{:.6f}  {:.6f}  {:.6f}\n".format(position[0], position[1], position[2])
        outfile.write(pos_string)

infile.close()
outfile.close()
