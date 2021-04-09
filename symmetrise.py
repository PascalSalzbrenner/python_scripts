# Python script using TBmodels to symmetrise (10.1103/PhysRevMaterials.2.103805)
# wannier90 (10.1016/j.cpc.2007.11.016 & 10.1088/1361-648X/ab51ff) tight-binding models
# written by Pascal Salzbrenner - pts28
# heavily indebted to the TBmodels tutorial http://z2pack.ethz.ch/tbmodels/doc/1.3/tutorial.html

# run in the directory where the w90 files already are

import tbmodels
import numpy as np
import pymatgen as mg
import symmetry_representation as sr

material_name = input("Please enter the name of the material you are modelling: ")
basis_size = int(input("Please enter the number of basis orbitals: "))

model = tbmodels.Model.from_wannier_files(
    hr_file="wannier90_hr.dat",
    wsvec_file="wannier90_wsvec.dat",
    xyz_file="wannier90_centres.xyz",
    win_file="wannier90.win"
)

coords = [[0,0,0], [0.25,0.25,0.25]]

orbitals = []

# for a spinful model of size basis_size, the first basis_size/2 positions in
# model.pos are spin-up, the second half spin-down - generate list according to that

spins = []

for i in range(basis_size/2):
    spins.append(sr.SPIN_UP)

for i in range(basis_size/2):
    spins.append(sr.SPIN_DOWN)

# orbital positions can be read out from model.pos

orbital_type_list = []

# add the orbitals for spin up and down

for i in range(2):

    # add the atomic s- and p-orbitals for both atoms
    for j in range(2):

        orbital_type_list.append(sr.WANNIER_ORBITALS["s"][0])

        for orbital in sr.WANNIER_ORBITALS["p"]:
            orbital_type_list.append(orbital)

    # (basis_size-16)/2 orbitals are random s-like orbitals for spin up and down each
    # add them

    for k in range((basis_size-16)/2):
        orbital_type_list.append(sr.WANNIER_ORBITALS["s"][0])

for i in range(basis_size):
    orbitals.append(sr.Orbital(position=model.pos[i], function_string=orbital_type_list[i],
    spin=spins[i]))

# time reversal symmetry

time_reversal = sr.get_time_reversal(orbitals=orbitals, numeric=True)

# space group symmetries
# not sure if I can get this to work with the randomly positioned orbitals

symmetries = []

structure = mg.Structure(lattice=model.uc, species=["{}".format(material_name[0:2]), "{}".format(material_name[2:4])],
coords=np.array(coords))

analyzer = mg.symmetry.analyzer.SpacegroupAnalyzer(structure)

sym_ops = analyzer.get_symmetry_operations(cartesian=False)

sym_ops_cart = analyzer.get_symmetry_operations(cartesian=True)

for sym, sym_cart in zip(sym_ops, sym_ops_cart):
    symmetries.append(sr.SymmetryOperation.from_orbitals(orbitals=orbitals,
    real_space_operator=sr.RealSpaceOperator.from_pymatgen(sym),
    rotation_matrix_cartesian=sym_cart.rotation_matrix,numeric=True))


# to symmetrise, all the orbitals must be passed to the orbitals list.
# They are instances of the sr.Orbital class, which takes as arguments position,
# orbital (as defined in the sr module, as functions of x, y, and z), and spin
# (again as defined in the sr module, as instances of sr.SPIN_UP or sr.SPIN_DOWN)
# in the tutorial, they achieve this via are very nifty iteration

# I hope I have some way of getting the random orbital in the right form to be processed by the sr module
# given that it's a well-defined option in w90, I would hope sr implements them
# just use the sr s orbital at the random position? I think that's what the
# orbital is anyways

# the WCCs aren't at the positions of the atoms in the crystal, so instead of
# using model.pos in the pymatgen symmetry determination, I'll have to pass it
# my lattice explicitly

# at the end, save the new model to hdf5 and hr files for further use
model.to_hdf5_file("{}.hdf5".format(material_name))
model.to_hr_file("{}_hr.dat".format(material_name))
