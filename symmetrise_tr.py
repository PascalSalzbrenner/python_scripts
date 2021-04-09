# Python script using TBmodels to apply TR symmetry (10.1103/PhysRevMaterials.2.103805) to
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
w90_path = input("Please enter the path to the directory in which the wannier90 input files are: ")

w90_path = w90_path.rstrip("/")

model = tbmodels.Model.from_wannier_files(
    hr_file="{}/wannier90_hr.dat".format(w90_path),
    wsvec_file="{}/wannier90_wsvec.dat".format(w90_path),
    xyz_file="{}/wannier90_centres.xyz".format(w90_path),
    win_file="{}/wannier90.win".format(w90_path))

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

model_tr = model.symmetrize([time_reversal])

# combine all information about the initial model in an hdf5 file
model.to_hdf5_file("{}.hdf5".format(material_name))
# save the new model to hdf5 and hr files for further use
model_tr.to_hdf5_file("{}_tr.hdf5".format(material_name))
model_tr.to_hr_file("{}_tr_hr.dat".format(material_name))
