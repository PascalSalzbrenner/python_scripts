# Python script using TBmodels to symmetrise (10.1103/PhysRevMaterials.2.103805)
# wannier90 (10.1016/j.cpc.2007.11.016 & 10.1088/1361-648X/ab51ff) tight-binding models
# written by Pascal Salzbrenner - pts28@cam.ac.uk
# heavily indebted to the TBmodels tutorial http://z2pack.ethz.ch/tbmodels/doc/1.3/tutorial.html

# run in the directory where the w90 files already are

import tbmodels
import numpy as np
import pymatgen as mg
import symmetry_representation as sr

### TODO: if this works, there are several things that could be generalised:
# read coords from wannier90.win file
# species for pymatgen can also be read here
# potentially use the code I have already written (supercell_hopping_average.py) to generate the model position (as well as the coords)
# would no longer require orbitals_per_atom, which is only used to generate the positions
# read orbital types from wannier90.win file (I could restrict, at least in the first instance, how they can be named to simplify matching to sr.WANNIER_ORBITALS)
# read whether or not there are spinors from wannier90.win file
# a lot of this is very similar to the wannier90.win parser I wrote for supercell_hopping_average.py

orbitals_per_atom = 8

model = tbmodels.Model.from_wannier_files(hr_file="configuration_average_hr.dat",
                                          win_file="configurations/static/wannier90_files/wannier90.win")

coords = [[0,0,0], [0,0,0.5], [0,0.5,0], [0,0.5,0.5], [0.5,0,0], [0.5,0,0.5], [0.5,0.5,0], [0.5,0.5,0.5], [0.125,0.125,0.125],
          [0.125,0.125,0.625], [0.125,0.625,0.125], [0.125,0.625,0.625], [0.625,0.125,0.125], [0.625,0.125,0.625], [0.625,0.625,0.125], [0.625,0.625,0.625]]

pos_list = []

for i in range(2):
    for pos in coords:
        for j in range(int(orbitals_per_atom/2)):
            pos_list.append(pos)

model.pos = np.array(pos_list)

# add the sp2 and pz orbitals
orbital_type_list = []

orbital_type_list.extend(sr.WANNIER_ORBITALS["sp2"])
orbital_type_list.append("z")

# generate list of sr-type orbitals
orbitals = []

for spin in (sr.SPIN_UP, sr.SPIN_DOWN):
    for pos in coords:
        for orbital_type in orbital_type_list:
            orbitals.append(sr.Orbital(position=pos, function_string=orbital_type, spin=spin))

# get and apply time reversal symmetry
time_reversal = sr.get_time_reversal(orbitals=orbitals, numeric=True)
model_tr = model.symmetrize([time_reversal])

# save tr symmetrised model to _hr.dat file
model.to_hr_file("configuration_average_tr_hr.dat")

# space group symmetries
symmetries = []

structure = mg.Structure(lattice=model.uc, species=["Cd","Cd","Cd","Cd","Cd","Cd","Cd","Cd","Te","Te","Te","Te","Te","Te","Te","Te"],
coords=np.array(coords))

analyzer = mg.symmetry.analyzer.SpacegroupAnalyzer(structure)

sym_ops = analyzer.get_symmetry_operations(cartesian=False)

sym_ops_cart = analyzer.get_symmetry_operations(cartesian=True)

for sym, sym_cart in zip(sym_ops, sym_ops_cart):
    symmetries.append(sr.SymmetryOperation.from_orbitals(orbitals=orbitals, real_space_operator=sr.RealSpaceOperator.from_pymatgen(sym),
    rotation_matrix_cartesian=sym_cart.rotation_matrix, numeric=True))

model_full_sym = model_tr.symmetrize(symmetries, full_group=True)

# save fully symmetrised model to _hr.dat file
model_full_sym.to_hr_file("configuration_average_full_sym_hr.dat")
