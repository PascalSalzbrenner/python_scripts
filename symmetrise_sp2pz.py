# Python script using TBmodels to symmetrise (10.1103/PhysRevMaterials.2.103805)
# wannier90 (10.1016/j.cpc.2007.11.016 & 10.1088/1361-648X/ab51ff) tight-binding models
# written by Pascal Salzbrenner - pts28
# heavily indebted to the TBmodels tutorial http://z2pack.ethz.ch/tbmodels/doc/1.3/tutorial.html

# run in the directory where the w90 files already are

import tbmodels
import numpy as np
import pymatgen as mg
import symmetry_representation as sr

orbitals_per_atom = 8

model = tbmodels.Model.from_wannier_files(
    hr_file="wannier90_hr.dat",
    wsvec_file="wannier90_wsvec.dat",
    win_file="wannier90.win")

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
model.to_hr_file("wannier90_tr_hr.dat")

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
model_full_sym.to_hr_file("wannier90_full_sym_hr.dat")
