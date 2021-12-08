#DOI: 10.1016/j.cpc.2019.107080

import pyprocar
import os

if not os.path.exists("./plots"):
    os.mkdir("./plots")

# check if the band structure calculation has been split
ls = os.listdir(".")
ls.sort()

if "split-01" in ls:
    # has indeed been split
    split_list = [el for el in ls if "split" in el and "KPOINTS" not in el]
    procar_list = ["{}/PROCAR".format(split) for split in split_list]
    pyprocar.cat(procar_list, "PROCAR")

pyprocar.repair("PROCAR","PROCAR-repaired")

#plot the projected contributions of the different orbitals, resolved by species

pyprocar.bandsplot("PROCAR-repaired", outcar="OUTCAR", mode="parametric",
elimit=[-7.5,12.5], knames=["$\\Gamma$", "$X$", "$W$", "$K$", "$\\Gamma$", "$L$", "$U$", "$W$", "$L$", "$K/U$", "$X$"],
kticks=[0,99,199,299,399,499,599,699,799,899,999], atoms=[1],cmap="hot_r", discontinuities=[899],
orbitals=[0],savefig="./plots/band_structure_Cd_s.pdf",figsize=(20,10))

pyprocar.bandsplot("PROCAR-repaired", outcar="OUTCAR", mode="parametric",
elimit=[-7.5,12.5], knames=["$\\Gamma$", "$X$", "$W$", "$K$", "$\\Gamma$", "$L$", "$U$", "$W$", "$L$", "$K/U$", "$X$"],
kticks=[0,99,199,299,399,499,599,699,799,899,999], atoms=[2], cmap="hot_r", discontinuities=[899],
orbitals=[0],savefig="./plots/band_structure_Te_s.pdf",figsize=(20,10))

pyprocar.bandsplot("PROCAR-repaired", outcar="OUTCAR", mode="parametric",
elimit=[-7.5,12.5], knames=["$\\Gamma$", "$X$", "$W$", "$K$", "$\\Gamma$", "$L$", "$U$", "$W$", "$L$", "$K/U$", "$X$"],
kticks=[0,99,199,299,399,499,599,699,799,899,999], atoms=[1],cmap="hot_r", discontinuities=[899],
orbitals=[1,2,3],savefig="./plots/band_structure_Cd_p.pdf",figsize=(20,10))

pyprocar.bandsplot("PROCAR-repaired", outcar="OUTCAR", mode="parametric",
elimit=[-7.5,12.5], knames=["$\\Gamma$", "$X$", "$W$", "$K$", "$\\Gamma$", "$L$", "$U$", "$W$", "$L$", "$K/U$", "$X$"],
kticks=[0,99,199,299,399,499,599,699,799,899,999], atoms=[2], cmap="hot_r", discontinuities=[899],
orbitals=[1,2,3],savefig="./plots/band_structure_Te_p.pdf",figsize=(20,10))

pyprocar.bandsplot("PROCAR-repaired", outcar="OUTCAR", mode="parametric",
elimit=[-7.5,12.5], knames=["$\\Gamma$", "$X$", "$W$", "$K$", "$\\Gamma$", "$L$", "$U$", "$W$", "$L$", "$K/U$", "$X$"],
kticks=[0,99,199,299,399,499,599,699,799,899,999], atoms=[1],cmap="hot_r", discontinuities=[899],
orbitals=[3],savefig="./plots/band_structure_Cd_px.pdf",figsize=(20,10))

pyprocar.bandsplot("PROCAR-repaired", outcar="OUTCAR", mode="parametric",
elimit=[-7.5,12.5], knames=["$\\Gamma$", "$X$", "$W$", "$K$", "$\\Gamma$", "$L$", "$U$", "$W$", "$L$", "$K/U$", "$X$"],
kticks=[0,99,199,299,399,499,599,699,799,899,999], atoms=[2], cmap="hot_r", discontinuities=[899],
orbitals=[3],savefig="./plots/band_structure_Te_px.pdf",figsize=(20,10))

pyprocar.bandsplot("PROCAR-repaired", outcar="OUTCAR", mode="parametric",
elimit=[-7.5,12.5], knames=["$\\Gamma$", "$X$", "$W$", "$K$", "$\\Gamma$", "$L$", "$U$", "$W$", "$L$", "$K/U$", "$X$"],
kticks=[0,99,199,299,399,499,599,699,799,899,999], atoms=[1],cmap="hot_r", discontinuities=[899],
orbitals=[1],savefig="./plots/band_structure_Cd_py.pdf",figsize=(20,10))

pyprocar.bandsplot("PROCAR-repaired", outcar="OUTCAR", mode="parametric",
elimit=[-7.5,12.5], knames=["$\\Gamma$", "$X$", "$W$", "$K$", "$\\Gamma$", "$L$", "$U$", "$W$", "$L$", "$K/U$", "$X$"],
kticks=[0,99,199,299,399,499,599,699,799,899,999], atoms=[2], cmap="hot_r", discontinuities=[899],
orbitals=[1],savefig="./plots/band_structure_Te_py.pdf",figsize=(20,10))

pyprocar.bandsplot("PROCAR-repaired", outcar="OUTCAR", mode="parametric",
elimit=[-7.5,12.5], knames=["$\\Gamma$", "$X$", "$W$", "$K$", "$\\Gamma$", "$L$", "$U$", "$W$", "$L$", "$K/U$", "$X$"],
kticks=[0,99,199,299,399,499,599,699,799,899,999], atoms=[1],cmap="hot_r", discontinuities=[899],
orbitals=[2],savefig="./plots/band_structure_Cd_pz.pdf",figsize=(20,10))

pyprocar.bandsplot("PROCAR-repaired", outcar="OUTCAR", mode="parametric",
elimit=[-7.5,12.5], knames=["$\\Gamma$", "$X$", "$W$", "$K$", "$\\Gamma$", "$L$", "$U$", "$W$", "$L$", "$K/U$", "$X$"],
kticks=[0,99,199,299,399,499,599,699,799,899,999], atoms=[2], cmap="hot_r", discontinuities=[899],
orbitals=[2],savefig="./plots/band_structure_Te_pz.pdf",figsize=(20,10))

pyprocar.bandsplot("PROCAR-repaired", outcar="OUTCAR", mode="parametric",
elimit=[-7.5,12.5], knames=["$\\Gamma$", "$X$", "$W$", "$K$", "$\\Gamma$", "$L$", "$U$", "$W$", "$L$", "$K/U$", "$X$"],
kticks=[0,99,199,299,399,499,599,699,799,899,999], atoms=[1],cmap="hot_r", discontinuities=[899],
orbitals=[0,1,2,3],savefig="./plots/band_structure_Cd_sp3.pdf",figsize=(20,10))

pyprocar.bandsplot("PROCAR-repaired", outcar="OUTCAR", mode="parametric",
elimit=[-7.5,12.5], knames=["$\\Gamma$", "$X$", "$W$", "$K$", "$\\Gamma$", "$L$", "$U$", "$W$", "$L$", "$K/U$", "$X$"],
kticks=[0,99,199,299,399,499,599,699,799,899,999], atoms=[2], cmap="hot_r", discontinuities=[899],
orbitals=[0,1,2,3],savefig="./plots/band_structure_Te_sp3.pdf",figsize=(20,10))

#plot the projected contributions of the different orbitals, not resolved by species

pyprocar.bandsplot("PROCAR-repaired", outcar="OUTCAR", mode="parametric",
elimit=[-7.5,12.5], knames=["$\\Gamma$", "$X$", "$W$", "$K$", "$\\Gamma$", "$L$", "$U$", "$W$", "$L$", "$K/U$", "$X$"],
kticks=[0,99,199,299,399,499,599,699,799,899,999],cmap="hot_r", orbitals=[0], discontinuities=[899],
savefig="./plots/band_structure_s.pdf",figsize=(20,10))

pyprocar.bandsplot("PROCAR-repaired", outcar="OUTCAR", mode="parametric",
elimit=[-7.5,12.5], knames=["$\\Gamma$", "$X$", "$W$", "$K$", "$\\Gamma$", "$L$", "$U$", "$W$", "$L$", "$K/U$", "$X$"],
kticks=[0,99,199,299,399,499,599,699,799,899,999], cmap="hot_r", orbitals=[1,2,3], discontinuities=[899],
savefig="./plots/band_structure_p.pdf",figsize=(20,10))

pyprocar.bandsplot("PROCAR-repaired", outcar="OUTCAR", mode="parametric",
elimit=[-7.5,12.5], knames=["$\\Gamma$", "$X$", "$W$", "$K$", "$\\Gamma$", "$L$", "$U$", "$W$", "$L$", "$K/U$", "$X$"],
kticks=[0,99,199,299,399,499,599,699,799,899,999], cmap="hot_r", orbitals=[3],
savefig="./plots/band_structure_px.pdf",figsize=(20,10))

pyprocar.bandsplot("PROCAR-repaired", outcar="OUTCAR", mode="parametric",
elimit=[-7.5,12.5], knames=["$\\Gamma$", "$X$", "$W$", "$K$", "$\\Gamma$", "$L$", "$U$", "$W$", "$L$", "$K/U$", "$X$"],
kticks=[0,99,199,299,399,499,599,699,799,899,999], cmap="hot_r", orbitals=[1],
savefig="./plots/band_structure_py.pdf",figsize=(20,10))

pyprocar.bandsplot("PROCAR-repaired", outcar="OUTCAR", mode="parametric",
elimit=[-7.5,12.5], knames=["$\\Gamma$", "$X$", "$W$", "$K$", "$\\Gamma$", "$L$", "$U$", "$W$", "$L$", "$K/U$", "$X$"],
kticks=[0,99,199,299,399,499,599,699,799,899,999], cmap="hot_r", orbitals=[2],
savefig="./plots/band_structure_pz.pdf",figsize=(20,10))

pyprocar.bandsplot("PROCAR-repaired", outcar="OUTCAR", mode="parametric",
elimit=[-7.5,12.5], knames=["$\\Gamma$", "$X$", "$W$", "$K$", "$\\Gamma$", "$L$", "$U$", "$W$", "$L$", "$K/U$", "$X$"],
kticks=[0,99,199,299,399,499,599,699,799,899,999], cmap="hot_r", orbitals=[0,1,2,3], discontinuities=[899],
savefig="./plots/band_structure_sp3.pdf",figsize=(20,10))
