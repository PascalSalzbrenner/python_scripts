# script to plot the output of a ddp-pes_all run, with the DFT at half opacity
# currently only works for a single element

# written by Pascal Salzbrenner, pts28@cam.ac.uk

import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# run in the directory containing the EDDP data

# define colour list
colours = ["#E6AB02", "#66A61E", "#8000C4", "#7570B3", "#E7298A", "#1E90FF"]

# define structure list & the corresponding labels
structures = ["dimer", "sc", "bcc", "fcc", "dc", "gh"]
labels_dict = {"bcc": "Bcc", "dc": "Diamond", "dimer": "Dimer", "fcc": "Fcc", "gh": "Graphene", "sc": "Simple cubic"}

# get xmin, xmax, ymin, ymax from commandline
xmin = sys.argv[1]
xmax = sys.argv[2]
ymin = sys.argv[3]
ymax = sys.argv[4]

# read element name
with open("elements.dat", "r") as element_file:
	element_name = element_file.readline().split()[0]

plt.xlabel("Interatomic Distance [A]")
plt.ylabel("Energy [eV/atom]")

plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)

for i in range(len(structures)):
	eddp_df = pd.read_csv("pes-{}-{}/pes.csv".format(element_name, structures[i]))
	eddp_distance = np.array(eddp_df["distance(A)"])
	eddp_energy = np.array(eddp_df["energy(eV/atom)"])

	dft_df = pd.read_csv("../DFT/pes-{}-{}/pes.csv".format(element_name, structures[i]))
	dft_distance = np.array(eddp_df["distance(A)"])
	dft_energy = np.array(eddp_df["energy(eV/atom)"])

	plt.plot(eddp_distance, eddp_energy, color=colours[i], linestyle="solid", label=labels_dict[structures[i]])
	plt.scatter(dft_distance, dft_energy, c=colours[i], marker=".", alpha=0.5)

plt.savefig("eddp_pes_w_dft.png", dpi=300)