# script to calculate the ratio between melting points without and with SOC from Debye temperatures within the Lindemann theory
# works with my debye_temperature.py scripts, and will calculate the ratio for any moment frequency Debye temperature present
# separately, a directory containing the volumes per atom at the relevant pressures must be provided
# format:
# pressure volume_no_soc volume_soc
# ...

# written by Pascal Salzbrenner, pts28@cam.ac.uk

import os

import numpy as np

# get all relevant files
debye_files = [debye_file for debye_file in os.listdir() if "debye_temperature" in debye_file]

# set up pressure and volume lists
pressures = []
volumes_no_soc = []
volumes_soc = []

# read pressures and volumes
with open("volumes.dat", "r") as volumefile:
	# read past first line
	volumefile.readline()

	for line in volumefile:
		data = line.split()
		pressures.append(int(data[0]))
		volumes_no_soc.append(float(data[1]))
		volumes_soc.append(float[data[2]])

volumes = [np.array(volumes_no_soc), np.array(volumes_soc)]

# set up debye_temperature dictionary
debye_temperature = {}

# iterate over all files
for debye_file in debye_files:
	# use filename to get relevant information

	split_name = debye_file.split("_")

	# structure name
	structure_name = "_".join(split_name[2:-4])

	# moment
	moment = int(split_name[split_name.index("moment")-1])

	# pressure
	pressure = int(split_name[3].split("p")[0])

	# whether or not the data includes SOC
	if "no_soc" in debye_file:
		soc_yn = False
	else:
		soc_yn = True

	# determine actual Debye temperature
	with open(debye_file, "r") as infile:
		# need the second line
		infile.readline()

		debye_temperature = float(infile.readline().split()[5])
















