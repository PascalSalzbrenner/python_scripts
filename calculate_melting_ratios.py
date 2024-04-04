# script to calculate the ratio between melting points without and with SOC from Debye temperatures within the Lindemann theory
# we divide soc/no_soc
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
debye_files.sort()

# set up pressure and volume ratio dictionary
pressures_volume_ratio = {}

# read pressures and volumes
with open("volumes.dat", "r") as volumefile:
	# read past first line
	volumefile.readline()

	for line in volumefile:
		data = line.split()
		pressures_volume_ratio[int(data[0])] = (float(data[2])/float(data[1]))**(2/3)

# set up melting_ratios dictionary
# structure {structure_name: {pressure: {moment: melting_ratio}}}
melting_ratios = {}

# iterate over all files
for debye_file in debye_files:
	# use filename to get relevant information

	split_name = debye_file.split("_")

	# whether or not the data includes SOC
	if "no_soc" in debye_file:
		soc_yn = False
	else:
		soc_yn = True

	if soc_yn:
		max_iter = -3
	else:
		max_iter = -4

	# structure name
	structure_name = split_name[2]

	for name_part in split_name[3:max_iter]:
		if "p0" in name_part:
			add_name = "_" + "-".join(name_part.split("-")[1:])
		else:
			add_name = "_" + name_part
		structure_name += add_name

	print(structure_name)

	# moment
	moment = int(split_name[split_name.index("moment")-1])

	# pressure
	pressure = int(split_name[3].split("p")[0])

	# determine actual Debye temperature
	with open(debye_file, "r") as infile:
		# need the second line
		infile.readline()

		debye_temperature = float(infile.readline().split()[5])

	# add to dictionary
	if structure_name in melting_ratios.keys():
		if pressure in melting_ratios[structure_name].keys():
			if moment in melting_ratios[structure_name][pressure].keys():
				# we already have either the point without or with SOC
				if soc_yn:
					melting_ratios[structure_name][pressure][moment] = debye_temperature / melting_ratios[structure_name][pressure][moment]
				else:
					melting_ratios[structure_name][pressure][moment] /= debye_temperature
			else:
				melting_ratios[structure_name][pressure][moment] = debye_temperature
		else:
			melting_ratios[structure_name][pressure] = {moment: debye_temperature}
	else:
		melting_ratios[structure_name] = {pressure: {moment: debye_temperature}}

print(melting_ratios)

# print one file per structure
for structure_name in melting_ratios.keys():

	header_written = False

	outfile = open("lindemann_melting_temperature_ratio_{}.dat".format(structure_name), "w")

	for pressure in melting_ratios[structure_name].keys():

		if not header_written:

			header = "# Pressure [GPa]"

			for moment in melting_ratios[structure_name][pressure].keys():

				header += "; Melting ratio, {} moment".format(moment)

			outfile.write("{}\n".format(header))
			header_written = True

		out_line = "{}".format(pressure)

		for moment in melting_ratios[structure_name][pressure].keys():
			out_line += " {}".format(pressures_volume_ratio[pressure]*(melting_ratios[structure_name][pressure][moment]**2))

		outfile.write("{}\n".format(out_line))

	outfile.close()
















