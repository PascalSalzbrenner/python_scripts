# script to compare the energies in the interpolated_free_energy.dat files of two different grid densities
# requires that these calculations be done with the same temperature step
# the first column contains the temperature, the third the vibrational energy in eV

# written by Pascal Salzbrenner, pts28@cam.ac.uk

# get names of the supercell directories

supercell_1 = input("What is the name of the directory containing the first calculation? ").rstrip("/")
supercell_2 = input("What is the name of the directory containing the second calculation? ").rstrip("/")

first_energies = open("{}/lte/interpolated_free_energy.dat".format(supercell_1), "r")
second_energies = open("{}/lte/interpolated_free_energy.dat".format(supercell_2), "r")

comparison_file = open("{}_{}_interpolated_free_energy_difference.dat".format(supercell_1, supercell_2), "w")

comparison_file.write("# difference in free energies calculated in {} and {}\n".format(supercell_1, supercell_2))
comparison_file.write("# T [K]; |energy difference| [meV]\n")

for first_line in first_energies:
    second_line = second_energies.readline()
    comparison_file.write("{} {}\n".format(first_line.split()[0], abs((float(second_line.split()[2])*1000-float(first_line.split()[2])*1000))))

first_energies.close()
second_energies.close()
comparison_file.close()
