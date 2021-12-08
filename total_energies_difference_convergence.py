# script to compare the energies in the output files of two of my convergence test runs
# requires that these calculations be done with the same code and pseudopotential, and over the same range
# the first column contains the convergence parameter, the second the energy in eV/atom

# written by Pascal Salzbrenner, pts28@cam.ac.uk

structure_1 = input("What is the name of the directory containing the first calculation? ").rstrip("/")
structure_2 = input("What is the name of the directory containing the second calculation? ").rstrip("/")

filename = input("What is the name of the convergence test data file? ").rstrip("/")

first_energies = open("{}/{}".format(structure_1, filename), "r")
second_energies = open("{}/{}".format(structure_2, filename), "r")

comparison_file = open("{}_{}_{}_energy_difference.dat".format(structure_1.split("/")[0], structure_2.split("/")[0], filename.split(".")[0]), "w")

comparison_file.write("# difference in DFT total energies calculated in {} and {}\n".format(structure_1, structure_2))
comparison_file.write("# convergence parameter; energy difference ({}-{}) [meV/atom]\n".format(structure_2, structure_1))

# read past the first line in both files

first_energies.readline()
second_energies.readline()

for first_line in first_energies:
    second_line = second_energies.readline()

    first_data = first_line.split()
    second_data = second_line.split()

    if (len(first_data) >= 3) and (len(second_data) >= 3):
        # if this is not fulfilled, something went wrong in the FP calculation and the data point should be skipped
        comparison_file.write("{} {}\n".format(first_data[0], (float(second_data[1])*1000-float(first_data[1])*1000)))
    else:
        continue

first_energies.close()
second_energies.close()
comparison_file.close()
