# script to calculate the linear thermal expansion coefficient 1/a * \delta a / \delta T
# requires as input a minimum_energy.txt file as output by my birch_murnaghan_eos_fit_thermal_expansion.py script
# written by Pascal Salzbrenner, pts28

infile = open("minimum_energy.txt", "r")
outfile = open("thermal_expansion_coefficient.txt", "w")

outfile.write("#T [K], linear thermal expansion coefficient [MK^-1]\n")

for line in infile:

    if line.startswith("#"):
        # first line - the next line will be the static line - ignored in this calculation
        infile.readline()
        # the following line is the first "proper" one
        first_T_line = infile.readline()
        prev_T = int(first_T_line.split()[0].split("_")[1].rstrip("K"))
        prev_lattice_param = float(first_T_line.split()[1])
    else:
        cur_T = int(line.split()[0].split("_")[1].rstrip("K"))
        cur_lattice_param = float(line.split()[1])
        thermal_expansion_coefficient = (cur_lattice_param*(10**6)-prev_lattice_param*(10**6))/((cur_T-prev_T)*prev_lattice_param)
        # multiply the numerator by 10**6 because the thermal expansion values tend to be quite small and we don't want to conflict with
        # potential numerical inaccuracies

        # write output
        outfile.write("{} {}\n".format(prev_T, thermal_expansion_coefficient))

        # output values
        prev_T = cur_T
        prev_lattice_param = cur_lattice_param

infile.close()
outfile.close()
