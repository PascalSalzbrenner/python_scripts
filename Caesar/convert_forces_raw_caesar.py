# file to take the forces as fetched from PROFESS and convert them to the format Caesar needs
# (should work for CASTEP as well but I haven't tested it)

# written by Pascal Salzbrenner, pts28@cam.ac.uk

import sys

atom = sys.argv[1]
disp = sys.argv[2]

input_file = open("forces_raw.dat", "r")
output = []

line_counter=1

for line in input_file:
    forces = line.split()

    for i in range(3):
        output.append("{} {} {} {} {}".format(atom, disp, line_counter, i+1, forces[i]))

    line_counter += 1

input_file.close()

with open("forces.dat", "w") as output_file:
    output_file.write("\n".join(output))
