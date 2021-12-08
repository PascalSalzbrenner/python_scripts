# script fill in a user-supplied amount of space around a given dataset in a ca -e *-enthalpy.agr file
# adds on the new sets for plotting at the end

# written by Pascal Salzbrenner, pts28@cam.ac.uk

import math

fileroot = input("What is the fileroot of the .agr file containing the data? ")
offset = float(input("How much would you like to add? (NB: this value is added above as well as below the line) "))
set_index = input("What is the index of the data set you want to add {} to? ".format(offset))
max_index = int(input("What is the highest index currently in use in {}.agr? ".format(fileroot)))

infile = open("{}.agr".format(fileroot), "r")
outfile = open("{}_filled.agr".format(fileroot), "w")

# spacing of the lines for plotting - this seems a good default for my purposes but can easily be changed
spacing = 0.25

# determine the number of lines we have to plot
num_lines = math.ceil(offset/spacing)

# set up container for lines related to set
set_lines = []

for line in infile:

    # every line is written to output; the changes are added on at the end
    outfile.write(line)

    # check if the line is related to the set we offsetting
    if "s{}".format(set_index) in line:
        # this will be true the first time for the first related line

        # don't add the legend entry again
        if not "legend" in line:
            set_lines.append(line)

        for set_line in infile:

            outfile.write(set_line)
            set_lines.append(set_line)

            if set_line.startswith("&"):
                # marks the end of the block of lines related to the data set of interest
                break

infile.close()

# add output

for i in range(num_lines):

    # positive
    for line in set_lines:

        if "s{}".format(set_index) in line:
            outfile.write(line.replace("s{}".format(set_index), "s{}".format(max_index+1+2*i)))
        elif "&" in line:
            outfile.write(line)
        else:
            # data
            data = line.split()
            outfile.write("{} {}\n".format(data[0], float(data[1])+offset-i*spacing))

    # negative
    for line in set_lines:

        if "s{}".format(set_index) in line:
            outfile.write(line.replace("s{}".format(set_index), "s{}".format(max_index+2+2*i)))
        elif "&" in line:
            outfile.write(line)
        else:
            # data
            data = line.split()
            outfile.write("{} {}\n".format(data[0], float(data[1])-offset+i*spacing))

outfile.close()
