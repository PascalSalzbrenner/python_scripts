# script to add a certain constant value to all points of a given dataset in a ca -e *-enthalpy.agr file
# adds on the new set for plotting at the end

# written by Pascal Salzbrenner, pts28@cam.ac.uk

fileroot = input("What is the fileroot of the .agr file containing the data? ")
offset = float(input("How much would you like to add? "))
set_index = input("What is the index of the data set you want to add {} to? ".format(offset))
max_index = int(input("What is the highest index currently in use in {}.agr? ".format(fileroot)))

infile = open("{}.agr".format(fileroot), "r")
outfile = open("{}_offset.agr".format(fileroot), "w")

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

            if set_line.startswith("&"):
                # marks the end of the block of lines related to the data set of interest
                set_lines.append(set_line)
                break
            elif set_line.startswith("@"):
                # some instructions, not a data point

                # don't add the legend entry again
                if not "legend" in set_line:
                    set_lines.append(set_line.replace("s{}".format(set_index), "s{}".format(max_index+1)))
            else:
                # data point
                data = set_line.split()
                set_lines.append("{} {}\n".format(data[0], float(data[1])+offset))

infile.close()

# add output
for line in set_lines:
    outfile.write(line)

outfile.close()
