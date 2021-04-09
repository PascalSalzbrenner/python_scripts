# the AIRSS ca -e command generates enthalpy plots with the enthalpy given in units of eV. This script converts to meV

# written by Pascal Salzbrenner, pts28@cam.ac.uk

fileroot = input("What is the fileroot of the .agr file? ")

infile = open("{}.agr".format(fileroot), "r")
outfile = open("{}_mev.agr".format(fileroot), "w")

for line in infile:

    if line.startswith("@") or line.startswith("&"):
        outfile.write(line)
    else:
        # every line that is not marked by one of these special characters contains data in the form "pressure enthalpy"
        data = line.split()
        outfile.write("{} {}\n".format(data[0], 1000*float(data[1])))

infile.close()
outfile.close()
