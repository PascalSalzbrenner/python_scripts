# script to take the output of a gstatement command run on CSD3 and sum the used core hours
# written by Pascal Salzbrenner, pts28@cam.ac.uk

filename = input("What is the name of the file containing the CSD3 'gstatement' data? ").strip()
time_file = open("{}".format(filename), "r")

core_hours = 0

# read the first line and check if it's the header - if not, we assume that the header was left out and all lines are job data
first_line = time_file.readline().split()

if "JobID" in first_line[0]:
    # header line - read past the next line also
    time_file.readline()
else:
    # assumed to give data on the first job
    core_hours += float(first_line[8])

# read rest of lines
for line in time_file:
    core_hours += float(line.split()[8])

time_file.close()

with open("total_time.dat", "w") as outfile:
    outfile.write("The total time for all calculations in {} was {} core hours.\n".format(filename, core_hours))
