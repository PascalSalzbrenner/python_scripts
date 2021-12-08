# script to read a column from a file and multiply or divide it by a user-specified number of factors of 10
# written by Pascal Salzbrenner, pts28

filename = input("Enter the filename: ")
operation = input("Do you want to multiply or divide: ")
factor = int(input("How many factors of 10 do you want to {} by: ".format(operation)))
column = int(input("Which column do you want to {} by {} factors of 10: ".format(operation, factor)))

if operation.lower().startswith("d"):
    # division requested - division by 10**n is equivalent to multiplication with 10**(-n)
    factor = -factor

infile = open("{}".format(filename), "r")
outfile =  open("{}_converted".format(filename), "w")

for line in infile:
    if line.startswith("#"):
        # comment lines are added on without processing
        outfile.write(line)
    else:
        elements = line.split()
        line_string = ""

        for i in range(len(elements)):
            if i == column-1:
                # the column we want to operate on
                line_string += "{} ".format(str(float(elements[i])*(10**(factor))))
            else:
                line_string += "{} ".format(elements[i])

        line_string = line_string.rstrip()
        line_string += "\n"

        outfile.write(line_string)

infile.close()
outfile.close()
