# short script to convert an energy, given in eV, to a temperature in Kelvin, via E=k_(B)T
# this is mostly so I don't have to look up the Boltzmann constant every time

# define Boltzmann constant [eV/K]
kb = 0.00008617333262

energy = float(input("What is the energy [eV]? "))

temp = energy / kb

print("The temperature is {} K.\n".format(temp))
