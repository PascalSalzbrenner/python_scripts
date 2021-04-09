# script to do various kinds of analysis and plots of the bonds in a crystal
# currently works using the bond analysis in a .castep output file of a CASTEP run
# I could write parsers for the structure input files of different DFT codes and then do the length / RDF type analyses for them

# Features:
# PDF: bins all bonds, divides by the number of atoms in the cell, and plots a histogram of frequency vs length
#      if a structure is given, essentially repeats the analysis of the .castep file
#      that is, PDF analysis is confined to unit cell
# RDF: proper RDF analysis - constructs 9x9x9 supercell from the structure (which can be read from the .castep or structure file)
#      and calculates proper RDF as defined in [reference here] by iterating over all atoms in the central cell and averaging
# bond length analysis: Takes the indices (bonds ranked by length from shortest to longest) of an arbitrary number of bonds and tracks their
#                       length(s) as a function of pressure
# bond population analysis: Likewise, but for the population. Only works from a .castep file

# written by Pascal Salzbrenner, pts28
