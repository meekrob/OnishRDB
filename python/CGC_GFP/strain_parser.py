#!/usr/bin/env python3
import sys,csv
from ParseFusionStr import *

# use the csv module to parse the quoted csv- parsing on commas 
# won't work otherwise

with open('fixed-cgc-search-1685649968.csv') as cgc_in:
    header = next(cgc_in)
    for line in cgc_in:
        try:
            cgc_line = CGCStrainEntry( line )
        except ValueError:
            print("couldn't parse %s" % line.strip(), file=sys.stderr)

        if cgc_line.species_name.find('elegans') < 0: continue

        print(cgc_line.strain_name, cgc_line.genotype, cgc_line.description.replace("\t", " "), sep="\t")
