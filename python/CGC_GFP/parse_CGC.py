#!/usr/bin/env python3

import sys,csv,re
from difflib import ndiff
from anyascii import anyascii
from ParseFusionStr import *

GENENAMES = {}

# put replacements in for typos or other things
FIXES = {}
FIXES['AX2164'] = ('GFP+', 'GFP +')
FIXES['TX1246'] = (' 771bp ','')

def parse_real_genotype(real_genotype): # the real_genotype is in brackets in the description field
    combo = map(lambda s: s.strip(), real_genotype.split('+'))
    is_genename = map(lambda g: g in GENENAMES, combo)
    

def parse_construct(s):
    if s.find('GFP') > -1:
        construct = s.split("::")
        ix = construct.index('GFP')

def cleanup(s):
    try:
        s = s.rstrip('p').rstrip('a')

        # find/remove a strain name in parentheses
        if s.startswith('('):
            scrape = ''.join( s.split(')').pop(0) )
        else:
            scrape = s.split('(')[0]
        return scrape.split()[0]
    except IndexError:
        print('fucked up', s, file=sys.stderr)
        raise

max_desc_length = 0
with open(sys.argv[1]) as CSV:
    lines = CSV.readlines()
    for line in csv.reader(lines, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL):
        fields = line
        strain_name,genotype,species_name,description = fields

        # Some Briggsae and E.coli are present in input
        if species_name.find('elegans') < 0: continue

        if strain_name in FIXES:
            search_str,replacement = FIXES[strain_name]
            description = description.replace( search_str, replacement )

        # pesky unicode dashes and ’
        """
        conv_description = anyascii(description.replace('’',"'").replace('—','-').replace('–','-') ) 
        if conv_description != description:
            diff = ndiff(description, conv_description)
            #print("changed:\n\t%s to\n\t%s" % (description,conv_description), file=sys.stderr)
            print("changed:", description, file=sys.stderr)
            for x in diff:
                if x.startswith('-') or x.startswith('+'): print(x, file=sys.stderr)
            description = conv_description
        """

        # truncate description for database        
        LIMIT = 60
        desc = description[:(LIMIT+1)]
        if len(description) > LIMIT:
            desc += '...'

        # find brackets that actually characterize the genotype
        try:
            leftBracket_i = description.index('[')
            rightBracket_i = description.index(']')
            real_genotype = description[(leftBracket_i+1):rightBracket_i]
        except ValueError:
            try: 
                leftBracket_i = genotype.index('[')
                rightBracket_i = genotype.index(']')
                real_genotype = genotype[(leftBracket_i+1):rightBracket_i]
            except ValueError:
                print("UNCLOSED BRACKETS FOR STRAIN", strain_name, "with DESCRIPTION:", description, file=sys.stderr)
                continue
        
        # ignore selection marker
        parts = real_genotype.split(' + ')
        print(strain_name + '>', real_genotype)
        for i,part in enumerate(parts):
            if part.find('::') > -1:

                fusionParts = []
                for bit in part.split('::'):
                    fusionParts.append( FusionPart( bit ) )

                try:
                    fusion = Fusion( fusionParts )
                    print("\t", i, fusion)
                except ValueError as v:
                    print("Error parsing strain", strain_name, file=sys.stderr)
                    print("Can't parse Fusion( %s )" % part, file=sys.stderr)
                    print("Description:", description, file=sys.stderr)
                    continue

            #else:
                #print("skipping part without '::'", strain_name, "with", part, file=sys.stderr) # a selection gene 


