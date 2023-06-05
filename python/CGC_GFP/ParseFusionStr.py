import sys,re,csv
import pickle
from anyascii import anyascii

# a dict of WBIDs keyed on geneNames
with open("geneNames.pickle", "rb") as pickin:
    GENENAMES = pickle.load(pickin)

PAREN_RE = re.compile('\(.+\)')
TAG_RE = re.compile('tagRFP|RFP|mRFP|mNeonGreen|wCherry|vhhGFP4|FLAG|tdTomato|mRuby|eGFP|GFP|YFP|mCherry', re.IGNORECASE)

class CGCStrainEntry:
    def __init__(self, line):
        self.line = line
        self.fields = next( csv.reader([line], quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL) )
        try:
            self.strain_name, self.genotype, self.species_name, self.description = self.fields
            conv_description = anyascii(self.description.replace('’',"'").replace('—','-').replace('–','-') ) 
            if conv_description != self.description:
                self.description = conv_description
        except ValueError:
            print(line, file=sys.stderr)
            print(self.fields, file=sys.stderr)
            raise
        self.parsed_genotype = ''

        self.parseErrors = []
        self.constructs = [] # contains '::' delimiters
        self.other = [] # doesn't contain '::' delimiters

    def parse_genotype(self):
        self.parse_description_genotype()
        self.parse_constructs()

    def parse_constructs(self):
        parts = self.parsed_genotype.split(' + ')
        for i,part in enumerate(parts):
            if part.find('::') > -1:
                fusionParts = []
                for bit in part.split('::'):
                    fusionParts.append( FusionPart( bit ) )

                try:
                    fusion = Fusion( fusionParts )
                    self.constructs.append( fusion )
                except ValueError as v:
                    self.parseErrors.append( "Error parsing Fusion( fusionParts )") 
                    continue
            else:
                self.other.append( part )

    def getParseErrors(self):
        for error in self.parseErrors: 
            yield error
        for fusion_obj in self.constructs: 
            for parseError in fusion_obj.parseErrors:
                yield parseError
            

    def parse_description_genotype(self):
        # this is the nomenclature enclosed in brackets in the description field
        description = self.description
        try:
            leftBracket_i = description.index('[')
            rightBracket_i = description.index(']')
            self.parsed_genotype = description[(leftBracket_i+1):rightBracket_i]
        except ValueError:
            try: 
                leftBracket_i = self.genotype.index('[')
                rightBracket_i = self.genotype.index(']')
                self.parsed_genotype = self.genotype[(leftBracket_i+1):rightBracket_i]
            except ValueError:
                self.parseErrors.append( " ".join(["Can't parse genotype: unclosed brackets in description"]))


    def __str__(self):
        construct_str = ""
        if len(self.constructs) > 0:
            for construct in self.constructs:
                for part in construct:
                    construct_str += "\t\t" + str(part) + "\n"
        
        return """STRAIN: %s
\tGENOTYPE: %s
\tGENOTYPE (parsed): %s
\tSPECIES: %s
\tCONSTRUCTS: %d
%s
\tOTHER: %s
\tDESCRIPTION: %s
\tPARSE ERRORS: %s
""" % (self.strain_name, self.genotype, self.parsed_genotype, self.species_name, 
        len(self.constructs), 
        construct_str, 
        ",".join(self.other), 
        self.description, 
        "\n".join( list( self.getParseErrors()) )
        )

class Fusion:
    def __init__(self, fusionParts): # call with Fusion( genotype.split('::') )
        self.promoterWBID = None
        self.translationalWBID = None
        self.translationalName = None
        self.tag = None
        self.fusionParts = fusionParts
        self.description = ""
        self.parseErrors = []
        marked_promoter = list( filter( lambda x: x.is_transcriptional_fusion(), self.fusionParts ) )

        if len(marked_promoter) == 0:
            self.promoter = None
        elif len(marked_promoter) == 1:
            self.promoter = marked_promoter[0]
            self.promoterWBID = self.promoter.wbid
        else:
            error = "Unexpected parsed contents (multiple promoters): " + " ".join(list(map(str, marked_promoter)))
            self.parseErrors.append(error)
            print(self.parseErrors[-1], file=sys.stderr)
            raise ValueError("Multiple promoters")

        self.tags = list( filter( lambda x: x.is_tag(), self.fusionParts ) )
        self.fusion_genes = list( filter( lambda x: x.is_translational_fusion(), self.fusionParts ) )

    def getParseErrors(self):
        return self.parseErrors

    def getNames(self):
        for f in self.fusionParts: yield f.name

    def getGeneNames(self):
        for f in self.fusionParts: yield f.geneName

    def getTags(self):
        for f in self.fusionParts: yield f.tag

    def getLabels(self):
        for f in self.fusionParts: yield f.label

    def getTranslationalFusions(self):
        for f in filter( lambda x: x.is_translational_fusion(), self.fusionParts ):
            yield f

    def getPromoter(self):
        for f in filter( lambda x: x.is_transcriptional_fusion(), self.fusionParts ):
            yield f

    def __str__(self):
        return ";;".join(map(str,self.fusionParts))

    def __len__(self):
        return len(self.fusionParts)

    def __getitem__(self, i):
        return self.fusionParts[i]

class FusionPart:
    def __init__(self, part_str):
        self.transcriptional_fusion = False
        self.is_c_elegans_gene = False
        self.label = part_str
        self.unknown = True
        self.tag = None
        self.name = None
        self.geneName = None
        self.wbid = None
        self.parse()

    def remove_parenthetical(s):
        return ''.join( list( filter( lambda x: bool(x), PAREN_RE.split(s) ) ) )
        
    def set_txn_fusion(self, setting=True):
        self.transcriptional_fusion = setting

    def is_transcriptional_fusion(self):
        return self.transcriptional_fusion

    def is_translational_fusion(self):
        return (not self.transcriptional_fusion) and  (not self.is_tag())

    def is_tag(self):
        return bool(self.tag)

    def parse(self):
        self.name = FusionPart.remove_parenthetical(self.label).strip()

        if bool(TAG_RE.match(self.name)):
            self.tag = self.name

        elif self.name.endswith('p'): 
            self.transcriptional_fusion = True # there are no geneNames in the promoter table that end with 'p'
            self.wbid = self.findGeneID( self.name.rstrip('p') )
            if self.wbid is not None:
                self.geneName = self.name.rstrip('p')

        else:
            self.wbid = self.findGeneID( self.name )
            if self.wbid is not None:
                self.geneName = self.name
                
            
    def get_status(self):
        if self.is_transcriptional_fusion():
            if self.wbid is None:
                return "UNKNOWN p"
            else:
                return self.wbid + ' p'
        if self.is_tag():
            return "tag"
        else:
            if self.wbid is None:
                return "UNKNOWN"
            else:
                return self.wbid

        return "???" 

    def __str__(self):
        return "%s:%s" % (self.name, self.get_status())
        

    def findGeneID(self, query):
        try:
            return GENENAMES[query]
        except KeyError:
            return None

if __name__ == '__main__':
    print( Fusion( list( map(lambda x: FusionPart( x ), sys.argv[1].split('::') ) ) ) )
