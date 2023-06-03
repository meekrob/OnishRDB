import sys,re
import pickle

with open("geneNames.pickle", "rb") as pickin:
    GENENAMES = GN = pickle.load(pickin)

PAREN_RE = re.compile('\(.+\)')
TAG_RE = re.compile('tagRFP|RFP|mRFP|mNeonGreen|wCherry|vhhGFP4|FLAG|tdTomato|mRuby|eGFP|GFP|YFP|mCherry', re.IGNORECASE)


class Fusion:
    def __init__(self, fusionParts): # call with Fusion( genotype.split('::') )
        self.fusionParts = fusionParts
        marked_promoter = list( filter( lambda x: x.is_txn_fusion(), self.fusionParts ) )
        if len(marked_promoter) == 0:
            self.promoter = None
        elif len(marked_promoter) == 1:
            self.promoter = marked_promoter[0]
        else:
            print("Unexpected parsed contents (multiple promoters):", list(map(str, marked_promoter)), file=sys.stderr)
            raise ValueError("Multiple promoters")

        self.tags = list( filter( lambda x: x.is_tag(), self.fusionParts ) )
        self.fusion_genes = list( filter( lambda x: x.is_translational_fusion(), self.fusionParts ) )

    def __str__(self):
        return ";;".join(map(str,self.fusionParts))

class FusionPart:
    def __init__(self, part_str):
        self.transcriptional_fusion = False
        self.is_c_elegans_gene = False
        self.part_str = part_str
        self.unknown = False
        self.tag = False
        self.name = None
        self.parse()

    def remove_parenthetical(s):
        return ''.join( list( filter( lambda x: bool(x), PAREN_RE.split(s) ) ) )
        
    def set_txn_fusion(self, setting=True):
        self.transcriptional_fusion = setting

    def is_txn_fusion(self):
        return self.transcriptional_fusion

    def is_translational_fusion(self):
        return (not self.transcriptional_fusion) and (not self.unknown) and (not self.is_tag())

    def is_tag(self):
        return bool(self.tag)

    def parse(self):
        part_str = FusionPart.remove_parenthetical(self.part_str)

        if part_str.endswith('p') and not bool(TAG_RE.match(part_str)):
            self.transcriptional_fusion = True
            self.name = self.findGeneName( part_str.rstrip('p') )
            if self.name is None:
                self.name = part_str
                self.unknown = True
        else:
            self.name = self.findGeneName( part_str )
            if self.name is None:
                if bool(TAG_RE.match(part_str)): 
                    self.tag = part_str
                else:
                    self.unknown = True
        
                self.name = part_str

            
    def get_status(self):
        if self.is_txn_fusion():
            if self.unknown:
                return "UNKp"
            else:
                return "PROM"
        if self.is_tag():
            return "tag"
        if self.unknown:
            return "UNK"

        return "TRN" 

    def __str__(self):
        return "%s(%s)" % (self.name, self.get_status())
        

    def findGeneName(self, query):
        try:
            if query in GENENAMES:
                return query
            else:
                return None
        except KeyError:
            return None

if __name__ == '__main__':
    print( Fusion( list( map(lambda x: FusionPart( x ), sys.argv[1].split('::') ) ) ) )
