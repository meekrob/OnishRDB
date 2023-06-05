import sys,re
import pickle

# a dict of WBIDs keyed on geneNames
with open("geneNames.pickle", "rb") as pickin:
    GENENAMES = pickle.load(pickin)

PAREN_RE = re.compile('\(.+\)')
TAG_RE = re.compile('tagRFP|RFP|mRFP|mNeonGreen|wCherry|vhhGFP4|FLAG|tdTomato|mRuby|eGFP|GFP|YFP|mCherry', re.IGNORECASE)

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
        self.name = FusionPart.remove_parenthetical(self.label)

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
        return "%s:%s(%s)" % (self.name, self.wbid, self.get_status())
        

    def findGeneID(self, query):
        try:
            return GENENAMES[query]
        except KeyError:
            return None

if __name__ == '__main__':
    print( Fusion( list( map(lambda x: FusionPart( x ), sys.argv[1].split('::') ) ) ) )
