from __future__ import generators
from frowns.mdl_parsers import sdfile, sditerator
from frowns.perception import RingDetection, BasicAromaticity, sssr
from frowns.FrownsError import FrownsError
import sys, re
try:
    import cStringIO as StringIO
except ImportError:
    import StringIO

class sdin:
    def __init__(self, file, transforms=[sssr.sssr,
                                         BasicAromaticity.aromatize],
                 stripHydrogens=1):
        """(file)->molecule
        Convert a smiles string into a molecule representation"""
        self.file = file
        self.reader = sditerator.reader(file, stripHydrogens=stripHydrogens)
        self.transforms = transforms
        self.lastError = None
        self.mol_text = None

    def next(self):
        """-> returns (frowns, error, mol_record)
        if frowns is None check the error message"""
        mol, text, error = self.reader.next()
        if mol:
            error = error or ""
            for transform in self.transforms:
                try:
                    mol = transform(mol)
                except FrownsError, reason:
                    if error: error += "\n"
                    error += "Transformation Failed:%s"%reason
                    mol = None                
                except Exception, reason:
                    if error: error += "\n"
                    error += "Transformation Failed Unexpectedly:%s"%reason
                    mol = None
                
        return (mol, error, text)

    def __iter__(self):
        return self

def scan(file, blocksize=10000):
    """Get the number of molecules in an sd file"""
    pat = re.compile("^[$][$][$][$]", re.MULTILINE)

    text = file.read(blocksize)
    count = 0
    while text:
        g = pat.findall(text)
        count += len(g)
        if text[-1] == "$" and text[-4]!='$':
            next = text[-6:]
            text = "".join([next, file.read(blocksize)])
        else:
            text = text = file.read(blocksize)

    return count

FIELDPATTERN=re.compile(">\s+<([^>]+)>\s+\(*([^)]*)")
ALTFIELDPATTERN = re.compile(">\s+<([^>]+)>")
def sdIterate(file):
    """iterate through an sdfile returning (text, fields)
    pairs.  Where text is the text of the record and fields
    is a dictionary of fields in the record

    X requires generators"""
    pat = re.compile("^[$][$][$][$]")
    fieldpat = re.compile(">  <([^>]*)> [(]([^)]*)[)]")
    
    text = []
    fields = {}
    lines = iter(file)
    id = ""
    for line in lines:
        text.append(line)
        if pat.match(line):
            yield "".join(text), fields, id
            text = []
            fields = {}
            id = ""
        elif line[0] == ">":
            res = fieldpat.match(line)
            if res:
                name, id = res.groups()
                line = lines.next()
                text.append(line)
                fields[name] = line.strip()
            else:
                res = ALTFIELDPATTERN.match(line)
                if res:
                    name = res.groups()[0]
                    line = lines.next()
                    text.append(line)
                    fields[name] = line.strip()


            
            
        
