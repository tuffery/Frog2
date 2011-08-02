"""file reader for
reading sdf (MDL) files
"""
from frowns.Atom import Atom
from frowns.Bond import Bond
from frowns.Molecule import Molecule
from frowns.perception import TrustedSmiles, RingDetection, BasicAromaticity
import re, os

# Atom format
#           111111111122222222223333333333444444444455555555556666666666
# 0123456789012345678901234567890123456789012345678901234567890123456789# 
#"xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee"
#"   -0.0187    1.5258    0.0104 C   0  0  0  0  0",
# x
# y
# z
#
# from the sdfile format distributed by mdl
# aaa atom symbol
# dd  mass difference -3,-2,-1,0,1,2,3,4 (0 if beyond these limits)
# ccc charge 0=uncharged or other than 1=+3, 2=+2, 3=+1, 4=doublet ^, 5=-1 6=-1 7=-3
# sss atom stereo parity 0 not stereo 1=odd, 2=even 3=either or unmarked stereo center
# hhh implicit hydrogen count + 1  1=HO, 2=H1, 3=H2, 4=H3, 5=H5 
# bbb
# vvv

def parse_atom(line):
    x, y, z = map(float, (line[0:10], line[10:20], line[20:30]))

    symbol =   line[31:34].strip()
    
    massdiff = line[34:36].strip()
    charge =   line[36:39].strip()
    stereo =   line[39:42].strip()
    hcount =   line[42:45].strip()

    if massdiff:
        massdiff = int(massdiff)
    else:
        massdiff = 0
        
    if charge:
        charge = int(charge)
        if charge: charge = 4 - charge
    else:
        charge = 0

    # if the hcount is set for the atom, 
    # we need to know this later
    #  i.e. the hcounts are explicit and can't be
    #  changed
    hcount_specified = 0
    if hcount:
        hcount = int(hcount)
        # an hcount of zero means it is not set
        # XXX FIX ME
        #  I don't have a way of forcing hcount to remain 0
        if hcount:
            hcount -= 1
            hcount_specified = 1
    else:
        hcount = 0
        
    # ignore stereo now
    stereo = 0
    return x, y, z, symbol, massdiff, charge, stereo, hcount, hcount_specified
#
# Bond Parser
#             11111111112 
#   012345678901234567890
#   111222tttsssxxxrrrccc
#
# 111 first atom index
# 222 second atom index
# ttt bond type 1=Single, 2=Double, 3=Triple, 4=Aromatic, 
#               Query types 5=Single or double, 6=Single or aromatic, 7=Double or aromatic, 8=Any
# sss bond stereo 0=Not stereo, 1=Up, 4=Either, 6=Down, Double Bonds
#                 0 = Use x-y-z coords from atom block to determine cis/trans
#                 3 = Cis or trans (either) double bond)

def parse_bond(line):
    atom1, atom2, type, remainder = line[0:3], line[3:6], line[6:9], line[9:]
    atom1 = int(atom1)
    atom2 = int(atom2)
    type = int(type)
    return atom1, atom2, type, remainder

BOND_SYMBOL = {1:('-', 1, 1, 1),
               2:('=', 2, 2, 1),
               3:('#', 3, 3, 1),
               4:('', 4, 4, 1)}

class MolReaderError(Exception):
    pass

class MolReader:
    def __init__(self, file, stripHydrogens=1):
        self.file = file
        self.iterator = iter(file)
        self._lastlines = []        # lastlines stores the original lines that made up the
                                            #  last molecule read
        self.stripHydrogens = stripHydrogens
        self._lastline = None
        
    def _readline(self, endOk=0):
        """internal readline function, if endOk is 0 then upon an end of
        line a MolReaderError is generated"""
        if self._lastline:
            res = self._lastline
            self._lastline = None
            return res

        try:
            line = self.iterator.next()
            self._lastlines.append(line)
        except StopIteration:
            line = None
            

        if not line and not endOk:
            raise MolReaderError, "Unexpected end of file"

        return line

    def _pushback(self, line):
        self._lastline = line
        
    def _clear(self):
        """Clear the _lastlines buffer"""
        self._lastlines = []

    def get_text(self):
        """->text that formed the last molecule read"""
        return "".join(self._lastlines)

    def get_lines(self):
        """->the lines of text that formed the last molecule read"""
        return self._lastlines
    
    def _read_to_next(self):
        readline = self._readline
        endOfMol = self._endOfMol
        
        while 1:
            line = readline(endOk=1)
            if not line:
                break

            if endOfMol(line):
                break


    def _endOfMol(self, line):
        """(line)-> return 1 if the line signifies the end of molecule
        0 otherwise"""
        if line[0:4] == "$$$$":
            return 1
        return 0

    def _readFields(self,
                    pattern=re.compile(">\s+<([^>]+)>\s+\(*([^)]*)")
                    ):
        """Read the database field component at the end of a molecule
        record.  Sets a dictionary of key->values"""
        readline = self._readline
        endOfMol = self._endOfMol
        
        fields = {}

        name = None
        while 1:
            # by setting endOk = 1 we can read mol files as
            # well as sdfiles
            line = readline(endOk=1)
            if not line:
                break
            
            if endOfMol(line):
                break
            elif line[0] == ">":
                # we have a data line so get the field
                #  and potentialID values
                if not endOfMol(line):
                    res = pattern.match(line)

                    if res:
                        field, potentialID = res.groups()
                    else:
                        field, potentialID = None, None

                if name is None:
                    name = potentialID
                elif name != potentialID:
                    name = "UNKNOWN (id clash)"
                    
                # read the data from the next line
                if field:
                    line = readline().strip()
                    data = []
                    while line:
                        data.append(line)
                        line = readline().strip()
                        
                    if not endOfMol(line):
                        fields[field] = os.linesep.join(data)
                    else:
                        break
                    
                if endOfMol(line):
                    break

        if not endOfMol(line):
            # by setting endok = 1 here we can read
            # mol files as well as sd files
            line = readline(endOk=1)
                 
        return fields, name

    def readMProps(self):
        readline = self._readline
        while 1:
            line = readline()
            if line[0] == ">":
                # need to push back the last line
                self._pushback(line)
                return
            
            if line[0:6] == "M  END":
                break

            if line[0:6] == "M  CHG":
                # parse the charge line and add charges
                # to the correct atoms
                groups = line[6:].split()[1:]
                index = 0
                while index < len(groups):
                    atomIndex = int(groups[index]) - 1
                    atom = self.atoms[atomIndex]
                    charge = int(groups[index+1])
                    self.atoms[atomIndex].charge = charge
                    index += 2
        
    def read_one(self):
        """Read one molecule from the sd file"""
        self._clear()
        readline = self._readline
        endOfMol = self._endOfMol
        try:
            name = readline().strip()
            userLine = readline()
            comment = readline()
            line = readline()
        except MolReaderError, msg:
            if str(msg) == "Unexpected end of file":
                return None            
            raise

        try:
            numAtoms, numBonds = map(int, (line[0:3], line[3:6]))
        except ValueError: # XXX FIX ME - trap exceptions and stuff
            print "cannot parse atom, bond line"
            self._read_to_next()
            return None
        atoms = self.atoms = []

        for index in range(numAtoms):
            line = readline()
            try:
                x,y,z,symbol,mass,charge,stereo,hcount,hcount_fixed = parse_atom(line)
            except: # XXX FIX ME - trap exceptions and stuff
                self._read_to_next()
                return None
            
            atom = Atom()
            atom._line = line
            atom.symbol = symbol
            atom.explicit_hcount = hcount
            atom.charge = charge

            #if hcount_fixed:                
            #symbol, hcount, charge, weight=0, aromatic=0)
            
            # XXX FIX ME
            # a really bad hack here.
            # ignore please!
            atom._line = line
            atom.x = x
            atom.y = y
            atom.z = z
            if hcount_fixed: atom.fixed_hcount = 1
            if mass:
                atom.weight = atom.mass + mass
            atom.index = len(atoms)
            atoms.append(atom)
            if vfgraph:
                insert_node(atom)

        bonds = []
        mappings = []
        bondCount = [0] * len(atoms)
        closures = {}
        for index in range(numBonds):
            line = readline()
            try:
                a1, a2, bondtype, remainder = parse_bond(line)
            except:
                self._read_to_next()
                return None
            a1 -= 1
            a2 -= 1

            symbol, bondorder, bondtype, fixed = BOND_SYMBOL[bondtype]
            atom1 = atoms[a1]
            atom2 = atoms[a2]
            if stripHydrogen:
                if atom1.symbol == "H":
                    atom2.hcount += 1
                    atom2.hcount_fixed = 1
                if atom2.symbol == "H":
                    atom1.hcount += 1
                    atom1.hcount_fixed = 1
            else:
                bond = Bond(symbol, bondorder, bondtype, fixed)
                bond._line = line
                # XXX FIX ME
                # a really bad hack here
                # ignore please!
                bond._line = remainder
                bond.atoms = [a1, a2]
                a1.bonds.append(bond)
                a2.bonds.append(bond)
                a1.oatoms.append(a1)
                a2.oatoms.append(a2)
                bonds.append(bond)
##            mappings.append((bond, a1, a2))
##                bondCount[a1] += 1
##                bondCount[a2] += 1

        self.readMProps()
        
        fields, potentialName = self._readFields()

        if not name:
            name = potentialName
        elif name != potentialName:
            # XXX FIX ME, what do I do here?
            pass
        
        # we've tokenized the molecule, now we need to build one
        # XXX FIX ME - Split this up into a builder and a tokenizer ?
        mol = Molecule()
        mol.name = name
        mol.fields = fields
#        for atom in atoms:
#            mol.add_atom(atom)
            
#        for bond, a1, a2 in mappings:
#            atom1, atom2 = atoms[a1], atoms[a2]

            # XXX FIX ME
            # does this format mean the atom's hcount can't
            #  change?
#            stripHydrogens = self.stripHydrogens
#            if not hasattr(atom1, "number"):
#                print atom
#                print atom.__dict__
##            if stripHydrogens and atom1.symbol == "H" and bondCount[a1] == 1:
##                atom2.hcount += 1
##                atom2.hcount_fixed = 1
##                bondCount[a2] -=1
##                mol.remove_atom(atom1)
##            elif stripHydrogens and atom2.symbol == "H" and bondCount[a2] == 1:
##                atom1.hcount += 1
##                atom1.hcount_fixed = 1
##                bondCount[a1] -= 1
##                mol.remove_atom(atom2)
##            else:
##                mol.add_bond(bond, atom1, atom2)
##                if bond.bondtype == 4:
##                    atom1.aromatic = 1
##                    atom2.aromatic = 1

##        # get rid of any non-bonded hydrogens
##        atomsToDelete = []
##        for atom in mol.atoms:
##            if atom.symbol == "H":
##                assert len(atom.bonds) == 0
##                atomsToDelete.append(atom)
##        for atom in atomsToDelete:
##            mol.remove_atom(atom)

        index = 0
        for atom in mol.atoms:
            assert atom.symbol != "H"
            if len(atom.bonds) > 1:
                atom._closure = 1
            atom.index = index
            index += 1

        return mol
                
                
            
def test():
    tests = [
	"   -0.0187    1.5258    0.0104 C   0  0  0  0  0",
	"    0.0021   -0.0041    0.0020 C   0  0  0  0  0",
	"    1.3951    2.0474   -0.0003 C   0  0  0  0  0",
	"    2.3236    1.2759   -0.0133 O   0  3  0  0  0",
	"    1.6511    3.5328   -0.0010 C   0  0  0  0  0",
	"    3.7321    1.7917   -0.0227 Cu  0  7  0  0  0",
	"    2.8773    3.8271   -0.8266 C   0  0  0  0  0",
	"    4.1685    2.0764    1.3837 O   0  3  0  0  0",
	"    4.6361    0.7604   -0.6303 O   0  3  0  0  0",
	"    3.7997    3.0530   -0.8316 O   0  3  0  0  0",
	"    2.9504    5.0907   -1.6446 C   0  0  0  0  0",
	"    5.3415    2.0641    1.6691 C   0  0  0  0  0",
	"    5.8107    0.7437   -0.3517 C  "
	]

    for line in tests:
        print parse_atom(line)

    bondTests = [
	"  1  2  1  0  0  0",
	"  1  3  1  0  0  0",
	"  1 20  1  0  0  0",
	"  1 21  1  0  0  0",
	"  2 22  1  0  0  0",
	"  2 23  1  0  0  0",
	]

    for line in bondTests:
	print parse_bond(line)


if __name__ == "__main__":
    test()

    file = open('../test/data/caffeine.mol')
    reader = MolReader(file)
    mol = reader.read_one()
    print mol.cansmiles()
