"""file reader for
reading sdf (MDL) files
"""
from __future__ import generators

from frowns.Atom import Atom
from frowns.Bond import Bond
from frowns.Molecule import Molecule
import re, os, sys, traceback
try:
    from cStringIO import StringIO
except:
    from StringIO import StringIO
    
FIELDPATTERN=re.compile(">\s+<([^>]+)>\s+\(*([^)]*)")
ALTFIELDPATTERN = re.compile(">\s+<([^>]+)>")
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
        if hcount:
            print symbol, hcount
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
# sss bond stereo 0=Not stereo, 1=Up, 4=Either, 6=Down,
#                 Double Bonds:
#                   0 = Use x-y-z coords from atom block to determine cis/trans
#                   3 = Cis or trans (either) double bond

SINGLE_BOND_LOOKUP_STEREO = {0:None, 1:"UP", 4:"EITHER", 6:"DOWN"}
DOUBLE_BOND_LOOKUP_STEREO = {0:None, 3:"EITHER"}
TRIPLE_BOND_LOOKUP_STEREO = {0:None}
AROMATIC_BOND_LOOKUP_STEREO = {0:None}
BOND_LOOKUP_STEREO = [SINGLE_BOND_LOOKUP_STEREO, DOUBLE_BOND_LOOKUP_STEREO,
                      TRIPLE_BOND_LOOKUP_STEREO, AROMATIC_BOND_LOOKUP_STEREO]

def parse_bond(line):
    atom1, atom2, type, stereo, remainder = line[0:3], line[3:6], line[6:9], line[9:12], line[12:]
    atom1 = int(atom1) - 1
    atom2 = int(atom2) - 1
    type = int(type)
    stereo = int(stereo)
    return atom1, atom2, type, stereo, remainder

BOND_SYMBOL = {1:('-', 1, 1, 1),
               2:('=', 2, 2, 1),
               3:('#', 3, 3, 1),
               4:('', 4, 4, 1)}

class MolReaderError(Exception):
    pass


class collector:
    def __init__(self, file):
        self.iterator = iter(file)
        self.current = ""
        self.last = []

    def next(self):
        line = self.current = self.iterator.next()
        self.last.append(line)
        return line

    def dump(self):
        text = "".join(self.last)
        self.last = []
        return text
    
def reader(file, stripHydrogens=1):
    lines = collector(file)

    while 1:
        try:
            fields = {}
            name = lines.next().strip()
            userLine = lines.next().strip()
            comment = lines.next().strip()
            molinfo = lines.next()
            numAtoms, numBonds = int(molinfo[0:3]), int(molinfo[3:6])

            atoms = []   # this is the full list of atoms
            _atoms = []  # this is the (potentially stripped list
                         # of atoms.  I.e. no hydrogens.)
            i = 0
            for index in range(numAtoms):
                line = lines.next()
                x,y,z,symbol,mass,charge,stereo,hcount,hcount_fixed = parse_atom(line)
                if symbol == "H" and stripHydrogens:
                    atoms.append(None)
                else:
                    atom = Atom()
                    atoms.append(atom)
                    _atoms.append(atom)
                    atom.set_symbol(symbol)# = symbol
                    atom.explicit_hcount = hcount
                    atom.charge = charge
                    atom._line = line
                    atom.x = x
                    atom.y = y
                    atom.z = z
                    if hcount_fixed:
                        print "hcount fixed"
                        atom.fixed_hcount = 1 # oops, we shouldn't use this
                        atom.has_explicit_hcount = True
                    if mass:
                        atom.weight = atom.mass + mass

                    atom.index = i
                    i = i + 1

            bonds = []
            for index in range(numBonds):
                line = lines.next()
                a1, a2, bondtype, stereo, remainder = parse_bond(line)
                symbol, bondorder, bondtype, fixed = BOND_SYMBOL[bondtype]

                atom1, atom2 = atoms[a1], atoms[a2]
                
                if atom1 is not None and atom2 is not None:
                    h1, h2 = atom1.handle, atom2.handle
                    bond = Bond(symbol, bondorder, bondtype, fixed)
                    bonds.append(bond)
                    bond._line = remainder
                    bond.index = index
                    bond.atoms = [atom1, atom2]
                    try:
                        bond.stereo = BOND_LOOKUP_STEREO[bondtype-1][stereo]
                    except KeyError:
                        raise MolReaderError("An SD record cannot have a bondtype of %s and a stereo value of %s"%(bondtype, stereo))
                    except IndexError:
                        print "*"*44
                        print line
                        print "bondtype, stereo", bondtype, stero
                        raise
                        
                    atom1.bonds.append(bond)
                    atom2.bonds.append(bond)
                    atom1.oatoms.append(atom2)
                    atom2.oatoms.append(atom1)
                    if atom1.symbol == "H": atom2.explicit_hcount += 1
                    if atom2.symbol == "H": atom1.explicit_hcount += 1
                else:
                    if atom1 is None and atom2 is not None:
                        atom2.explicit_hcount += 1
                    elif atom2 is None and atom1 is not None:
                        atom1.explicit_hcount += 1
                        
            ##############################################################
            # read the mprops if necessary
            line = lines.next().strip()
            while 1:
                if line and line[0:6] == "M  END":
                    line = lines.next().strip()
                    break
                elif line == "M  CHG":
                    groups = line[6:].split()[1:]
                    index = 0
                    while index < len(groups):
                        atomIndex = int(groups[index]) - 1
                        atom = self.atoms[atomIndex]
                        charge = int(groups[index+1])
                        self.atoms[atomIndex].charge = charge
                        index += 2
                    line = lines.next().strip()
                elif line and line[0] == ">":
                    break
                elif line[0:4] == "$$$$":
                    break
                line = lines.next().strip()
                # What about end of mol?

            #############################################################
            # read the fields if necessary
            
            while line != "$$$$":
                if line and line[0] == ">":
                    res = FIELDPATTERN.match(line)
                    if res: field, potentialID = res.groups()
                    else:
                        res = ALTFIELDPATTERN.match(line)
                        if res:
                            field = res.groups()[0]
                            potentialID = None
                        else:
                            field, potentialID = None, None
                            
                    if name is None: name = potentialID

                    if field:
                        data = []
                        line = lines.next().strip()
                        while line and line != "$$$$":
                            data.append(line)
                            line = lines.next().strip()
                            
                        fields[field] = os.linesep.join(data)
                    
                line = lines.next().strip()
            mol = Molecule(_atoms, bonds)
            mol.name = name
            mol.fields = fields
            mol.name = name
            yield mol, lines.dump(), None
            
        except StopIteration:
            break
        except Exception:
            line = lines.current.strip()
            
            while line[0:4] != "$$$$":
                line = lines.next().strip()                
                
            stdout, stderr = sys.stdout, sys.stderr
            sys.stdout = sys.stderr = io = StringIO()
            traceback.print_exc()
            sys.stdout = stdout
            sys.stderr = stderr
            yield None, lines.dump(), io.getvalue()
            

