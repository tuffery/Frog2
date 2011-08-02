"""Experimental SDFILE Writer

No support for fields like
M  CHG
M  STY
and the like
"""

import os
# float format
#xxxxx.xxxx

class SDFloat:
    """Class to print out a float in the format the sd file
    expects"""
    def __init__(self, f):
        self.f = f

    def __str__(self):
        f = self.f
        if f>=0:
            sign=""
        else:
            sign="-"


        i = int(f)
        r = int(abs(f-i) * 10000)
        out = "%s%s.%04d"%(sign,abs(i),r)
        padding = len("xxxxx.xxxx") - len(out)
        assert padding >=0
        out = " "*padding + out
        
        return out


def writeAtom(atom):
    """atom->return the string that represents the atom in an sdfile"""
    return atom._line.rstrip()

##    symbol = atom.symbol
##    try:
##        x,y,z = map(SDFloat, (atom.x, atom.y, atom.z))
##    except AttributeError:
##        x,y,z = map(SDFloat,(0.0, 0.0, 0.0))

##    weight = atom.weight
##    charge = atom.charge
##    stereo = 0
##    hcount = atom.hcount

##    massdiff = 0 # XXX FIX ME
##                 # don't know how to compute yet

##    if charge:
##        charge = 4-charge

##    hcount = atom.hcount + 1
    
##    return "%s%s%s%3s%2i%3i%3i%3i%3i%3i"%\
##           (x,y,z,symbol,massdiff,charge,0,hcount,0,0)

def writeBond(bond, a1, a2):
    return "%3i%3i%3i%s"%(a1,
                           a2,
                           bond.bondorder,
                           bond._line.rstrip())

def writeMolecule(m, atoms=None, bonds=None, fields=1, index=1):
    if atoms is None:
        atoms = m.atoms
    if bonds is None:
        bonds = m.bonds
        
    # first write a header
    result = [""]
    result.append("  -ISIS-  03150009252D  ")
    result.append("")

    numatoms = len(atoms)
    numbonds = len(bonds)
    result.append("%3i%3i  0  0  0  0  0  0  0  0999 V2000"%(numatoms, numbonds))

    atomOutput = []

    for atom in atoms:
        result.append(writeAtom(atom))
        
    for bond in bonds:
        atom1, atom2 = bond.atoms
        a1, a2 = atoms.index(atom1), atoms.index(atom2)
        result.append(writeBond(bond, a1+1, a2+1))

    charges = []
    for atom in atoms:
        if atom.charge != 0:
            charges.append("%4s%4s"%(atoms.index(atom)+1, atom.charge))
            
    if charges:
        result.append("M  CHG%3s%s"%(len(charges), "".join(charges)))
        
    result.append("M  END")

    if fields and hasattr(m, "fields") and not m.fields.has_key("ID"):
        result.append(">  <ID> (%s)"%index)
        result.append(m.name)
        result.append("")
    
    if fields and hasattr(m, "fields"):
        for key, val in m.fields.items():
            result.append(">  <%s> (%s)"%(key, index))
            result.append(val)
            result.append("")
    else:
        result.append("")
    
    result.append("$$$$")
    return "\n".join(result)

    
if __name__ == "__main__":
    from frowns import MDL
    reader = MDL.mdlin(open("../test/bad.sdf"))
    m = reader.next()
    #print writeMolecule(m)

        
    print writeMolecule(m, m.atoms[0:1], [])
