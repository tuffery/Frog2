"""Linear Paths

Generate linear paths for a molecule.

For example, generate all linear paths up to depth 5
paths = generatePaths(molecule, maxdepth=5)

These paths can be used for a variety of cases, but we are using
them for the purposes of fingerprinting molecules.
See Fingerprint.py
"""
from Fingerprint import *

# Once again we are using a depth first search approach to walking
# a molecule.  Each linear span is converted into a string value
# this string value is used to create the fingerprints.

#
# Modify name_atom and name_bond to change how the
# hashing works.  Do you want charges? aromaticity?
# anything.

# XXX FIX ME
# A simple optimization is to cache all the names before
# the dfs walk.
# There are more optimizations for later...
def name_atom(atom):
    if atom.aromatic:
        if atom.symbol == "N" and atom.imp_hcount == 0 and atom.hcount == 1:
            return "nH"
        else:
            return atom.symbol[0].lower() + atom.symbol[1:]
    return atom.symbol

def name_bond(bond, lookup={1:'-',2:'=',3:'#',4:'~'}):
    return lookup[bond.bondtype]
    
def _dfswalk(atom, visitedAtoms, path, paths, depth, maxdepth,
             name_atom, name_bond):
    if depth >= maxdepth:
        return

    for bond in atom.bonds:
        oatom = bond.xatom(atom)
        if not visitedAtoms.has_key(oatom.handle):
            path.append("%s%s"%(name_bond(bond), name_atom(oatom)))
            # only keep the path if the first character in the head of the path
            # is less than the last character in the end of the path
            if path[0][-1] <= path[-1][-1]:
                p = (depth+1, "".join(path))
                paths[p] = 1                
            visitedAtoms[atom.handle] = 1
            _dfswalk(oatom, visitedAtoms, path, paths, depth+1, maxdepth,
                     name_atom, name_bond)
            path.pop()
            del visitedAtoms[atom.handle]

def generatePaths(molecule, maxdepth=5,
                  name_atom=name_atom, name_bond=name_bond):
    """(molecule, maxdepth, *name_atom, *name_bond) -> linear paths
    Generate all linear paths through a molecule up to maxdepth
    change name_atom and name_bond to name the atoms and bonds
    in the molecule

    name_atom and name_bond must return a stringable value"""
    paths = {}
    for atom in molecule.atoms:
        _dfswalk(atom, {atom:1}, [name_atom(atom)], paths, 1, maxdepth,
                 name_atom, name_bond)
    return paths.keys()

class SplitFingerprintGenerator:
    def __init__(self, maxdepth=7, integersPerAtoms=[4]*6):
        self.maxdepth = maxdepth
        self.integersPerAtoms = integersPerAtoms
        
        assert maxdepth-1 == len(integersPerAtoms)


    def createFP(self, molecule):
        p = SplitFingerprint(self.maxdepth, self.integersPerAtoms)

        paths = generatePaths(molecule, maxdepth=self.maxdepth)
        paths.sort()
        for length, s in paths:
            p.addPath(length, s)

        return p

    
    
if __name__ == "__main__":
    from frowns import Smiles
    mol = Smiles.smilin("CCCc1cc[nH]c1")
    mol2 = Smiles.smilin("c1cc[nH]c1")

    paths = generatePaths(mol)
    pathLengths = {}
    for p in paths:
        l, s = p
        pathLengths[l] = pathLengths.get(l, []) + [s]

    generator = SplitFingerprintGenerator()
    sp = generator.createFP(mol)
    sp2 = generator.createFP(mol2)

    assert sp in sp
    assert sp2 in sp

    
            
    print "".join(map(str,sp.to_list()))
