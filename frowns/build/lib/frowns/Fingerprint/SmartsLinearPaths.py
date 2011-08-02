"""Linear Paths

Generate linear paths for a molecule.

For example, generate all linear paths up to depth 5
paths = generatePaths(molecule, maxdepth=5)

These paths can be used for a variety of cases, but we are using
them for the purposes of fingerprinting molecules.
See Fingerprint.py

Smarts molecules are a little different than standard
molecule representations.  Sometimes, Smarts molecules
have complicated logical expressions that are used to
generate the smarts atoms and bonds.

The following code is pretty much a hack.  It tries
to determine whether smarts atoms and bonds rely
on any logical operators.

It an atom or bond does not have this property, then
a path is not made

"""
# Once again we are using a breadth first search approach to walking
# a molecule.  Each linear span is converted into a string value
# this string value is used to create the fingerprints.

#
# Modify name_atom and name_bond to change how the
# hashing works.  Do you want charges? aromaticity?
# anything.

# XXX FIX ME
# A simple optimization is to cache all the names before
# the bfs walk.
# There are more optimizations for later...

ONE_STRING_PER_PATH = 1

def name_atom(atom):
    # this extracts names from smarts atoms --> not
    # really used yet
    atom = atom.__dict__.get("value", atom)
    
    if atom.aromatic:
        return atom.symbol[0].lower() + atom.symbol[1:]
    
    return atom.symbol

def name_bond(bond):
    return bond.symbol

def countRings(object):
    """count the number of rings associated with an object,
    return the number of rings and the number of aromatic rings"""
    nrings = len(object.rings)
    aromaticRings = 0
    for ring in object.rings:
        if object.aromatic:
            aromaticRings += 1
    return nrings, aromaticRings
            
def name_ring_atom(atom):
    if atom.aromatic:
        symbol = atom.symbol[0].lower() + atom.symbol[1:]
    else:
        symbol = atom.symbol

    nrings, arings = countRings(atom)

    if not nrings:
        return symbol
    else:
        return "%s<%s:%s>"%(symbol, nrings, arings)

def name_ring_bond(bond):
    nrings, arings = countRings(bond)
    if not nrings:
        return bond.symbol
    else:
        return "%s<%s:%s>"%(bond.symbol, nrings, arings)
    
def _bfswalk(atom, visitedAtoms, path, rpath, paths, depth, maxdepth,
             name_atom, name_bond, name_ring_atom, name_ring_bond,
             make_rpath=1):
    if depth > maxdepth:
        return
    
    for bond in atom.bonds:
        oatom = bond.xatom(atom)
        if not visitedAtoms.has_key(oatom.handle):
            ##path.append("%s%s"%(name_bond(bond), name_atom(oatom)))
            path.append(name_bond(bond))
            path.append(name_atom(oatom))
            p1 = "".join(path)
            if ONE_STRING_PER_PATH:
                path.reverse()
                p2 = "".join(path)
                path.reverse()
                if p1 <= p2:
                    paths[p1] = 1
            else:            
                paths["".join(path)] = 1
            
            if make_rpath:
                rpath.append("%s%s"%(name_ring_bond(bond), name_ring_atom(atom)))
                paths["".join(rpath)] = 1            


            visitedAtoms[atom.handle] = 1
            _bfswalk(oatom, visitedAtoms, path, rpath, paths, depth+1, maxdepth,
                     name_atom, name_bond, name_ring_atom, name_ring_bond, make_rpath)
            path.pop()
            del visitedAtoms[atom.handle]
            
        
def generatePaths(molecule, maxdepth=5,
                  name_atom=name_atom, name_bond=name_bond,
                  name_ring_atom=name_ring_atom, name_ring_bond=name_ring_bond,
                  make_rpath=0):
    """(molecule, maxdepth, *name_atom, *name_bond) -> linear paths
    Generate all linear paths through a molecule up to maxdepth
    change name_atom and name_bond to name the atoms and bonds
    in the molecule

    name_atom and name_bond must return a stringable value"""
    paths = {}
    for atom in molecule.atoms:
        _bfswalk(atom, {atom:1}, [name_atom(atom)], [name_ring_atom(atom)], paths, 1, maxdepth,
                 name_atom, name_bond, name_ring_atom, name_ring_bond, make_rpath)
    return paths.keys()

if __name__ == "__main__":
    from frowns import Smiles
    mol = Smiles.smilin("CCNc1cccc1")
    paths = generatePaths(mol)
    print len(paths)
    print paths
    
    ONE_STRING_PER_PATH = 1
    paths = generatePaths(mol)
    print len(paths)
    print paths
            
