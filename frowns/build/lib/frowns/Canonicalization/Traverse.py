"""Recursively traverse a molecule building up a canonical
representation.

Each atom of a Molecule or Graph must have a attribute 'symorder' which
is a unique number.  This number guarantees only one traversal for
the graph.

Additionally each bond must have an attribute equiv_class which is
a unique value for each different type of bond.  This guarantees
proper canonicalization of bonds as well as atoms.

canonical_string = Traverse.draw(molecule)

Advanced usage:
canonical_string = Traverse.draw(molecule, TraversalType)

TraversalType controls how the traversal is represented.
SmilesTraversal is the default TraversalType but this can be
subclassed for different representations.  For example it is
easy to create a subclass to generate Tripos Line formats.
"""

from SmilesTraversal import SmilesTraversal, SmartsTraversal, \
     IsomericSmilesTraversal

def _traverse(atom, traverse, prevAtom,
              visitedAtoms, visitedBonds,
              atoms, bonds, Traversal, bondIndex=0):
    visitedAtoms[atom] = 1
    traverse.addAtom(atom)
    atoms.append(atom)

    bondsToTraverse = []
    traversals = []
    
    for bond in atom.bonds:
        oatom = bond.xatom(atom)
        if prevAtom is not None and oatom == prevAtom:
            # we are traversing back the way we came!
            # so don't...
            pass
        elif visitedAtoms.has_key(oatom):
            # a closure!
            traverse.addClosure(atom, oatom, bond)
            bonds.append(bond)
            visitedBonds[bond] = 1
        else:
            bondsToTraverse.append((oatom.symorder,
                                    bond.equiv_class,
                                    bondIndex,
                                    oatom,
                                    bond))
        bondIndex += 1

    if not bondsToTraverse:
        # dead end, return
        return

    bondsToTraverse.sort()

    for symorder, bondEclass, index, oatom, obond in bondsToTraverse:
        if visitedAtoms.has_key(oatom):
            # somehow, we've seen this atom so skip it
            continue
        nextTraverse = Traversal(traverse)
        traversals.append(nextTraverse)
        nextTraverse.addBond(obond)
        bonds.append(obond)
        visitedBonds[obond] = 1

        _traverse(oatom, nextTraverse, atom,
                  visitedAtoms, visitedBonds,
                  atoms, bonds, Traversal, bondIndex)

    for t in traversals[:-1]:
        traverse.addBranch()
        traverse.append(t)
        traverse.addBranchEnd()
        
    for t in traversals[-1:]:
        traverse.append(t)

def _get_lowest_symorder(atoms):
    best = atoms[0]
    for atom in atoms[1:]:
        if atom.symorder < best.symorder:
            best = atom
    return best

def draw(molecule, TraversalType=SmilesTraversal):
    """(molecule)->canonical representation of a molecule
    Well, it's only canonical if the atom symorders are
    canonical, otherwise it's arbitrary.

    atoms must have a symorder attribute
    bonds must have a equiv_class attribute"""
    result = []
    atoms = allAtoms = molecule.atoms

    visitedAtoms = {}
    #
    # Traverse all components of the graph to form
    # the output string
    while atoms:
        atom = _get_lowest_symorder(atoms)
        visitedAtoms[atom] = 1

        visitedBonds = {}
        nextTraverse = TraversalType()
        atomsUsed, bondsUsed = [], []
        _traverse(atom, nextTraverse, None,
                  visitedAtoms, visitedBonds,
                  atomsUsed, bondsUsed, TraversalType)
        atoms = []
        for atom in allAtoms:
            if not visitedAtoms.has_key(atom):
                atoms.append(atom)
        assert nextTraverse.atoms == atomsUsed
        assert nextTraverse.bonds == bondsUsed, "%s %s"%(
            nextTraverse.bonds, bondsUsed)
        

        result.append((str(nextTraverse),
                       atomsUsed, bondsUsed))

    result.sort()
    fragments = []
    for r in result:
        fragments.append(r[0])

    return ".".join(fragments), result

def drawSmarts(molecule):
    return draw(molecule, TraversalType=SmartsTraversal)

def drawIsomeric(molecule):
    return draw(molecule, TraversalType=IsomericSmilesTraversal)
