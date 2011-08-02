"""Routines to find all atoms and bonds that are
in rings"""
# ----------------------------------------------------------------------
#  Mark all the atoms and bonds that are in rings
# ----------------------------------------------------------------------
def _recurseMarkRingAtomsAndBonds(atom, atomsVisited, bondsVisited,
                                  path, atomsInRings, bondsInRings):
    if atomsVisited.has_key(atom.handle):
        # full circle, complete the path
        j = len(path)-1
        bond = path[j]
        bondsInRings[bond.handle] = bond

        j -= 1
        while (j >= 0):
            bond = path[j]; j-=1
            bondsInRings[bond.handle] = bond

            shouldBreak = 0
            for matom in bond.atoms:
                atomsInRings[matom.handle] = matom

                if matom.handle == atom.handle:
                    shouldBreak = 1
            if shouldBreak:
                break
    else:
        atomsVisited[atom.handle] = 1
        for bond in atom.bonds:
            if not bondsVisited.has_key(bond.handle):
                bondsVisited[bond.handle] = 1
                nextatom = bond.xatom(atom)
                assert nextatom is not None
                path.append(bond)
                _recurseMarkRingAtomsAndBonds(nextatom, atomsVisited,
                                              bondsVisited, path,
                                              atomsInRings, bondsInRings)

def getAtomsBondsInRings(graph):
    """(graph)->mark all the atoms that can possibly be in a ring
    system.  These are marked with atom._inring = 1 and bond._inring = 1
    otherwise atom._inring = 0 and bond._inring  = 0
    atoms and bonds must have an _inring attribute"""
    path = []
    bondsVisited = {}
    atomsVisited = {}
    atomsInRings = {}
    bondsInRings = {}
    

    for atom in graph.atoms:
        if not atomsVisited.has_key(atom.handle):
            _recurseMarkRingAtomsAndBonds(atom, atomsVisited,
                                          bondsVisited, path,
                                          atomsInRings, bondsInRings)

    return atomsInRings.values(), bondsInRings.values()


def test():
    from frowns import Smiles
    import time
    smiles, id = "OC1=C(Cl)C=CC=C1[N+]C2CCC3CC3CCC2", 1
    mol = Smiles.smilin(smiles, transforms=[])

    t1 = time.time()
    for i in range(100):
        atoms, bonds = getAtomsBondsInRings(mol)
    print (time.time() - t1)/100


if __name__ == "__main__":
    test()
