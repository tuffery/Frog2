# ----------------------------------------------------------------------
#  Compute the frerejaque number for a graph
# ----------------------------------------------------------------------
class Component:
    def __init__(self, atoms, bonds):
        """(atoms, bonds) -> store the atoms and bonds that make
        a connected component"""
        self.atoms = atoms
        self.bonds = bonds
        
def getConnectedComponents(graph):
    """(graph)->return a list of connected subgraphs of a molecule"""
    atoms, bonds = {}, {}
    for atom in graph.atoms:
        atoms[atom.index] = atom
    for bond in graph.bonds:
        bonds[bond.index] = bond
        
    components = []
    while atoms:
        startAtom = atoms.values()[0]

        boundary = {startAtom.index:startAtom}
        connectedAtoms = {startAtom.index:startAtom}
        connectedBonds = {}

        while boundary:
            for atom in boundary.values():
                del boundary[atom.index]
                del atoms[atom.index]
                for oatom in atom.oatoms:
                    if not connectedAtoms.has_key(oatom.index):
                        boundary[oatom.index] = oatom
                        connectedAtoms[oatom.index] = oatom
                        for bond in oatom.bonds:
                            if bonds.has_key(bond.index):
                                connectedBonds[bond.index] = bond
                                del bonds[bond.index]

                
        components.append(Component(connectedAtoms.values(),
                                    connectedBonds.values()))

    # sanity check, make sure that we have used up all
    #  atoms and bonds
    assert not atoms
    assert not bonds

    return components

def frerejaqueNumber(graph):
    """(graph)->return the number of rings in a graph"""
    # step 1, get all the components from a graph
    components = getConnectedComponents(graph)
    result = 0
    for subgraph in components:
        result += len(subgraph.bonds) - len(subgraph.atoms) + 1

    return result

def test():
    from frowns import Smiles
    import time
    smiles, id = "OC1=C(Cl)C=CC=C1[N+]C2CCC3CC3CCC2", 1
    mol = Smiles.smilin(smiles, transforms=[])
    assert frerejaqueNumber(mol) == 3
    t1 = time.time()
    for i in range(100):
        frerejaqueNumber(mol)
    print (time.time() - t1)/100    

if __name__ == "__main__":
    test()
