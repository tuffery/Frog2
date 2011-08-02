import Components

def computeMW(mol):
    """mol -> mw of largest fragment, mw of total molecule"""
    fragments = []
    total = 0.0
    if not mol.atoms: return 0.0, 0.0
    for atoms, bonds in Components.components(mol):
        mw = 0.0
        for atom in atoms:
            mw += atom.mass + atom.hcount
            fragments.append(mw)
        total += mw

    largest = max(fragments)
    return largest, total

if __name__ == "__main__":
    from frowns import Smiles
    mol = Smiles.smilin("CCC.NNN.Br")
        
    print computeMW(mol)
