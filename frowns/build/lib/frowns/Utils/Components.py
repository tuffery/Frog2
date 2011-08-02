"""components(molecule) -> [(atoms, bonds), (atoms, bonds)]
return the connected components of a molecule
"""

def _recurse(atom, visited):
    for bond in atom.bonds:
        for next in bond.atoms:
            if not visited.has_key(next):
                visited[next] = 1
                _recurse(next, visited)
    
                
def components(molecule):
    atoms = {}
    for atom in molecule.atoms:
        atoms[atom] = 1

    components = []

    while atoms:
        start = atoms.keys()[0]
        visited = {start:1}
        _recurse(start, visited)
        bonds = {}
        for atom in visited.keys():
            del atoms[atom]
            for bond in atom.bonds:
                bonds[bond] = 1
                
        components.append((visited.keys(), bonds.keys()))

    return components

if __name__ == "__main__":
    from frowns import Smiles
    mol = Smiles.smilin("CCC.NNN")
    print components(mol)
        
    print mol.cansmiles()
