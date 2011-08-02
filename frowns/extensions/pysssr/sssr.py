import _pysssr
_sssr = _pysssr.sssr
from frowns.Cycle import Cycle
from frowns.perception import figueras

class PySSSRError(Exception): pass

def toposort(initialAtoms):
    """topologically sort the atoms and return the bonds as well.
    Assumes that the atoms are the smallest set of rings"""

    indices = {}
    for atom in initialAtoms:
        indices[atom.index] = 1

    
    next = initialAtoms[0]
    atoms = [next]
    bonds = []
    del indices[next.index]
    
    while len(atoms) < len(initialAtoms):
        for neighbor in next.oatoms:
            if neighbor.index in indices:
                del indices[neighbor.index]
                bond = next.findbond(neighbor)
                bonds.append(bond)
                atoms.append(neighbor)                
                next = neighbor
                break
        else:
            raise PySSSRError("Toposort Error: bad atoms for ring")
        
    bonds.append(atoms[-1].findbond(atoms[0]))
    return atoms, bonds
        
def sssr(molecule):
    connections = []
    for bond in molecule.bonds:
        a1, a2 = bond.atoms
        connections.append((a1.index, a2.index))

    try:
        indices = _sssr(len(molecule.atoms), connections)
    except:
        return figueras.sssr(molecule)
    rings = []
    cycles = []
    for r in indices:
        atoms = [molecule.atoms[x] for x in r]
        atoms, bonds = toposort(atoms)        
        rings.append((atoms, bonds))
        cycles.append(Cycle(atoms, bonds))

    molecule.rings = rings
    molecule.cycles = cycles
    return molecule

    
if __name__ == "__main__":
    from frowns import Smiles
    import time
    
    data = open("../Working/frowns/test/data/NCI_small").read()
    count = 0
    i = 0
    for line in data.split("\n"):
        i += 1
        smiles = line.split()[0]
        m = Smiles.smilin(smiles)
        t1 = time.time()
        m = sssr(m)
        t2 = time.time()
        count += (t2-t1)
        
        assert len(m.cycles) == len(m.rings), "failed %s in %f average"%(float(count)/i)

    print "average sssr", float(count)/i


