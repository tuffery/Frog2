"""A new ring detection algorithm inspired by talking
to various people"""
from frowns.Cycle import Cycle
from frerejaque import frerejaqueNumber
from InRings import getAtomsBondsInRings

class Path:
    """A span is really a path back to a central atom.  Consider it
    a type of linked list that can compute a linear path from
    an atom to a root atom.

    span = Path(atom, parent)
    path = []
    span.update(path)
    makes path a linear path from span.atom to the root span
    """
    def __init__(self, atom, rings,
                 last=None, parent=None):
        self.atom = atom
        self.last = last
        self.parent = parent
        self.rings = rings
        path = self.path = [atom]
        hpath = self.hpath = [atom.handle]
        self.bits = 1L << atom.index
        
        if parent:
            path += parent.path
            hpath += parent.hpath
            self.bits |= parent.bits
            
        self.ring = parent and atom is parent.path[-1]
        # do we have a ring?
        if self.ring:
            lookup = parent.hpath
            lookup.sort()
            rings[tuple(lookup)] = parent.path
        
    def next(self):
        """bfs traverse forming new paths don't go the way we came"""
        paths = []
        rings = self.rings
        last = self.last
        for atom in self.atom.oatoms:
            if atom is last: continue
            elif atom is self.path[-1] or not atom.handle in self.hpath:
                # XXX FIX Me, we need to do this
                # to complete the ring
                paths.append(Path(atom, rings, self.atom, self))

        return paths
    
    def __repr__(self):
        return "%s Path(%s)"%(self.ring, self.path)


    
class RingFinder:
    def __init__(self, start, rings):
        self.start = start
        self.rings = rings

        self.paths = [Path(start, rings)]        
        
    def next(self):
        new = []
        for path in self.paths:
            newpaths = path.next()
            for path in newpaths:
                if not path.ring:
                    new.append(path)
        self.paths = new
        
    def len(self):
        return len(self.rings)

class Fusion:
    def __init__(self, path, rings):
        self.path = path
        self.rings = rings

    def updatePath(self, newpath, ring):
        self.path |= newpath
        self.rings.append(ring)
 
def uniqueifyRings(sets, numRings):
    rings = []
    for lookup, ring in sets.items():
        path = 0L
        for atom in ring:
            path |= 1L<<atom.index

        rings.append((len(ring), path, ring, lookup))
    rings.sort()

    # -----------------------------------------------------
    # new code enabling ring fusion detection

    # a partition will store a set of fused rings
    partitions = []
    for ring in rings:
        rsize, path, ring, key = ring
        if not partitions:
            partitions.append(Fusion(path, [ring]))
            continue

        for partition in partitions:
            # see if this path intersects one
            # of the partitions we already have
            if path & partition.path != 0:
                partition.updatePath(path, ring)
                break
        else:
            # no partition intersections were found,
            #  so add a new one
            partitions.append(Fusion(path, [ring]))

    # now from each partition, remove largest rings that are completly
    # contained within the partition
    while rings:# and (len(sets) > numRings):
        for partition in partitions:
            # hmm, can a partition not have a ring???
            if not partition.rings:
                continue
            
            checkSize, checkPath, checkRing, checkKey = rings.pop()
            smallerRings = 0L
            for size, path, ring, key in rings:
                if size <= checkSize:
                    smallerRings |= path
            if (smallerRings & checkPath) == checkPath:
                del sets[checkKey]
                break
            if not rings:
                break
        
    return sets

def sssr(graph):
    targetRings = frerejaqueNumber(graph)
    ringAtoms, ringBonds = getAtomsBondsInRings(graph)

    # prepare the generators
    generators = []
    rings = {}
    for atom in ringAtoms:
        generators.append(RingFinder(atom, rings))

    while len(rings) < targetRings and generators:
        for generator in generators:
            # advance the rings
            generator.next()
            if not generator.paths:
                generators.remove(generator)
        #print len(rings), targetRings
        #if len(rings) >= targetRings:
        #    print "uniqueifing"
        #    uniqueifyRings(rings, targetRings)
        #    print len(rings)
        if not generators:
            print "out of generators"

    #assert len(rings) == targetRings
    
    cycles = graph.cycles = []
    for atoms in rings.values():
        bonds = []
        for atom1, atom2 in zip(atoms, atoms[1:]):
            bond = atom1.findbond(atom2)
            assert bond
            bonds.append(bond)
        bond = atoms[0].findbond(atoms[-1])
        assert bond
        bonds.append(bond)
        cycles.append(Cycle(atoms, bonds))

    return graph

def testPath():
    from frowns import Smiles
    mol = Smiles.smilin("C1CCCC1CCCC2CCCC3CCC3CCC2")
    rings = {}
    spans = [Path(mol.atoms[0], rings)]
    while spans:
        new = []
        for span in spans:
            new += span.next()
            if len(rings) == 1:
                break
        if len(rings) == 1:
            break
        spans = new
    print rings

def test():
    from frowns import Smiles
    import RingDetection
    import sssr as sssr2
    import time
    from frowns.Canonicalization import Disambiguate, Traverse
    
    ##mol = Smiles.smilin("C12C3C4C1C5C4C3C25CCCCCC12C3C4C1C5C4C3C25")
    ##mol1 = Smiles.smilin("C1O[CH]1C2=CC=CC=C2",
    ##                    transforms = [])
    ##print mol1.cansmiles()
    smiles = "COC(=O)[CH]1C2CCC(CC2)[CH]1C(=O)OCC"
    smiles = "CC1(C)[CH]2CC[CH](CC(O)=O)[CH]1C2"
    smiles = "C[CH]1C(=O)O[CH]2CCN3CC=C(COC(=O)[C](C)(O)[C]1(C)O)[CH]23"
    #smiles = "CC1=CC(=O)C=CC1=O"
    smiles = "NC1=CC2=C(C=C1)C(=O)C3=C(C=CC=C3)C2=O"
    smiles = "COC1=CC2=C(C=C1OC)[C]34CCN5CC6=CCO[CH]7CC(=O)N2[CH]3[CH]7[CH]6C[CH]45.OS(O)(=O)=O"

    mol1 = Smiles.smilin(smiles,
                         transforms = [])
    
    print mol1.cansmiles()
    NUM = 1
    
    t1 = time.time()
    for i in range(NUM):
        sssr(mol1)
    t2 = time.time()
    print (t2-t1)/NUM
    symbols = []
    for atom in mol1.atoms:
        symbols.append(atom.symbol)

    for atom in mol1.atoms:
        atom.symbol += ".%s."%atom.index
        atom.symorder = mol1.atoms.index(atom)
        
    index = 0
    for cycle in mol1.cycles:
        print cycle.atoms
        print cycle.bonds
        for atom in cycle.atoms:
            atom.symbol += "<%s>"%index
        index += 1

    print Traverse.draw(mol1)[0]

    for atom, symbol in zip(mol1.atoms, symbols):
        atom.symbol = symbol

    for atom in mol1.atoms:
        atom.symbol += ".%s."%atom.index
        atom.symorder = mol1.atoms.index(atom)
        
    t1 = time.time()
    for i in range(NUM):
        RingDetection.sssr(mol1)
    t2 = time.time()
    print (t2-t1)/NUM
    index = 0
    for cycle in mol1.cycles:
        print cycle.atoms
        print cycle.bonds
        for atom in cycle.atoms:
            atom.symbol += "<%s>"%index
        index += 1

    print Traverse.draw(mol1)[0]
    
if __name__ == "__main__":
    testPath()
    test()
