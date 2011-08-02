"""
Smallest set of smallest rings detection code inspired by
Pat Walter's and Matt Stahl's Babel rings.c with plenty of
changes!

I don't like it very much but it appears to work all right.

Figueras is slightly faster but has some problems.
This implementation forms a spanning tree up to a given depth
(20) which means that 40+ atom rings might not be found
correctly.  see maxDepth below.

implementation Notes:
 This could be sped up slightly by making a class based
  object out of it.

 The one problem you may run into is the bond._closure property
  that may or may not be defined.  Essentially this indicates
  whether a bond was used as a closure bond.  This ring detection
  algorithm can leverage off the smiles parser to start ring
  searches only from closure bonds.

  It's kind of ugly and won't work with standard daylight
  molecules because they don't have this property

  But my guess is look at GraphGenerator and somewhere just set
  closureBond._closure = 1


 BUGS:
 The following code won't work directly on daylight (PyDaylight)
 molecules yet for the following reasons.
 
   atom.index needs to be the atom's index into mol.atoms
   bond.index needs to be the bond's index into mol.bonds
   atom._inring is used to indicate whether an atom is
                in a ring
   bond._inring is used to indicate whether an bond is
                in a ring

   bond._closure optionally indicates whether the bond was
     used in the parsing routine as a closure.  This can
     reduce the number of atoms searched as ring starts.

   These data place holders will be refactored in a future
   version of this code for better (and easier) checking
   against daylight molecules.  That is, it would be nice
   if you could use either a daylight or one of our own
   molecules for testing purposes.  Although we still
   have to worry about what to do with the ring information
   for the daylight molecule case...
                

The license for the rings.c was use at your own risk and
as-is.  I guess that applies here :)
"""
from frowns import Cycle
from frerejaque import frerejaqueNumber

# ----------------------------------------------------------------------
#  Mark all the atoms and bonds that are in rings
# ----------------------------------------------------------------------
def _recurseMarkRingAtomsAndBonds(graph, atom, atomsVisited, bondsVisited, path):
    if atomsVisited.has_key(atom.index):
        # full circle, complete the path
        j = len(path)-1
        path[j]._inring = 1
        j -= 1
        while (j >= 0):
            bond = path[j]; j-=1
            bond._inring = 1
            shouldBreak = 0
            for matom in bond.atoms:
                matom._inring = 1
                if matom.index == atom.index:
                    shouldBreak = 1
            if shouldBreak:
                break
    else:
        atomsVisited[atom.index] = 1
        for bond in atom.bonds:
            if not bondsVisited.has_key(bond.index):
                bondsVisited[bond.index] = 1
                nextatom = bond.xatom(atom)
                assert nextatom is not None
                path.append(bond)
                _recurseMarkRingAtomsAndBonds(graph, nextatom, atomsVisited,
                                              bondsVisited, path)

def markRingAtomsAndBonds(graph):
    """(graph)->mark all the atoms that can possibly be in a ring
    system.  These are marked with atom._inring = 1 and bond._inring = 1
    otherwise atom._inring = 0 and bond._inring  = 0"""
    path = []
    bondsVisited = {}
    atomsVisited = {}

    for atom in graph.atoms:
        if not atomsVisited.has_key(atom.index):
            _recurseMarkRingAtomsAndBonds(graph, atom, atomsVisited,
                                          bondsVisited, path)

# ----------------------------------------------------------------------
#  Finds all linear paths from a root atom to all other atoms
#   This is the time consuming part
# ----------------------------------------------------------------------
class Span:
    """A span is really a path back to a central atom.  Consider it
    a type of linked list that can compute a linear path from
    an atom to a root atom.

    span = Span(atom, parent)
    path = []
    span.update(path)
    makes path a linear path from span.atom to the root span
    """
    def __init__(self, atom, parent=None):
        self.atom = atom
        self.parent = parent
        path = self.path = [atom]
        if parent:
            path += parent.path

    def update(self, path):
        path.append(self.atom)
        if self.parent: self.parent.update(path)

    def __repr__(self):
        return "Span(%s)"%self.path

def breadthFirstSpans(root, spans, visited, maxDepth=20):
    """(atom, spans, visited, maxDepth=20)->generate the breadthFirst spans
    eminating from atom to all other accessible atoms.
    spans must be a list of None's equal to the number of atoms
    in the molecule so that spans[atom.index] should indicate a
    particular atom's span to the root

    The maxDepth flag is fairly important, but we're dealing with
    small molecules here right?  You can always set this to
    the number of atoms in the molecule for a complete solution.
    """
    spans[root.index] = Span(root)
    visited = visited.copy()
    boundary = {root.index:root}
    depth = 0
    while 1:
        nextBoundary = {}
        for atom in boundary.values():
            for oatom in atom.oatoms:
                # only go to atoms that can be in rings
                # you might want to remove oatom._inring because
                # I'm not sure if it works yet
                if oatom._inring and not visited.has_key(oatom.index):
                    nextBoundary[oatom.index] = oatom
                    visited[oatom.index] = oatom
                    spans[oatom.index] = Span(oatom, spans[atom.index])
        if not nextBoundary:
            break

        boundary = nextBoundary
        depth += 1
        if depth > maxDepth:
            break
# ----------------------------------------------------------------------
#  From a given root bond make all the spans from its atoms
#   then recombine them to form rings
# ----------------------------------------------------------------------
def gotOne(lring, rring, sets):
    """(ring, sets)->add ring to the sets dictionary"""
    # we want the atoms to be in order of traversal
    # so we need to reverse one path and then add it
    # to the other path
    ring = lring[:]
    ring.reverse()
    ring.extend(rring)

    # use the indices to get a unique identifier for the path
    indices = [atom.index for atom in ring]
    indices.sort()
    lookup = tuple(indices)
    if not sets.has_key(lookup):
        sets[lookup] = ring
        
def getRing(graph, bond, sets):
    """(graph, bond, sets)->find all rings starting from bond and
    add them to the sets member"""
    # the terminology left and right is arbitrary here
    #  just so's you know
    spans_left = [None] * len(graph.atoms)
    spans_right = [None] * len(graph.atoms)
    latom, ratom = bond.atoms

    a1, a2 = bond.atoms

    breadthFirstSpans(latom, spans_left, {ratom.index:ratom, latom.index:latom})
    breadthFirstSpans(ratom, spans_right,  {ratom.index:ratom, latom.index:latom})

    for lspan in spans_left:
        if not lspan:
            continue
        
        # can we close this path from the other side?
        rspan = spans_right[lspan.atom.index]
        if rspan:
            lpath, rpath = lspan.path, rspan.path

            pathintersects = 0
            # remember the first atom in lpath or rpath
            #  should be the root atom!
            assert lpath[0] == rpath[0] #sanity check
            
            lring = [lpath[0]]             # left ring
            for latom in lpath[1:]:
                lring.append(latom)
                
                rring = []                 # right ring
                for ratom in rpath[1:]:
                    if latom.index == ratom.index:
                        if len(lring) + len(rring) > 2:
                            gotOne(lring, rring, sets)
                        pathintersects = 1
                        break
                    rring.append(ratom)

                    #
                    # findbond is called a LOT
                    # this is a relatively slow function so
                    # SPEED IT UP!  (I did locally with great results)
                    if ratom.findbond(latom) and len(lring)+len(rring) > 2:
                        gotOne(lring, rring, sets)
                        
                # if paths intersect at the same atom
                #  we can't continue on this path (it won't be a smallest ring)
                if pathintersects:
                    break

# ----------------------------------------------------------------------
#  Now we have to remove all large rings that completely
#  contain smaller rings
# ----------------------------------------------------------------------
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
    while rings and (len(sets) > numRings):
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

    # -------------------------------------------------------
    #  old code (no fusion ring detection)
    #  uncomment to see differences
    ##while rings and (len(sets) > numRings):
    ##    checkSize, checkPath, checkRing, checkKey = rings.pop()
    ##    smallerRings = 0L
    ##    for size, path, ring, key in rings:
    ##        if size <= checkSize:
    ##            smallerRings |= path
    ##    if (smallerRings & checkPath) == checkPath:
    ##        del sets[checkKey]
        
    return sets

# ----------------------------------------------------------------------
#  Actual Ring routine
# ----------------------------------------------------------------------
from frowns import Cycle
def getBondsInOrder(atoms):
    bonds = []
    for a1, a2 in zip(atoms, atoms[1:] + [atoms[0]]):
        bond = a1.findbond(a2)
        assert bond
        bonds.append(bond)
    return bonds

def sssr(graph, closuresMarked=1):
    """(molecule, closuresMarked=1) -> determine the sssr (smallest set of smallest rings) for a graph
    if closuresMarked is set then it is assumed that bonds used in closures
    are marked with an attribute _closure"""

    # yuck, atoms are being modified with attributes
    # this has to be replaced, perhaps with a dictionary
    # lookup.  This won't work with daylight atoms and
    # molecules which is a bug
    for atom in graph.atoms:
        atom._inring = 0
        atom.rings = []

    for bond in graph.bonds:
        bond._inring = 0
        bond.rings = []
        
    frj = frerejaqueNumber(graph)

    if not frj:
        graph.cycles = []
        graph.rings = []
        return graph
 
    markRingAtomsAndBonds(graph)
    
    sets = {}
    closures = []

    # we can speed this up slightly by starting only from atoms
    # that are deemed to be closures during the parsing routines
    if closuresMarked:    
        for bond in graph.bonds:        
            if hasattr(bond, "_closure"):
                closures.append(bond)
    else:
        # otherwise add all bonds that are connected
        # to atoms that have more than one bond
        for bond in graph.bonds:
            for atom in bond.atoms:
                if len(atom.bonds) > 1:
                    closures.append(bond)
                    break


    for bond in closures:
        getRing(graph, bond, sets)

    if len(sets) > frj:
        uniqueifyRings(sets, frj)

    # what do you do if len(sets) > frj?  This can be
    # true for cubane for instance!

    graph.rings = sets.values()
    index = 0
    cycles = graph.cycles = []
    for ring in graph.rings:
        bonds = getBondsInOrder(ring)
        cycles.append(Cycle.Cycle(ring, bonds))        

    return graph

def sssr_no_closures(graph):
    return sssr(graph, 0)
