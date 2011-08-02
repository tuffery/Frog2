"""Ring Perception Algorithm ala figueras

This code has been tested against a slower version of ring detection
similar to the ones that the Babel converter uses.

This code is MUCH faster although their appears to be deficencies
in this ring detection code that need investigating.

On the plus side, it generates the ringsets that can be used for
aromaticity detection very quickly (tested on 300k NCI compounds,
the canonicalization using aromaticity detection was the same for
both ring detectors).

On the minus side, there may be some bugs in the code that may limit
its use as a full smallest set of smallest rings finder.

This is also the ring detection code that CDK uses for their ring
perception so that might be a branch to investigate.
"""
from frowns.Cycle import Cycle
#from CheckFiguerasRings import checkRings

def sssr(molecule):
    """molecule -> generate the molecule.cycles that contain
    the smallest set of smallest rings"""
    results = {}
    lookup = {}
    fullSet = {}
    oatoms = {}
    # XXX FIX ME
    # copy atom.oatoms to atom._oatoms
    # atom._oatoms will be modified my the routine

    for atom in molecule.atoms:
        atom.rings = []
	fullSet[atom.handle] = 1
	lookup[atom.handle] = atom
        oatoms[atom.handle] = atom.oatoms[:]

    for bond in molecule.bonds:
        bond.rings = []
        
    trimSet = []

    while fullSet:
	nodesN2 = []
        minimum, minimum_degree = None, 100000


	# find the N2 atoms and remove atoms with degree 0
	for atomID in fullSet.keys():
	    atom = lookup[atomID]
	    degree = len(oatoms[atom.handle])
	    if degree == 0:
		del fullSet[atomID]
		#fullSet.remove(atomID)
	    elif degree == 2:
		nodesN2.append(atom)

	    # keep track of the minimum degree
	    if (degree > 0) and ( (not minimum) or 
		                (degree < minimum_degree)):
		minimum, minimum_degree = atom, degree

	if not minimum:
	    # nothing to do!  (i.e. can't have a ring)
	    break

	if minimum_degree == 1:
	    # these cannot be in rings so trim and remove
	    # my version of trimming
            for oatom in oatoms[minimum.handle]:
		oatoms[oatom.handle].remove(minimum)
            oatoms[minimum.handle] = []
	    del fullSet[minimum.handle]

	elif minimum_degree == 2:
	    # find the rings!
	    startNodes = []
	    for atom in nodesN2:
		ring, bonds = getRing(atom, fullSet, lookup, oatoms)

		if ring:
                    rlookup = ring[:]
                    rlookup.sort()
                    rlookup = tuple(rlookup)
		    if (not results.has_key(rlookup)):# not in results):
			results[rlookup] = ring, bonds
			startNodes.append(atom)

	    # in case we didn't get a ring remove the head of the nodesN2
	    startNodes = startNodes or [nodesN2[0]]
	    for atom in startNodes:
		# again, my version of trimming
                if oatoms[atom.handle]:
		    oatom = oatoms[atom.handle].pop()
		    oatoms[oatom.handle].remove(atom)

	elif minimum_degree > 2:
	    # no N2 nodes so remove the "optimum" edge to create
	    # N2 nodes in the next go-around.
	    ring, bonds = getRing(minimum, fullSet, lookup, oatoms)
            if ring:
                key = ring[:]
                key.sort()
                key = tuple(key)
                if not results.has_key(key):
                    results[key] = ring, bonds
                    atoms = map(lookup.get, ring)
                    atoms, bonds = toposort(atoms, bonds)
                    checkEdges(atoms, lookup, oatoms)
            else:
                del fullSet[minimum.handle]
 	else:
	    raise ShouldntGetHereError

    # assign the ring index to the atom
    rings = []
    index = 0

    # transform the handles back to atoms
    for result, bonds in results.values():
	ring = []
        for atomID in result:
	    atom = lookup[atomID]
            assert atom.handle == atomID
	    ring.append(atom)
	rings.append((ring, bonds))
        index = index + 1

    molecule.rings = rings
    potentialCycles = []
    index = 0
    for atoms, bonds in rings:
        # due to the dictionaries used in getRing
        # the atoms are not in the order found
        # we need to topologically sort these
        # for the cycle
        atoms, bonds = toposort(atoms, bonds)
        potentialCycles.append((atoms, bonds))

    rings = potentialCycles#checkRings(potentialCycles)
    molecule.rings = rings
    molecule.cycles = [Cycle(atoms, bonds) for atoms, bonds in rings]
    return molecule

def toposort(initialAtoms, initialBonds):
    """initialAtoms, initialBonds -> atoms, bonds
    Given the list of atoms and bonds in a ring
    return the topologically sorted atoms and bonds.
    That is each atom is connected to the following atom
    and each bond is connected to the following bond in
    the following manner
    a1 - b1 - a2 - b2 - ... """
    atoms = []
    a_append = atoms.append
    bonds = []
    b_append = bonds.append

    # for the atom and bond hashes
    # we ignore the first atom since we
    # would have deleted it from the hash anyway
    ahash = {}
    bhash = {}
    for atom in initialAtoms[1:]:
        ahash[atom.handle] = 1
        
    for bond in initialBonds:
        bhash[bond.handle] = bond

    next = initialAtoms[0]
    a_append(next)

    # do until all the atoms are gone
    while ahash:
        # traverse to all the connected atoms
        for atom in next.oatoms:
            # both the bond and the atom have to be
            # in our list of atoms and bonds to use
            # ugg, nested if's...  There has to be a
            # better control structure
            if ahash.has_key(atom.handle):
                bond = next.findbond(atom)
                assert bond
                # but wait! the bond has to be in our
                # list of bonds we can use!
                if bhash.has_key(bond.handle):
                    a_append(atom)
                    b_append(bond)
                    del ahash[atom.handle]
                    next = atom
                    break
        else:
            raise RingException("Atoms are not in ring")

    assert len(initialAtoms) == len(atoms)
    assert len(bonds) == len(atoms) - 1
    lastBond = atoms[0].findbond(atoms[-1])
    assert lastBond
    b_append(lastBond)
    return atoms, bonds
    
def getRing(startAtom, atomSet, lookup, oatoms):
    """getRing(startAtom, atomSet, lookup, oatoms)->atoms, bonds
    starting at startAtom do a bfs traversal through the atoms
    in atomSet and return the smallest ring found

    returns (), () on failure
    note: atoms and bonds are not returned in traversal order"""

    path = {}
    bpaths = {}
    for atomID in atomSet.keys():
	# initially the paths are empty
	path[atomID] = None
        bpaths[atomID] = []
    
    q = []
    handle = startAtom.handle
    for atom in oatoms[handle]:
	q.append((atom, handle))
        path[atom.handle] = {atom.handle:1, handle:1}
        bpaths[atom.handle] = [startAtom.findbond(atom)]
            
    qIndex = 0
    lenQ = len(q)

    while qIndex < lenQ:	
	current, sourceHandle = q[qIndex]
        handle = current.handle
	qIndex += 1

        for next in oatoms[handle]:
	    m = next.handle

	    if m != sourceHandle:
		if not atomSet.has_key(m):
		    return (), ()
                
		if path.get(m, None):
		    intersections = 0
		    for atom in path[handle].keys():
			if path[m].has_key(atom):
			    intersections = intersections + 1
			    sharedAtom = atom

		    if intersections == 1:
			del path[handle][sharedAtom]
			path[handle].update(path[m])
			result = path[handle].keys()
                        bond = next.findbond(current)
                        # assert bond not in bpaths[handle] and bond not in bpaths[m]
                        bonds = bpaths[handle] + bpaths[m] + [bond]
			return result, bonds
		else:
		    path[m] = path[handle].copy()
		    path[m][m] = 1
                    bond = next.findbond(current)
                    # assert bond not in bpaths[m] and bond not in bpaths[handle]
                    bpaths[m] = bpaths[handle] + [next.findbond(current)]
		    q.append((next, handle))
		    lenQ = lenQ + 1

    return (), ()


def checkEdges(ringSet, lookup, oatoms):
    """atoms, lookup -> ring
    atoms must be in the order of traversal around a ring!
    break an optimal non N2 node and return the largest ring
    found
    """
    bondedAtoms = map( None, ringSet[:-1], ringSet[1:] )
    bondedAtoms += [ (ringSet[-1], ringSet[0]) ]

    # form a lookup for the ringSet list
    atomSet = {}
    for atomID in ringSet:
	atomSet[atomID] = 1
    results = []
    
    # for each bond in the ring, break it and find the smallest
    # rings starting on either side of the bond
    # keep the largest but rememeber to add the bond back at the
    # end
    for atom1, atom2 in bondedAtoms:
        # break a single edge in the ring
        handle1 = atom1.handle
        handle2 = atom2.handle
        oatoms1 = oatoms[handle1]
        oatoms2 = oatoms[handle2]
        index1 = oatoms1.index(atom2)
        index2 = oatoms2.index(atom1)

        # break the bond
        del oatoms1[index1]
	del oatoms2[index2]

	ring1 = getRing(atom1, atomSet, lookup, oatoms)
	ring2 = getRing(atom2, atomSet, lookup, oatoms)
	
	# keep the larger of the two rings
	if len(ring1) > len(ring2):
	    results.append((len(ring1),
                            handle1, handle2,
                            ring1))
	else:
    	    results.append((len(ring2),
                            handle2, handle1,
                            ring2))

        # retie the bond
        oatoms1.insert(index1, atom2)
        oatoms2.insert(index2, atom1)


    if not results:
	return None

    #  find the smallest ring
    size, incidentHandle, adjacentHandle, smallestRing = min(results)
    # dereference the handles
    incident, adjacent = lookup[incidentHandle], lookup[adjacentHandle]

    # break the bond between the incident and adjacent atoms
    oatomsI = oatoms[incidentHandle]
    oatomsA = oatoms[adjacentHandle]
    assert incident in oatomsA
    assert adjacent in oatomsI
    
    oatomsI.remove(adjacent)
    oatomsA.remove(incident)
