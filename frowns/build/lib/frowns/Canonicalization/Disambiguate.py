"""Disambiguate

FreedDisambiguate(graph)

given atom and bond equiv_classes use Freed's technique to
find symmettry classes and symmetry orders for the graph

These were initial described in the following papers published in
JCICS
   1. WEININGER D,
     "SMILES, A chemical language and information-system. 1.
      Introductions to methodology and encoding rules"

      JOURNAL OF CHEMICAL INFORMATION AND COMPUTER SCIENCES,28(#1),
      1988,31-36 
   2. WEININGER D, WEININGER A, WEININGER JL,
     "SMILES 2. Algorithm for generation of unique SMILES notation"

     JOURNAL OF CHEMICAL INFORMATION AND COMPUTER SCIENCES,29(#2),
     1989,97-101 

They do not talk about canonicalizing using bond types or stereo
chemistry though.  See notes below for details.
"""
import Primes

class FreedDisambiguate:
    def __init__(self, graph):
        """(graph)
        
        given atom and bond equiv_classes use Freed's technique to
        find symmettry classes and symmetry orders for the graph
        """
        self.graph = graph
        atoms = graph.atoms
        self.range = range(len(atoms))
        offsets = self.offsets = []

        # symmetry classes are the symmetrical points in a graph
        #  so if two vertices share the same value, they are symmetric
        symclasses = self.symclasses = []

        # symmetry orders describe the traversal order to canonicalize
        # the graph.  All of these should be unique values
        symorders = self.symorders = [0] * len(atoms)
        indices = {}

        # cache the atom indices in a lookup table
        index = 0
        for atom in atoms:
            indices[atom.handle] = index
            index += 1

        # initialize the necessary data
        for atom in atoms:
            # For each atom, get the neighbors offsets and the connected
            # bondtype for a given atom
            # N.B. the oatom (adjacent atom)
            # list and the bond list must be kept in sync
            # bondtype must be an integer!!!!
            
            # the weininger algorithm is modified by also
            # multiplying by not only the neighboring symorder
            #  but by the bondtype used to reach the neighbor.
            #  The traversal algorithm will then traverse to
            #  the lowest symorder atom with the lowest bondtype
            #  first.
            neighbors = zip([indices[x.handle] for x in atom.oatoms],
                            [x.bondtype for x in atom.bonds])

            # store the symorder index of the neighboring atoms
            # and the bondtype used to traverse to each neighbor
            
            # This forms a list of the form
            # (neighbor_symorder_index, bondtype),
            #   ...(neighbor_symorder_index, bondtype) ...
            offsets.append(neighbors)

            # the symmetry classes start out as the atom's equivalence
            # class
            symclasses.append(atom.equiv_class)

        #########################################################
        # find the symmetry classes
        # rank the current symmetry orders 
        self.rank()

        # Find the sym classes invariant for the current
        self.symclasses = self.findInvariant(self.symclasses)

        # now find the symmetry orders from the symmetry classes
        self.symorders = self.symclasses[:]
        
        self.findInvariantPartitioning()

        # Give them back to the atoms and bonds
        symclasses = self.symclasses
        symorders = self.symorders
        i = 0
        for atom in atoms:
            atom.symclass = symclasses[i]
            atom.symorder = symorders[i]
            i = i + 1

    def disambiguate(self, symclasses):
        """Use the connection to the atoms around a given vertex
        as a multiplication function to disambiguate a vertex"""
        offsets = self.offsets
        result = symclasses[:]
	for index in self.range:
	    try:
		val = 1
		for offset, bondtype in offsets[index]:
		    val *= symclasses[offset] * bondtype
	    except OverflowError:
                # Hmm, how often does this occur?
		val = 1L
		for offset, bondtype in offsets[index]:
		    val *= symclasses[offset] * bondtype
            result[index] = val
	return result

    def rank(self):
        """convert a list of integers so that the lowest integer
        is 0, the next lowest is 1 ...
        note: modifies list in place"""
        # XXX FIX ME, should the lowest value be 1 or 0?
        symclasses = self.symclasses
        stableSort = map(None, symclasses, range(len(symclasses)))
        stableSort.sort()

        last = None
        x = -1
        for order, i in stableSort:
            if order != last:
                x += 1
                last = order
            symclasses[i] = x

    def breakRankTies(self, oldsym, newsym):
        """break Ties to form a new list with the same integer ordering
        from high to low

        Example
        old = [ 4, 2, 4, 7, 8]  (Two ties, 4 and 4)
        new = [60, 2 61,90,99]
        res = [ 4, 0, 3, 1, 2]
                *     *        This tie is broken in this case
        """
        stableSort = map(None, oldsym, newsym, range(len(oldsym)))
        stableSort.sort()

        lastOld, lastNew = None, None
        x = -1
        for old, new, index in stableSort:
            if old != lastOld:
                x += 1
                # the last old value was changed, so update both
                lastOld = old
                lastNew = new
            elif new != lastNew:
                # break the tie based on the new info (update lastNew)
                x += 1
                lastNew = new
            newsym[index] = x
    
    def findLowest(self, symorders):
        """Find the position of the first lowest tie in a
        symorder or -1 if there are no ties"""
        _range = range(len(symorders))
        stableSymorders = map(None, symorders, _range)

        # XXX FIX ME
        # Do I need to sort?
        stableSymorders.sort()
        lowest = None
        for index in _range:
            if stableSymorders[index][0] == lowest:
                return stableSymorders[index-1][1]
            lowest = stableSymorders[index][0]
        return -1

    def findInvariant(self, symclasses):
        """(symclasses) -> converge the disambiguity function
        until we have an invariant"""
        get = Primes.primes.get
        disambiguate = self.disambiguate
        breakRankTies = self.breakRankTies

        while 1:
            newSyms = map(get, symclasses)
            newSyms = disambiguate(newSyms)
            breakRankTies(symclasses, newSyms)
            if symclasses == newSyms:
                return newSyms
            symclasses = newSyms

    def findInvariantPartitioning(self):
        """Keep the initial ordering of the symmetry orders
        but make all values unique.  For example, if there are
        two symmetry orders equal to 0, convert them to 0 and 1
        and add 1 to the remaining orders

          [0, 1, 0, 1]
        should become
          [0, 2, 1, 3]"""
        
        symorders = self.symorders[:]
        _range = range(len(symorders))
        while 1:
            pos = self.findLowest(symorders)
            if pos == -1:
                self.symorders = symorders
                return
            for i in _range:
                symorders[i] = symorders[i] * 2 + 1
            symorders[pos] = symorders[pos] - 1
            
            symorders = self.findInvariant(symorders)


