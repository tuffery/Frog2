import vflib

def to_argedit(molecule):
    """molecule -> vflib editable graph"""
    G = vflib.ARGEdit()
    lookup = {}
    for atom in molecule.atoms:
        lookup[atom.handle] = G.InsertNode(atom)

    for bond in molecule.bonds:
        atom1, atom2 = bond.atoms
        index1, index2 = lookup[atom1.handle], lookup[atom2.handle]
        G.InsertEdge(index1, index2, bond)
        # we have to make the graph UNDIRECTED
        G.InsertEdge(index2, index1, bond)
    return G

class Path:
    def __init__(self, atoms, bonds):
        """contains atoms and bonds for a single match
        of a pattern against a molecule.
        path.atoms conatain atoms
        path.bonds contain bonds"""
        self.atoms = atoms
        self.bonds = bonds

class PathSet:
    def __init__(self, atoms, bonds, paths):
        """Contain all the matches of a target against a molecule
        PathSet.atoms contain all the atoms hit
        PathSet.bonds contain all the bonds hit
        PathSet.paths contain all the paths hit"""
        self.atoms = atoms
        self.bonds = bonds
        self.paths = paths

    def __getitem__(self, index):
        return self.paths[index]

    def __len__(self):
        return len(self.paths)
    
class Matcher:
    def __init__(self, pattern):
        """(argedit, undirected=1)->initialize the matching object
        argedit MUST be a vflib.ARGEdit object
        undirected specifies whether the initial graph was
        directed or undirected.  This is used for pruning the
        returned results"""
        self._pattern = pattern
        ag = self._argedit = to_argedit(pattern)
        self._g = vflib.GraphMatcher(ag)

    def dump(self):
        return self._pattern.dump()
    
    def match(self, molecule, limit=-1):
        """(h, limit=-1)-> return all matches against the
        target molecule

        
        The match function returns up to limit matches of the
        matcher topology on the matchableGraph topology.
        This returns all matches found in all permutations.
        """
        #molecule._graph = h = molecule._graph or to_argedit(molecule)
        h = to_argedit(molecule)   
        results = self._g.matchVF2Mono(h, limit)
        # we need to prune multiple edge results
        # from undirected graphs

        pruned = []
        atoms = []
        bonds = []
        for nodes, edges in results:
            found = {}
            prunedEdges = []
            for edge in edges:
                if not found.has_key(edge.handle):
                    prunedEdges.append(edge)
                found[edge.handle] = 1
            pruned.append(Path(nodes, tuple(prunedEdges)))
            atoms += list(nodes)
            bonds += prunedEdges
            
        results = pruned
        if results:
            return PathSet(atoms, bonds, results)

        return None

    def umatch(self, molecule, limit=-1):
        """(h, limit=-1)-> return all the unique matches against
        the molecule.
        
        The match function returns up to limit matches of the
        matcher topology on the matchableGraph topology.
        This returns all the unique matches on a graph.
        Permutations of the same match are not returned.
        """
        molecule._graph = h = molecule._graph or to_argedit(molecule)

        results = self.match(h, limit)
        # gather up the unique matches
        # each object must have a unique handle attribute!
        matches = []
        groups = []
        hash = {}
        atoms = []
        bonds = []
        for nodes, edges in results:
            hash = {}
            for object in nodes + edges:
                hash[object.handle] = 1

            if hash not in matches:
                matches.append(hash)
                groups.append((nodes, edges))
                atoms += list(nodes)
                bonds += list(edges)
        if groups:        
            return PathSet(atoms, bonds, groups)
        return None


