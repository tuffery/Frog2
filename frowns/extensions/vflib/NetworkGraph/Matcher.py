import vflib

class Matcher:
    def __init__(self, argedit, undirected=1):
        """(argedit, undirected=1)->initialize the matching object
        argedit MUST be a vflib.ARGEdit object
        undirected specifies whether the initial graph was
        directed or undirected.  This is used for pruning the
        returned results"""
        self._g = vflib.GraphMatcher(argedit)
        self.undirected = undirected

    def match(self, h, limit=-1):
        """(h, limit=-1)-> return all matches against h
        h MUST be a vflib.ARGEdit object
        
        The match function returns up to limit matches of the
        matcher topology on the matchableGraph topology.
        This returns all matches found in all permutations.
        """
        results = self._g.matchVF2Mono(h, limit)
        # we need to prune multiple edge results
        # from undirected graphs
        if self.undirected:
            pruned = []
            for nodes, edges in results:
                found = {}
                prunedEdges = []
                for edge in edges:
                    if not found.has_key(edge.handle):
                        prunedEdges.append(edge)
                    found[edge.handle] = 1
                pruned.append((nodes, tuple(prunedEdges)))
            results = pruned
        return results
                        
        return results

    def umatch(self, h, limit=-1):
        """(h, limit=-1)-> return all the unique matches against h
        h MUST be a vflib ARGEDit object
        
        The match function returns up to limit matches of the
        matcher topology on the matchableGraph topology.
        This returns all the unique matches on a graph.
        Permutations of the same match are not returned.
        """
        results = self.match(h, limit)
        # gather up the unique matches
        # each object must have a unique handle attribute!
        matches = []
        groups = []
        hash = {}
        for nodes, edges in results:
            hash = {}
            for object in nodes + edges:
                hash[object.handle] = 1

            if hash not in matches:
                matches.append(hash)
                groups.append((nodes, edges))

        return groups



        
