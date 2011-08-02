import vflib
import random, types

"""vflib  positive testing
 This test checks the Monomorphism substructure functions.
 A random core G is built
 Then this core is randomly extended to become H

 G should still be contained in H (using monomorphism)

 This should also give a good idea how to build your own
 graphs.
"""

class Core:
    """Core -> stores the core of a graph"""
    def __init__(self, nodes=None, edges=None, undirected=0):
        self.nodes = nodes or []
        self.edges = edges or {}
        self.undirected = undirected # is the graph directed
                                     #  or undirected

    def InsertNode(self, data):
        """(data) -> node index
        Add a node to a graph returning the index of the new node"""
        self.nodes.append(data)
        return len(self.nodes)

    def InsertEdge(self, start, end, data):
        """(start, end, data) -> Add an edge to a graph from start to end.
        The edge has type data"""
        assert start >= 0 and start < len(self.nodes), "start node error %s (max %s)"%(start, len(self.nodes))
        assert end >= 0 and end < len(self.nodes), "end node error %s (max %s)"%(end, len(self.nodes))
        assert type(start) == types.IntType
        assert type(end) == types.IntType
        self.edges[(start, end)] = data
        if self.undirected:
            self.edges[(end, start)] = data

    def to_graph(self):
        """()->ARGEdit()
        return an ARGEdit object that can be used as a match graph
        or a target graph"""
        G = vflib.ARGEdit()
        for node in self.nodes:
            G.InsertNode(node)

        for key, data in self.edges.items():
            start, end = key
            G.InsertEdge(start, end, data)

        return G

    def clone(self):
        """()->return a clone of this core"""
        edges = {}
        edges.update(self.edges)
        return Core(self.nodes[:], edges)

    def to_connection_table(self):
        """return a connection table graph"""
        result = []
        result.append("%s %s"%(len(self.nodes), len(self.edges)))
        for node in range(len(self.nodes)):
            result.append("N%s"%(node+1))
        output = {}
        for edge in self.edges:
            node1, node2 = edge

            if self.undirected and output.has_key((node2, node1)):
                continue

            output[(node1, node2)] = 1
            result.append("generic %s %s"%(node1+1, node2+1))
        return "\n".join(result)
        

    def to_graphviz(self):
        """return a dotty compatible graph"""
        result = ["graph Test {"]
        output = {}
        for edge in self.edges:
            node1, node2 = edge
            if self.undirected and output.has_key((node2, node1)):
                continue

            output[(node1, node2)] = 1
            result.append("\tN%s -- N%s"%(node1+1, node2+1))
        result.append("}")
        return "\n".join(result)

def grow_core(maxnodes = 100, maxedges = 1000, probNew=0.5, undirected=1):
    G = Core(undirected=undirected)

    G.InsertNode(None)
    G.InsertNode(None)
    G.InsertEdge(0, 1, None)
    numnodes = 2
    numedges = 1
    
    # add up to maxedges to the graph
    # if the generated probability is below probNew
    # then add a new node and connect it to a random
    # node from the graph
    # otherwise randomly connect two nodes in the graph
    for edge in range(int(random.random() * maxedges)):
        p = random.random()
        if p < probNew:
            # add a new node and connect
            G.InsertNode(None)
            numnodes += 1
            connection = random.randrange(numnodes-1)
            G.InsertEdge(numnodes-1, connection, None)
        else:
            edge1 = edge2 = 0
            count = 0
            while edge1 == edge2:
                edge1 = random.randrange(numnodes)
                edge2 = random.randrange(numnodes)
                count += 1
                if count > 10000: raise "Too many iterations in while loop!"
            G.InsertEdge(edge1, edge2, None)

    return G
            
def make_core(maxnodes=100, maxedges=1000):
    G = Core()
    maxnodes = max(3, int(random.random() * maxnodes)+1)
    for i in range(maxnodes):
        G.InsertNode(i)
    lastnode = i
    
    for j in range(int(random.random() * maxedges)+1):
        n1 = n2 = 0
        count = 0
        while n1 == n2:
            n1 = int(random.random() * lastnode)
            n2 = int(random.random() * lastnode)
            count += 1
            if count > 10000: raise "Too many iterations in while loop!"
             
        G.InsertEdge(n1, n2, 1)
        
    return G

def add_to_core(G, maxnodes=100, maxedges=1000):
    start = len(G.nodes)
    maxnodes = max(3, int(random.random() * maxnodes)+1)
    for i in range(maxnodes):
        G.InsertNode(i+start)
    lastnode = len(G.nodes)
    
    for j in range(int(random.random() * lastnode)+1):
        n1 = n2 = 0
        count = 0
        while n1 == n2:
            n1 = int(random.random() * lastnode)
            n2 = int(random.random() * lastnode)
            count += 1
            if count > 10000: raise "Too many iterations in while loop!"
            G.InsertEdge(n1, n2, 1)

def test(num_cores=5, num_extensions=5):
    for i in range(num_cores):
        G = make_core()

        graph = G.to_graph()
        matcher = vflib.GraphMatcher(graph)
        assert matcher.matchVF(graph, -1) != []

        for i in range(num_extensions):
            H = G.clone()
            assert H.nodes == G.nodes
            assert H.edges == G.edges
            add_to_core(H)
            matchers = ["matchVFMono", "matchVF2Mono"]
            graph2 = H.to_graph()
            for functionName in matchers:
                func = getattr(matcher, functionName)
                assert func(graph2, -1) != [], functionName
    print "random extension tests passed"

def test_reference_count_bug():
    graph = make_core()
    G = graph.to_graph()
    H = graph.to_graph()

    matcher = vflib.GraphMatcher(G)
    matcher.matchVF(H, -1)

    del G
    try:
        matcher.matchVF(H, -1)
    except RuntimeError:
        print "Caught runtime error from deallocated reference"
    else:
        print "no deallocated reference!"

if __name__ == "__main__":
    test_reference_count_bug()
    test()

    # grow a random core
    #g = grow_core()
    #print g.to_connection_table()

    
