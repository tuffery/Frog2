"""Graph
   a class for graph manipulation

"""
import vflib
from Matcher import Matcher

# for now these molecules will be immutable
def is_remove(list, object):
    """remove objects from a list using is equivalence instead
    of comparison equivalence"""
    result = []
    for e in list:
        if object is not e:
            result.append(e)
    return result

class Graph:
    def __init__(self, nodes=None, edges=None, undirected=1):                
        """(undirected=1)->initialize a graph structure
        undirected=1 makes an undirected graph
        undirected=0 makes a directed graph"""
        self.nodes = []
        self.edges = []
        self.undirected = undirected
        
        self._nodes = {} # hash the nodes to their edges
        self._edges = {} # hash the edges to their nodes
        
        self._findedges = {}
        self._xnodes = {}
        self.dirty = 1
        self.name = ""


    def add_node(self, node):
        """(node)-> add a node instance to the graph"""
        if self._nodes.has_key(node.handle):
            raise "Already have node"

        # set the parent of the node to be this object
        # set up all the relevant datastructures for
        # the node
        node.set_parent(self)
        node.index = len(self.nodes)
        self.nodes.append(node)
        self._nodes[node.handle] = []
        self._findedges[node.handle] = {}
        self.dirty = 1

    def add_edge(self, edge, node1, node2):
        """(edge, node1, node2) -> add an edge between node1 and node2.
        Edge must be an edge instance and node1 and node2 must be node
        instances"""
        if not self._nodes.has_key(node1.handle) or \
           not self._nodes.has_key(node2.handle):
            raise "Nodes not in molecule, need to add first"

        edge.set_parent(self)
        edge.index = len(self.edges)
        self.edges.append(edge)
        self._edges[edge.handle] = (node1, node2)
        self._xnodes[edge.handle] = {node1.handle:node2,
                                     node2.handle:node1}
        self._findedges[node1.handle][node2.handle] = edge
        self._findedges[node2.handle][node1.handle] = edge

        self._nodes[node1.handle].append(edge)
        self._nodes[node2.handle].append(edge)
        self.dirty = 1

    def clone(self, ignoreNodes=None, ignoreEdges=None):
        """->return a clone of this graph
        Note that the actual node and edge objects
        are not cloned!  So affecting the node properties of
        one graph might change the node of the cloned
        graph.  The topology's are distinct though."""
        ignoreNodes = ignoreNodes or []
        ignoreEdges = ignoreEdges or []
        g = Graph(undirected=self.undirected)
        oldNodeToNewNode = {}

        # clone the edges and add to the new graph
        # keep a mapping from old node to cloned node
        for node in self.nodes:
            # skip nodes to be ignored
            if node in ignoreNodes:
                continue
            
            clonedNode = node.clone()
            oldNodeToNewNode[node] = clonedNode
            g.add_node(clonedNode)

        for edge in self.edges:
            # skip edges to be ignored
            if edge in ignoreEdges:
                continue
            
            node1, node2 = self._edges[edge.handle]
            if node1 in ignoreNodes or node2 in ignoreNodes:
                continue
            
            clonedEdge = edge.clone()

            newNode1, newNode2 = oldNodeToNewNode[node1], \
                                 oldNodeToNewNode[node2]
            g.add_edge(clonedEdge, newNode1, newNode2)

        return g

    def dump(self):
        """Print out the nodes topology in human readable form"""
        print "Topology"
        print "Nodes"
        for node in self.nodes:
            print "\t",node

        print
        print "Edges"
        for edge in self.edges:
            print "\t", edge, edge.nodes
        
    def has_node(self, node):
        """return 1 if the graph has the component node"""
        return self._nodes.has_key(node.handle)

    def has_edge(self, edge):
        """return 1 if the graph has the component edge"""
        return self._edges.has_key(edge.handle)
    
    def remove_edge(self, edge):
        """(edge)->remove edge from the graph"""
        # if we don't have the edge, return
        #  XXX FIX ME should we raise an exception???
        if not self.has_edge(edge):
            return
        # we need to use the is operator here since we
        # are overloading __eq__ for the edges and nodes
        self.edges = is_remove(self.edges, edge)

        a1, a2 = self._edges[edge.handle]
        del self._edges[edge.handle]
        del self._xnodes[edge.handle]

        handle1, handle2 = a1.handle, a2.handle
        self._nodes[handle1] = is_remove(self._nodes[handle1], edge)
        self._nodes[handle2] = is_remove(self._nodes[handle2], edge)
        self.dirty = 1
        edge.set_parent(None)

    def remove_node(self, node):
        """(node)->remove node from the graph"""
        # if we don't have the node, return
        #  XXX FIX ME should we raise an exception?
        if not self.has_node(node):
            return
        for edge in self._nodes[node.handle]:
            self.remove_edge(edge)

        self.nodes = is_remove(self.nodes, node)
        del self._nodes[node.handle]
        self.dirty = 1
        node.set_parent(None)


    def to_graph(self):
        """->return an ARGEdit object from this graph"""
        G = vflib.ARGEdit()
        lookup = {}
        for node in self.nodes:
            print "adding node", node
            lookup[node.handle] = G.InsertNode(node)

        for edge in self.edges:
            node1, node2 = self._edges[edge.handle]
            index1, index2 = lookup[node1.handle], lookup[node2.handle]
            G.InsertEdge(index1, index2, edge)
            if self.undirected:
                G.InsertEdge(index2, index1, edge)

        return G

    def to_matcher(self):
        """->create a matching graph from this graph"""
        return Matcher(self.to_graph(), undirected=self.undirected)
            

