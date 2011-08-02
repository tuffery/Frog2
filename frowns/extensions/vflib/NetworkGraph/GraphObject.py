class HandleGenerator:
    def __init__(self, start=0):
        self.start = start

    def next(self):
        self.start += 1
        return self.start

# base class for all nodes and edges.  Just provides a unique
# handle for each object 
class GraphObject:
    def __init__(self, label="",
                 handlegenerator=HandleGenerator()):
        """() -> graph object base clase"""
        self.handle = handlegenerator.next()
        self.parent = None
        self.label = label

    def set_parent(self, parent):
        if parent is None:
            self.parent = None
        else:
            assert self.parent is None, "%s already belongs to a graph!"%\
                   (self.__class__.__name__)
            self.parent = parent

    def __eq__(self, other):
        """Graph objects are equivalent if the labels
        are equivalent"""
        return self.label == other.label
        
    def __hash__(self):
        return self.handle
    
    def __repr__(self):
        return "%s(%s)"%(self.__class__.__name__, `self.label`)

    def clone(self):
        return GraphObject(label=self.label)
    

class GraphNode(GraphObject):
    def findedge(self, otherNode):
        return self.parent._findedges[self.handle].get(otherNode.handle, None)
        
    def __getattr__(self, key):
        if key == "edges":
            return self.parent._atoms[self.handle]
        raise AttributeError, key

    def findedge(self, otherNode):
        """(otherNode)->returns the edge between this node and the otherNode)
        None if there is no edge"""
        return self.parent._findedges[self.handle].get(otherNode.handle, None)

    def clone(self):
        return GraphNode(self.label)
    
    def destroy(self):
        """remove this node from the parent graph"""
        self.parent.remove_node(self)
 
class GraphEdge(GraphObject):
    def xnode(self, node):
        """(node)->return the node at the other end of this edge
        or None if node is not part of this edge"""
        return self.parent._xnodes[self.handle].get(node.handle, None)

    def __getattr__(self, key):
        if key == "nodes":
            return self.parent._edges[self.handle]

        raise AttributeError, key

    def clone(self):
        return GraphEdge(self.label)
    
    def destroy(self):
        """remove this edge from the parent graph"""
        self.parent.remove_edge(self)        
    
