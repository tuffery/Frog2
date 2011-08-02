from GraphObject import GraphEdge

class Edge(GraphEdge):
    def __init__(self, cable):
        GraphEdge.__init__(self)
        self.cable = cable

    def __eq__(self, edge):
        return 1
    

    def __repr__(self):
        return "%s(%s)"%(self.__class__.__name__,
                         `self.cable`)

    def clone(self):
        return Edge(self.cable)
    
