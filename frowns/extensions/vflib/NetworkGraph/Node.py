from GraphObject import GraphNode

class Node(GraphNode):
    def __init__(self, machine):
        GraphNode.__init__(self)
        self.machine = machine

    def __eq__(self, node):
        # all nodes match all other nodes
        print "*"*44
        print node.machine[0], self.machine[0]
        return node.machine[0] == self.machine[0]
        #return 1

    def __repr__(self):
        return "%s(%s)"%(self.__class__.__name__,
                         `self.machine`)
    def clone(self):
        return Node(self.machine)
