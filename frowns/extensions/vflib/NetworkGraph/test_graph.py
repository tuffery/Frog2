from Graph import Graph
from Node import Node
from Edge import Edge

node1 = Node('darkstar.wi.mit.edu')
node2 = Node('fang.harvard.edu')
node3 = Node('www.wi.mit.edu')
edge1 = Edge('t3')
edge2 = Edge('t1')

network = Graph()
network.add_node(node1)
network.add_node(node2)
network.add_node(node3)
network.add_edge(edge1, node1, node2)
network.add_edge(edge2, node1, node3)

print "*"*44
print "current graph"
network.dump()

print
print "*"*44
print "now removing", edge2
network.remove_edge(edge2)

network.dump()

print "*"*44
print "now removing", node1
network.remove_node(node1)
network.dump()

