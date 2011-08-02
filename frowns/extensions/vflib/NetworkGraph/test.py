from NetworkGraph.GraphObject import GraphObject
N = GraphObject()
print N.handle
M = GraphObject()
print M.handle
print N == M
print N is M
print N.handle == M.handle
from NetworkGraph.Graph import Graph
g = Graph()
from NetworkGraph.GraphObject import GraphNode
node1 = GraphNode()
node2 = GraphNode(label="blue")
print node1 == node2
g.add_node(node1)
g.add_node(node2)
print g.has_node(node1)
from NetworkGraph.GraphObject import GraphEdge
edge1 = GraphEdge(label="my dog has fleas")
g.add_edge(edge1, node1, node2)
print g.has_edge(edge1)
print edge1.nodes
n = edge1.xnode(node1)
print n is node2
print n.handle == node2.handle
g.dump()
g.remove_edge(edge1)
g.dump()
g.add_edge(edge1, node1, node2)

h = g.to_graph()
matcher = g.to_matcher()
results = matcher.match(h)
for nodes, edges in results:
    print nodes
    print edges
clone = g.clone()
for node in g.nodes:
    assert not clone.has_node(node)
for original, cloned in zip(g.nodes, clone.nodes):
  assert original == cloned
  assert original is not cloned
node3 = GraphNode("I am a clone!")
edge2 = GraphEdge("new edge")
clone.add_node(node3)
n1 = clone.nodes[0]
clone.add_edge(edge2, node3, n1)
matchableClone = clone.to_graph()
results = matcher.umatch(matchableClone)
nodes, edges = results[0]
partialClone = clone.clone(ignoreNodes=nodes, ignoreEdges=edges)
partialClone.dump()
