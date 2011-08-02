from Graph import Graph
from Node import Node
from Edge import Edge

def topology_to_network(nodes, edges):
    """(nodes, edges) -> Graph
    Given a simple topology description, generate a graph network
    of the computers and cable types created by the network
    
     network_nodes is a list of computer addresses
     network cables is a list of tuples of the form
      (cabletype, computer1, computer2)"""
    
    network = Graph()
    nodeLookup = {}
    for name in nodes:
        node = Node(name)
        nodeLookup[name] = node
        network.add_node(node)
        
    for cabletype, computer1, computer2 in edges:
        node1 = nodeLookup[computer1]
        node2 = nodeLookup[computer2]
        edge = Edge(cabletype)
        network.add_edge(edge, node1, node2)

    return network

