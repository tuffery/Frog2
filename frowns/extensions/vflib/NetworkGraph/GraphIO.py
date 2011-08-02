"""Graph related input and output

 Connection table format
 numnodes, numedges
 machine_address
 .
 .
 cable_type node_index1 node_index2
 .
 .

There should be as many nodes as numnodes and edges as numedges.

An example of a connection table is:
3 2
darkstar.wi.mit.edu
fang.harvard.edu
www.wi.mit.edu
t3 1 2
t2 2 3

Note that the node indices start with 1 and not 0
"""
from Graph import Graph
from Node import Node
from Edge import Edge

def read_connection_table(file):
    """file->graph"""
    line = file.readline()
    if not line:
        return None

    try:
        numnodes, numedges = map(int, line.split())
    except ValueError:
        raise AssertionError("Can't parse header into integer values")
    except IndexError:
        raise AssertionError("Header has wrong number of elements")

    result = Graph()
    nodes = []
    # read the nodes
    for i in range(numnodes):
        line = file.readline()
        if not line:
            raise AssertionError("Need %s nodes"%numnodes)

        # just get the machine name
        machine = line.strip()
        node = Node(machine)
        result.add_node(node)
        nodes.append(node)

    # read the edges
    edges = []
    for i in range(numedges):
        line = file.readline()
        if not line:
            raise AssertionError("Need %s edges!"%numedges)

        line = line.strip()
        groups = line.split()
        if len(groups) < 3:
            raise AssertionError("Format for edge is cable node_index1 node_index2\ngot %s"%
                                 line)

        if len(groups) > 3:
            raise AssertionError("Too many columns for edge")
        

        cable, index1, index2 = groups
        try:
            index1 = int(index1)
        except ValueError:
            raise AssertionError("%s is not an integer for edge %s"%(index1, line))
        try:
            index2 = int(index2)
        except ValueError:
            raise AssertionError("%s is not an integer for edge %s"%(index2, line))

        node1 = nodes[index1-1]
        node2 = nodes[index2-1]
        
        edge = Edge(cable)
        result.add_edge(edge, node1, node2)

    return result

def graph_to_connection_table(graph):
    result = [ "%s %s"%(len(graph.nodes), len(graph.edges))]
    # we need to make a lookup table of the nodes
    # just incase their __eq__ is overloaded and
    # list.index operations don't work
    indices = {}
    index = 0
    for node in graph.nodes:
        result.append(node.machine)
        indices[node] = index
        index += 1

    
    for edge in graph.edges:
        node1, node2 = edge.nodes
        index1 = indices[node1]
        index2 = indices[node2]
        result.append("%s %s %s"%(edge.cable, index1, index2))

    return "\n".join(result)

# This data file has 3 nodes and 2 edges
test = """3 2
darkstar.wi.mit.edu
fang.harvard.edu
www.wi.mit.edu
t3 1 2
t2 2 3
"""

if __name__ == "__main__":
    from StringIO import StringIO
    # make a fake file
    file = StringIO(test)
    graph = read_connection_table(file)
    # dump it out
    graph.dump()

    # make the connection table
    table = graph_to_connection_table(graph)
    # make it into a fake file
    newfile = StringIO(table)
    graph2 = read_connection_table(newfile)
    graph2.dump()
    
