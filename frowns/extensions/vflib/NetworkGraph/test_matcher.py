from topology import *

# Example code
network_nodes = [
    'darkstar.wi.mit.edu',
    'fang.harvard.edu',
    'www.wi.mit.edu'
    ]

network_cables = [ ('t1', 'darkstar.wi.mit.edu', 'fang.harvard.edu'),
                   ('t3', 'fang.harvard.edu', 'www.wi.mit.edu'),
                   ('10base2', 'www.wi.mit.edu', 'darkstar.wi.mit.edu')]

# network is the original graph
# the example will match the network against it self
# and remove the matched nodes from the original network.
network = topology_to_network(network_nodes, network_cables)

#####################################################
# make a couple of clones to play with later
network1 = network.clone()

# This is an independed graph, the nodes and edges are different!
# they are not the same objects!
for node in network.nodes:
    assert network1.has_node(node) == 0

# although they are equivalent to each other
for node1, node2 in zip(network.nodes, network1.nodes):
    assert node1 == node2

#####################################################
# dump the toplogy
network.dump()

print
# make the matching graph
matcher = network.to_matcher()

# generate a graph to be matched against
graph = network.to_graph()

#####################################################
# collect all the results
match =  matcher.match(graph)
print "number of results", len(match)
for nodes, edges in match:
    for edge in edges:
        print '\t', edge, edge.nodes
    print

#####################################################
# now collect the unique results
umatch = matcher.umatch(graph)
nodes, edges = umatch[0]
print "number of unique results", len(umatch)

#####################################################
# remove the matched edges

partialClone = None
for nodes, edges in umatch:
    for edge in edges:
        print '\t', edge, edge.nodes
    print

    # try the partial cloning procedure
    # we'll just do this for the first match
    if partialClone is None:
        partialClone = network.clone(ignoreNodes=nodes,
                                     ignoreEdges=edges)

#####################################################
# let's check the nodes and edges of the graph
# the network topology should be empty
print "Should not be modified"
network.dump()
print

print "partial clone should be empty"
partialClone.dump()
print

print "original clone should be the original network"
network1.dump()
