import vflib

# every node matches every other node
# each node has an id so you know which node you
# are dealing with.
class Node:
    def __init__(self, id):
        self.id = id

    def __eq__(self, other):
        return 1
    
    def __repr__(self):
        return "Node(%s)"%repr(self.id)


# make the matching graph
G = vflib.ARGEdit()
gnodes = [Node(1), Node(2), Node(3), Node(4)]
gedges = [(0,1), (1,2), (2,3)]
# Graph is:
#  1-2-3-4

for node in gnodes:
    G.InsertNode(node)

for source, dest in gedges:
    G.InsertEdge(source, dest, None)

####################################################
# make the target graph
# the target graph is undirected so there should
# be two matches returned by the matching function
H = vflib.ARGEdit()
# Graph is
#  A-B-C-D
hnodes = [Node('A'), Node('B'), Node('C'), Node('D')]
hedges = [(0,1), (1,2), (2,3)]

for node in hnodes:
    H.InsertNode(node)

for source, dest in hedges:
    # need to add both edges for an undirected graph
    H.InsertEdge(source, dest, None)
    H.InsertEdge(dest, source, None)

###################################################
# Match G against H and show all the matches
matcher = vflib.GraphMatcher(G)
result =  matcher.matchVF2Mono(H, -1)
index = 1
for nodes, edges in result:
    print "match", index
    index += 1
    for gnode, hnode in zip(gnodes, nodes):
        print "\t", gnode, "->", hnode

# output should be
#match 0
#	Node(1) -> Node(A)
#	Node(2) -> Node(B)
#	Node(3) -> Node(C)
#	Node(4) -> Node(D)
#match 1
#	Node(1) -> Node(D)
#	Node(2) -> Node(C)
#	Node(3) -> Node(B)
#	Node(4) -> Node(A)
