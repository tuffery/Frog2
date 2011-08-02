import vflib

class Wrapper:
    def __init__(self, data, nodeid):
        self.data = data
        self.nodeid = nodeid

    def __eq__(self, other):
        return self.data == other.data

    def __repr__(self):
        return "Wrapper(data=%s,nodeid=%s)"%(self.data,self.nodeid)


G = vflib.ARGEdit()
G.InsertNode(Wrapper(1, nodeid=0))
G.InsertNode(Wrapper(2, nodeid=1))
G.InsertEdge(0,1, "first edge")

H = vflib.ARGEdit()
H.InsertNode(Wrapper(2, nodeid=2))
H.InsertNode(Wrapper(1, nodeid=1))
H.InsertEdge(1,0, "first edge")

matcher = vflib.GraphMatcher(G)
print matcher.matchVF(H, -1)
