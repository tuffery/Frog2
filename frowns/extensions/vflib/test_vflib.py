from vflib import *

q = ARGEdit()
a = "A"
b = "B"

n1 = q.InsertNode(a)
n2 = q.InsertNode(b)
q.InsertEdge(n1, n2, None)


g = ARGEdit()
n1 = g.InsertNode(a)
n2 = g.InsertNode(b)
n3 = g.InsertNode(a)
print 1
g.InsertEdge(n1, n2, None)
g.InsertEdge(n2, n1, None)
print 2
g.InsertEdge(n2, n3, None)
g.InsertEdge(n3, n2, None)

match = GraphMatcher(q)
print match.matchVF2Mono(g, -1)


# test the exception handler
class F:
    def __eq__(self, o):
        print "comparing to o"
        print "this should fail"
        a = lafkjasdlfjasdf

q2 = ARGEdit()
n1 = q2.InsertNode(F())

match = GraphMatcher(q2)
res = match.matchVF2Mono(g, -1)
print res
for nodes, edges in res:
    print nodes, edges
