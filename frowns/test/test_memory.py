from frowns import Smiles

class T:
    def __del__(self):
        print "deleting"

t = T()
del t

import weakref
class P:
    def __init__(self, t):
        self.t = weakref.ref(t)

t = T()
p = P(t)
del t
        
# watch the memory this process consumes to see if there is a memory leak
#  anywhere
# This, obvious, should not be used in standard regression
# testing :)
while 1:
    mol = Smiles.smilin("CCNCc1ccccc1")
