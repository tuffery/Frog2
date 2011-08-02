import time
from frowns import Smiles
from frowns import Smarts
mol = Smiles.smilin("CCCCC(C)CCCC(C)CCCCCC(C)C(C)CCCCCC(C)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCN")
pat = Smarts.compile("C*")
print len(pat.match(mol))

pat2 = Smarts.compile("[#6]")
t1 = time.time()
print len(pat2.match(mol))
t2 = time.time()
print t2-t1
t1 = time.time()
print len(pat2.match(mol))
t2 = time.time()
print t2-t1

print mol.vfgraph
