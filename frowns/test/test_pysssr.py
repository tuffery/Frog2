import sys
sys.path.append("/Users/stockwel/Working")
from frowns import Smiles
import time
from frowns.perception import RingDetection, BasicAromaticity, figueras, pysssr

transforms = [figueras.sssr,
              BasicAromaticity.aromatize]

tests = [x.split()[0] for x in open("data/NCI_small")][0:100]

t1 = time.time()

for smiles in tests:
    m = Smiles.smilin(smiles, transforms)

t2 = time.time()

time1 = (t2-t1)/len(tests)
print "start", t1
print "end", t2
print "elapsed", t2-t1
print "number of smiles", len(tests)
print "average molecule generation", time1


transforms = [pysssr.sssr,
              BasicAromaticity.aromatize]

t1 = time.time()
for smiles in tests:
    m = Smiles.smilin(smiles, transforms)

t2 = time.time()

time2 = (t2-t1)/len(tests)
print "start", t1
print "end", t2
print "elapsed", t2-t1
print "number of smiles", len(tests)
print "average molecule generation", time2

print "speed up is ", time1/time2
