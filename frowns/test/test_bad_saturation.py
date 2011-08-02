from frowns.MDL import mdlin
from frowns.Utils import SaturatedRings
from frowns.perception import RingDetection, BasicAromaticity, sssr
#transforms = [RingDetection.sssr, BasicAromaticity.aromatize]

file = open("data/really_bad_saturation.mol")
reader = mdlin(file)
m = reader.next()

print len(m.cycles)
#file = open("data/bad_saturation_molecule.mol")
#reader = mdlin(file)
#m2 = reader.next()

s = SaturatedRings.getSaturationScore(m)
print "saturation", s
