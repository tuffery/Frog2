import _pysssr

atoms = 3

connections = [(0,1), (1,2), (2,0)]

res = _pysssr.sssr(atoms, connections)
print "results is", res
