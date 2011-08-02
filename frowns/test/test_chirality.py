from frowns import Smiles

m = Smiles.smilin("N[C@](C)(F)C(=O)O")
print "initial smiles:", m.arbsmiles(isomeric=1)
print "canonical smiles:", m.cansmiles(isomeric=1)

print
m = Smiles.smilin("[C@H](C)(F)(O)")
print "initial smiles:", m.arbsmiles(isomeric=1)
print "canonical smiles:", m.cansmiles(isomeric=1)
