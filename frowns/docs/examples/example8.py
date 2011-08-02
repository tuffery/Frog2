from frowns import Smiles

mol = Smiles.smilin("c1ccccc1CCC1CC1")

index = 0
for cycle in mol.cycles:
    print "cycle", index
    print "\t", cycle.atoms
    print "\t", cycle.bonds
    index = index + 1
