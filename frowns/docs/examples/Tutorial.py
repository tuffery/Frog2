from frowns import Smiles

mol = Smiles.smilin("c1ccccc1")

print mol.cansmiles()

print "atoms"
for atom in mol.atoms:
    print atom.symbol, atom.hcount, atom.aromatic

print "bonds"
for bond in mol.bonds:
    print bond.bondorder, bond.bondtype

