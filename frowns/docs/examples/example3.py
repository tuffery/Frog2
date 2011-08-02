from frowns import Smiles

print "*"*33
print "benzene with default transforms"
mol = Smiles.smilin("c1ccccc1")
print mol.cansmiles()

print "atoms"
for atom in mol.atoms:
    print "\tsymbol %s hcount %s aromatic %s"%(
        atom.symbol, atom.hcount, atom.aromatic)

print "bonds"
for bond in mol.bonds:
    print "\tbondorder %s, bondtype %s fixed %s"%(
        bond.bondorder, bond.bondtype, bond.fixed)

print "*"*33
print "benzene with no transforms"
mol = Smiles.smilin("c1ccccc1", transforms=[])
print mol.cansmiles()

print "atoms"
for atom in mol.atoms:
    print "\tsymbol %s hcount %s aromatic %s"%(
        atom.symbol, atom.hcount, atom.aromatic)

print "bonds"
for bond in mol.bonds:
    print "\tbondorder %s, bondtype %s fixed %s"%(
        bond.bondorder, bond.bondtype, bond.fixed)

print "*"*33
print "benzene with no transforms but bonds fully specified"
mol = Smiles.smilin("c1-c=c-c=c-c=1", transforms=[])
print mol.cansmiles()

print "atoms"
for atom in mol.atoms:
    print "\tsymbol %s hcount %s aromatic %s"%(
        atom.symbol, atom.hcount, atom.aromatic)

print "bonds"
for bond in mol.bonds:
    print "\tbondorder %s, bondtype %s fixed %s"%(
        bond.bondorder, bond.bondtype, bond.fixed)
