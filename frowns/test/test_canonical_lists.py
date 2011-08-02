from frowns import Smiles

mol = Smiles.smilin("c1cccc1")

print len(mol.bonds)
print mol.cansmiles()
print map(len, mol.canonical_list[0])
