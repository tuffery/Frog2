from frowns import Smiles

mol = Smiles.smilin("C(\C)(C)=C(/C)(C)")
print mol.cansmiles()

mol = Smiles.smilin("F/C=C/F")
print mol.cansmiles()

for bond in mol.bonds:
    print bond.symbol
    print bond.equiv_class
