from frowns import Smiles, Smarts
for smiles in ['OC1C=CC(=O)C=C1', 'c1ccccc1']:
    mol = Smiles.smilin(smiles)
    print mol.cansmiles()

    print "with square brackets"    
    pattern = Smarts.compile('OC(C)[A]')
    match = pattern.match(mol)
    if match:
        for path in match:
            print path.atoms

    print
    print "without square brackets" 
    pattern2 = Smarts.compile('O=C(C)A')
    match2 = pattern2.match(mol)
    if match2:
        for path2 in match2:
            print path2.atoms

m = Smiles.smilin("c1ccccc1")
p = Smarts.compile("c1ccccc1")
assert p.match(m)


