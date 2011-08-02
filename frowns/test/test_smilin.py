from test_smiles import smilesStrings
from frowns.Smiles import smilin
import time
#smilesStrings = ["Cc1[nH]ccc1"]
#smilesStrings = ["c1ccc[nH]1"]
##smilesStrings = ["c1ccccc1",
#smilesStrings = ["C=1=CC=CC=1"]
#smilesStrings = ["c1ccccc1"]

t1 = time.time()
for smile in smilesStrings:
    print "="*44
    print smile
    mol = smilin(smile)
    out = mol.arbsmiles()
    can = mol.cansmiles()
    for bond in mol.bonds:
        print bond.symbol, bond.bondorder, bond.bondtype, bond.fixed
    for atom in mol.atoms:
        print atom, atom.sumBondOrders()
        
    print smile, out, can

    
t2 = time.time()
print (t2-t1)/len(smilesStrings)
print len(smilesStrings)/(t2-t1)

