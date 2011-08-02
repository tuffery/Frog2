from frowns.perception import sssr, TrustedSmiles
from frowns import Smiles

for smiles in ["c1ccccc1", "c1cc1CCc2cc2"]:
    m = Smiles.smilin(smiles, transforms=[sssr.sssr, TrustedSmiles.trusted_smiles])

    print m.cansmiles()
    print m.arbsmarts()
    for bond in m.bonds:
        print bond.fixed, bond.symbol, bond.bondorder
    
fusedrings = "c12c3c(cc(c(Cl)c3oc2c(Cl)ccc1Cl)Cl)Cl"
for smiles in [fusedrings, 'c12[nH]c(c3ccc(CC)cc3)[nH0]c2[nH0]ccc1']:
    print "*"*44
    m = Smiles.smilin(smiles, transforms=[sssr.sssr, TrustedSmiles.trusted_smiles])

    print m.cansmiles()
    print m.arbsmarts()

    m = Smiles.smilin(smiles)
    print m.cansmiles()
    print m.arbsmarts()
