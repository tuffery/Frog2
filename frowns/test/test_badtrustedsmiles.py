from frowns import Smiles
from frowns.perception import TrustedSmiles, sssr

smiles = "CC(=O)OC1C(N(c2c(CC1c1ccc(OC)cc1)cc(cc2)Oc1ccccc1)CCN(C)C)=O"
#smiles = "c1ccccc1-c2ccccc2"
smiles = "c12[nH0]c(N)[nH0]c(c2[nH0]c[nH0]1C1OC(C(C1)N=[N+]=[N-])CO)N"

mol = Smiles.smilin(smiles, transforms=[sssr.sssr,
                                        TrustedSmiles.trusted_smiles])

        
mol = Smiles.smilin(smiles)

for bond in mol.bonds:
    print bond, bond.bondtype, bond.bondorder
for cycle in mol.cycles:
    if cycle.aromatic:
        print "*"*44

        for bond in cycle.bonds:
            print bond, bond.bondtype, bond.bondorder

print "*"*44

