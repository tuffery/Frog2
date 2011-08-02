from frowns.perception import sssr, TrustedSmiles
from frowns import Smiles

s = 'c12[nH]c(c3ccc(CC)cc3)[nH0]c2[nH0]ccc1'
m = Smiles.smilin(s, transforms=[sssr.sssr, TrustedSmiles.trusted_smiles])

print m.arbsmarts()
for atom in m.atoms:
    print atom.symbol, atom.hcount, atom.explicit_hcount
