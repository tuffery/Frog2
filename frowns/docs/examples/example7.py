import frowns.Fingerprint
from frowns import Smiles

pattern = "CCN"
targets = ["CCN", "CCNCC", "c1cccc1CCN", "CC"]

pattern_molecule = Smiles.smilin(pattern)
pfp = frowns.Fingerprint.generateFingerprint(pattern_molecule)

for target in targets:
    mol = Smiles.smilin(target)
    molfp = \
          frowns.Fingerprint.generateFingerprint(mol)
    # pfp must be "in" molfp for test to pass
    if pfp in molfp:
        print "%s hits target %s"%(pattern, target)
    else:
        print "%s does not hit target %s"%(pattern, target)

