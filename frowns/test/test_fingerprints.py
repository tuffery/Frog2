import frowns.Fingerprint
from frowns import Smiles

patterns_and_targets = (
    ("CCN",
     ("CCN", "CCNCC", "c1cccc1CCN")),
    )

for pattern, targets in patterns_and_targets:
    pattern_molecule = Smiles.smilin(pattern)
    pfp = frowns.Fingerprint.generateFingerprint(pattern_molecule)

    for target in targets:
        mol = Smiles.smilin(target)
        molfp = \
              frowns.Fingerprint.generateFingerprint(mol)
        # pfp must be "in" molfp for test to pass
        assert pfp in molfp, "match %s is not in target %s!"%(pattern, target)

