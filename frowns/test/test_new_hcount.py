from frowns import Smiles
m = Smiles.smilin("[CH0]")
assert m.cansmiles() == "[C]"
m = Smiles.smilin("C")
assert m.cansmiles() == "C"
