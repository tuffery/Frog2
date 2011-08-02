from frowns import Smiles
from frowns import Smarts

mol = Smiles.smilin("CCN")
pattern = Smarts.compile("CCN")

# simple match
match = pattern.match(mol)
assert match
index = 1
for path in match:
    print "match", index
    print "\tatoms", path.atoms
    print "\tbond", path.bonds
    index = index + 1

print "*"*33
# more complicated match
pattern = Smarts.compile("C*")
match = pattern.match(mol)
assert match
index = 1
for path in match:
    print "match", index
    print "\tatoms", path.atoms
    print "\tbond", path.bonds
    index = index + 1

print "*"*33
pattern = Smarts.compile("[!N]-[!C]")
match = pattern.match(mol)
assert match
index = 1
for path in match:
    print "match", index
    print "\tatoms", path.atoms
    print "\tbond", path.bonds
    index = index + 1
