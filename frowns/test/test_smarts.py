# Just some simple postive testing for smarts
# negative testing is a little harder to come by

from frowns import Smiles
from frowns import Smarts

mol = Smiles.smilin("C1CCCCC1CCN(=O)OCCN=C")

# search for ring atoms
matcher = Smarts.compile("[R1]")
print matcher.dump()
pathset = matcher.match(mol)
assert pathset
print pathset.atoms
print pathset.bonds

# search for CCN=O or CCN=C
matcher = Smarts.compile("CCN=[O,C]")
print matcher.dump()
pathset = matcher.match(mol)
assert pathset
for path in pathset:
    print path.atoms, path.bonds


# search for a wildcard
matcher = Smarts.compile("*=O")
print matcher.dump()
pathset = matcher.match(mol)
assert pathset
for path in pathset:
    print path.atoms, path.bonds

# from Identification of Biologically Active Profiles
# J. Cham. Ing. Comput. Sci., Vol 38, No. 2, 1998

HBD = "[!#6;!H0]"

matcher = Smarts.compile(HBD)
print matcher.dump()
pathset = matcher.match(mol)
assert pathset
for path in pathset:
    print path.atoms, path.bonds

# this will fail with  RunTimeError because
# the recursive matching is not yet implemented
HBA = "[$([!#6;+0]);!$([F,Cl,Br,I]);!$([o,s,nX3]);!$([Nv5,Pv5,Sv4,Sv6])]"
matcher = Smarts.compile(HBA)
print matcher.dump()
pathset = matcher.match(mol)
assert pathset
for path in pathset:
    print path.atoms, path.bonds
