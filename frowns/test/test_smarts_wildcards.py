from frowns import Smiles
from frowns import Smarts
from frowns.perception import sssr
from frowns.perception import RingDetection

## '*' is ok
mol=Smiles.smilin("SCCCCSCCCC",transforms=[RingDetection.sssr])
#pattern=Smarts.compile("S*CC")

##!! '~' has problem
pattern=Smarts.compile("S~C~C")

##!! pattern can't match itself
mol=Smiles.smilin("C[CH](C1=CC=CC=C1)[C](C)(C#N)C2=CC=CC=C2")
pattern = Smarts.compile(mol.arbsmarts())
print mol.arbsmarts()
pattern=Smarts.compile("CC(C1=CC=CC=C1)C(C)(C#N)C2=CC=CC=C2")

match=pattern.match(mol)
assert match
index=1
for path in match:
    print "match",index
    print "\tatoms",path.atoms
    print "\tbond",path.bonds
    index=index+1
