from frowns import Smiles
from frowns.perception import figueras

badsmi = "C1=CC(=C2C=CC=C(N=CC34C5=C6C7=C3[Fe]456789%10%11C%12C8=C9C%10=C%11%12)C2=C1)N=CC%13%14C%15=C%16C%17=C%13[Fe]%14%15%16%17%18%19%20%21C%22C%18=C%19C%20=C%21%22"

m1 = Smiles.smilin(badsmi, transforms=[figueras.sssr])
m2 = Smiles.smilin(badsmi)

print len(m1.cycles)
print len(m2.cycles)

