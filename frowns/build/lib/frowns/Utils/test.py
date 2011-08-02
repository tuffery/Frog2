# TIC library selection criteria
from frowns import utils, Smarts
# we want to remove steroids
steriod = Smarts.compile("C1CCC2C(C1)CCC3C2CC=C4CCCCC34")

# collect the molecules that we should keep
keep = file("MoleculesToKeep.sdf")
for molecule, error, record in MDLSDIN(file("PotentialMolecules.sf")):    
    if frowns.utils.saturationScore(molecule) < 0.5 and \
       frowns.utils.getStereoCenters(molecule) > 0 and \
       not steroid.match(molecule):
        # keep this molecule
        keep.write(record)
        
       
