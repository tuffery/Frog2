from bsddb3 import dbshelve
from frowns import Smiles, Smarts, Fingerprint

# ################################################################
# load database (this only needs to be done once
# you can comment this out at a later data
# ################################################################
# this is a persistent disk cache
cache = dbshelve.open("FOO.shelve")
molecules = ["CCCN", "c1ccccc1CCN", "C#C(N)CCc1ccccc1CCc1ccccc1"]

for m in molecules:
    mol = Smiles.smilin(m)
    fp = Fingerprint.generateFingerprint(mol)
    # cache the fingerprint using the canonical smiles
    # as the key and the fingerprint object as the value
    # shelves can store any picklable object!
    cache[mol.cansmiles()] = fp

cache.close()

# ##################################################################
# Iterate through database to see which molecules we need to
# do a full substructure search on
# ##################################################################
cache = dbshelve.open("FOO.shelve", 'r')

query = "CCCN" # this should reject molecule #2 above

testFP = Fingerprint.generateFingerprint(Smiles.smilin(query))

# this is how we iterate through a bsddb3 database
cursor = cache.cursor() # get a cursor
pair = cursor.first() # set it to the first key, value pair

while pair:           # pair will be None at end of iteration
    smiles, fp = pair # get the smiles and fingerprint
    if testFP in fp:
        print "need to substructure search", smiles    

    pair = cursor.next() # go to the next key, value pair

cache.close()
