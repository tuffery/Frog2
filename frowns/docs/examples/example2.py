from frowns import Smiles

listOfSmiles = ["CCN", "NCC", "CCC"]

duplicates = {}
for smile in listOfSmiles:
      mol = Smiles.smilin(smile)
      canonicalString = mol.cansmiles()
      if duplicates.has_key(canonicalString):
         print "found duplicate molecule", smile
      else:
         duplicates[canonicalString] = 1

print len(duplicates), "unique molecules found"
