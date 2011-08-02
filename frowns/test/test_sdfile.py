from frowns.MDL import sdin
from frowns.FrownsError import FrownsError

import sys

def getSymbolCounts(mol):
    counts = {}
    for atom in mol.atoms:
        counts[atom.symbol] = counts.get(atom.symbol, 0) + 1
    symbols = counts.keys()
    symbols.sort()
    l = []
    for s in symbols:
        num = counts[s]
        l.append("%s%s"%(s,num))
    return "".join(l)

def test(filename):
    print "starting!"

    file = open(filename)
    reader = sdin(file)

    error_file = open(filename + '.err', 'w')
    smiles_error_file = open(filename + '.smiles.err', 'w')

    empty_mols = open(filename + 'empty', 'w')

    dupes = {}

    num = 0
    numberBad = 0
    numEmpty = 0
    count = 0
    generateIDS = 0
    for mol, error, text in reader:
        num += 1
        if num % 100 == 0:
            print num, "done", numberBad, "bad", numEmpty, "empty"
            break

        if mol:
            smiles = mol.cansmiles()

            if not smiles:
                numEmpty += 1
                empty_mols.write(reader.reader.get_text())
            else:
                counts = getSymbolCounts(mol)
                smiles = (len(mol.atoms), len(mol.bonds), counts, smiles)
                id = mol.fields.get('ID', generateIDS)
                generateIDS += 1
                dupes[smiles] = dupes.get(smiles,[]) + [id]
                
        else:
            print error


    print 'number of compounds', len(dupes)
    numdupes = 0
    keys = dupes.keys()
    keys.sort()
    lastA, lastB, lastCounts = None, None, None
    for key in keys:
        ids = dupes[key]

        a, b, counts, smiles = key
        if (lastA, lastB, lastCounts) != (a, b, counts):
            print '*'*44
        lastA, lastB, lastCounts = (a, b, counts)
        
        if len(ids) > 1:
            print '*',
            numdupes += 1
        print key, ids


    print 'number of duplicates', numdupes


if __name__ == "__main__":
    test('data/test_6.sdf')
    

        


    
    
