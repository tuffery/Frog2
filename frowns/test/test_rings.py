from test_smiles import smilesStrings
import frowns.perception.RingDetection

from frowns import Smiles

def test():
    hardSmiles = [
        "C12C3C4C1C5C4C3C25",
        "C12C3C4C1C5C4C3C25CCCCCC12C3C4C1C5C4C3C25"]
    
    for smiles in smilesStrings + hardSmiles:
        print smiles,

        mol = parse.parse(smiles)
        mol = frowns.perception.RingDetection.sssr(mol, 0)
        rings = {}
        for ring in mol.rings:
            rings[len(ring)] = rings.get(len(ring), 0) + 1

        keys = rings.keys()
        keys.sort()
        
        for k in keys:
            print "%s:%s"%(k, rings[k]),
        print

def test2():

    smiles, id = "OC1=C(Cl)C=CC=C1[N+]", 1
    mol = Smiles.smilin(smiles)
    print mol.cansmiles()
    for atom in mol.atoms:
        print atom.index

def test3():
    # test NCI data
    file = open("data/NCI_aug00_SMI")
    import time
    max = None
    t1 = time.time()
    index = 0
    for line in file.xreadlines():
        smiles, id = line.split()
        mol1 = smilin(smiles)
        mol2 = smilin(smiles, 0)

        print smiles, len(mol1.rings), len(mol2.rings)
        if len(mol1.rings) != len(mol2.rings):
            print "%s %s!=%s"%(smiles, len(mol1.rings), len(mol2.rings))
            raise "Ring error"

        index += 1
        if max is not None and index > max:
            break
    t2 = time.time()
    print "tot time", t2-t1
    print "compound time", (t2-t1)/index
    print "compounds/sec", index/(t2-t1)
    
test2()
