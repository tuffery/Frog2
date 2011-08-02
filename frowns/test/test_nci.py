from frowns.Smiles import smilin
from frowns.perception import RingDetection, BasicAromaticity, figueras, sssr, rings

transform1 = [
    rings.sssr,
    BasicAromaticity.aromatize
    ]

transform2 = [RingDetection.sssr,
              BasicAromaticity.aromatize]

lastSmiles = None
def test():
    global lastSmiles
    # test NCI data
    file = open("data/NCI_aug00_SMI")
    import time

    t1 = time.time()
    index = 0
    for line in file.xreadlines():
        smiles, id = line.split()
        lastSmiles = smiles
        mol1 = smilin(smiles, transform1)
        mol2 = smilin(smiles, transform2)

        if mol1.cansmiles() != mol2.cansmiles():
            "%s different smiles->%s %s"%(smiles,
                                             mol1.cansmiles(), mol2.cansmiles())
        
        # check to see if the rings are the same
        if len(mol1.cycles) != len(mol2.cycles):
            print "%s number of rings are different"%smiles

        rings1 = {}
        
        keys = {}
        cycles1 = mol1.cycles
        for cycle in mol1.cycles:
            rings1[len(cycle)] = rings1.get(len(cycle), 0) +1
            atoms, bonds = cycle.atoms, cycle.bonds
            atoms.sort()
            bonds.sort()
            keys[(tuple(atoms), tuple(bonds))] = 1

        mol1 = RingDetection.sssr(mol1)
        rings2 = {}
        for cycle in mol1.cycles:
            rings2[len(cycle)] = rings1.get(len(cycle), 0) +1
            atoms, bonds = cycle.atoms, cycle.bonds
            atoms.sort()
            bonds.sort()
            key = tuple(atoms), tuple(bonds)
            if not keys.has_key(key):
                print smiles, "->rings are different"
                break        

        lengths1 = [len(cycle) for cycle in cycles1]
        lengths2 = [len(cycle) for cycle in mol1.cycles]
        lengths1.sort()
        lengths2.sort()
        if lengths1 != lengths2:
            print "%s forms rings of different sizes"%smiles
            print "\t",lengths1
            print "\t",lengths2

        index += 1
        if max and index >= max:
            break


    t2 = time.time()
    print "tot time", t2-t1
    print "compound time", (t2-t1)/index
    print "compounds/sec", index/(t2-t1)

test()

############################################3
## uncomment to profile
##import profile
##profile.run("test()", 'parser_profiler')

##import pstats
##p = pstats.Stats('parser_profiler')
##p.sort_stats('cum')
##p.print_stats()

