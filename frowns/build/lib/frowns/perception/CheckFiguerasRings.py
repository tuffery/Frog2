# This code is hosed a bit... bummer
try:
    from kjbuckets import kjSet
except:
    import kjbuckets0 as kjbuckets
    
def checkRings(cycles):
    """This fixes a bug in figuras algorithm.  Sometimes when
    traversing through fused rings embedded in a larger ring
    two paths are returned when only one should have been
    This removes (most) of these cases"""
    # XXX I've encountered this bug before in figueras algorithm
    #     I really need to investigate if further and see
    #     if CDK has the same issues.
    res = []

    set = kjSet()
    sameSizes = {}
    i = 0
    for atoms, bonds in cycles:
        size = len(bonds)
        sameSizes[size] = sameSizes.get(size, []) + [(atoms, bonds)]

    keys = sameSizes.keys()
    keys.sort()
    last = kjSet()
    for size in keys:
        cycles = sameSizes[size]
        nextCycles = kjSet()

        i = 0
        if len(cycles) == 1:
            res.append(cycles[0])
            print "adding ring of size", cycles[0][1]
            continue

        dontAdd = {}
        for atoms1, bonds1 in cycles:
            s1 = kjSet(bonds1)
            j = i + 1
            for atom2, bonds2 in cycles[j:]:
                s2 = kjSet(bonds2)
                check = (s1-s2) + (s2-s1)
                l = len(check)

                if not (l and len(last & check) == l):
                    print "adding ring of size", len(bonds1)
                    res.append((atoms1, bonds1))
                else:
                    if i not in dontAdd:
                        dontAdd[j] = 1
                        res.append((atoms1, bonds1))
                        print "adding ring of size", len(bonds1)
                        print "dropping ring of size", len(bonds2)
        
            i += 1
        
                    
            nextCycles += s1
        last += nextCycles
    return res
