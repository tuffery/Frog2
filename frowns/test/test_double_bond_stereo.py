UP = "UP"
DOWN = "DOWN"

def getDirectionRelativeToStereoCenter(bond, center):
    """(bond, center)->given a bond and the stereo center
    atom return the bond's direction relative to the stereo
    center."""
    
    a1, a2 = bond.atoms
    stereo = bond.stereo
    assert stereo
    assert center in bond.atoms
    
    if a1 is center:
        return stereo
    elif stereo == UP:
        return DOWN
    return UP


def getStructure(center, centerBond, stereoList):
    """(center, centerBond, stereoList)->return the stereo structure
    around the atom center.

    returned value is atom1, bond1, direction1, atom2, bond2, direction2
    where
    direction1 is the direction of  bond1 to atom1 relative to center and
    direction2 is the direction of bond2 to atom2 relative to center
    """
    a1, b1 = stereoList[0]
    s1 = getDirectionRelativeToStereoCenter(b1, center)
                
    # figure out the directions for the other bonds if necessary
    if len(stereoList) == 2:                
        a2, b2 = stereoList[1]
        s2 = getDirectionRelativeToStereoCenter(b2, center)
        assert s1 != s2, "Stereo atoms going in same direction"
    else:
        # get the other atom and get a direction for it
        if s1 == UP:
            s2 = DOWN
        else:
            s2 = UP
            
        for bond in center.bonds:
            if bond not in [centerBond, b1]:
                potentialAtom1, potentialAtom2 = bond.atoms
                if potentialAtom1 is center:
                    a2 = potentialAtom2
                else:
                    a2 = bontentialAtom1
                b2 = bond
                break
            else:
                a2, b2, s2 = None, None, s2

    return a1, b1, s1, a2, b2, s2

def assignStereo(m):
    """Assign stereo centers to the bond in m and figure out cis and trans
    atoms"""
    stereoAtoms = {}
    stereoCenters = {}

    # find the stereo bonds
    for b in m.bonds:
        if b.stereo:
            a1, a2 = b.atoms
            stereoAtoms[a1] = stereoAtoms.get(a1, []) + [(a2, b)]            
            stereoAtoms[a2] = stereoAtoms.get(a2, []) + [(a1, b)]

    # XXX FIX ME: this is nasty, nasty code
    for b in m.bonds:
        if b.bondorder == 2:
            atom1, atom2 = b.atoms
            if atom1 in stereoAtoms and atom2 in stereoAtoms:
                stereo1 = stereoAtoms[atom1]
                stereo2 = stereoAtoms[atom2]

                # we will now get four atoms and four bonds around the
                # stereo center
                #
                # like so:
                #  a1                 a3
                #    b1(s1)         b3(s3)
                #      atom1 = atom2
                #    b2(s2)         b4(s4)
                #  a2                 a4
                #
                # where s1, s2, s3 and s4 are the relative directions
                # from the stereo center atoms to the cis/trans atoms
                
                a1, b1, s1, a2, b2, s2 = getStructure(atom1, b, stereo1)
                a3, b3, s3, a4, b4, s4 = getStructure(atom2, b, stereo2)
                cis = []
                trans = []
                # now assign cis/trans
                if s1 == s3:
                    cis.append((a1, a3))
                    cis.append((a2, a4))
                    trans.append((a1, a4))
                    trans.append((a2, a3))
                else:
                    print "\t", a1, "trans", a3
                    print "\t", a1, "cis", a4
                    print "\t", a2, "cis", a3
                    print "\t", a2, "trans", a4
                 

                


def test():
    from frowns import Smiles
    strings = ["Cl/C(F)=C(Br)/I",
               #"Cl\C)(F)=C(Br)\I", <- fix parsing bug
               #"Cl\C(\F)=C(/Br)\I",
               "Cl\C(F)=C(Br)\I",
               "Cl\C(F)=C(Br)/I",
               "Cl/C(F)=C(Br)\I"
               ]
    for s in strings:
        print s
        m = Smiles.smilin(s)
        assignStereo(m)

if __name__ == "__main__":
    test()
