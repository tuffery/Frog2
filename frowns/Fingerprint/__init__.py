"""Codes for creating a daylight like fingerprint
for a molecule"""

import Fingerprint, Fingerlist, LinearPaths

# XXX FIX ME
# make the linear paths a generator eventually
# this should save memory...
def generateFingerprint(molecule, numInts=32, pathLength=7):
    """(molecule, numInts=32)->Fingerprint
    given a molecule and the number of integers to use to generate
    the fingerprint, return the fingerprint of the molecule)"""
    
    paths = LinearPaths.generatePaths(molecule, maxdepth=pathLength)
    fp = Fingerprint.Fingerprint(numIntegers=numInts)

    for path in paths:
        fp.addPath(path)

    return fp

def generateFingerlist(molecule, numInts=1024, pathLength=7):
    """(molecule, numInts=32)->Fingerprint
    given a molecule and the number of integers to use to generate
    the fingerprint, return the fingerprint of the molecule)"""
    
    paths = LinearPaths.generatePaths(molecule, maxdepth=pathLength)
    fp = Fingerlist.Fingerlist(numIntegers=numInts)

    for path in paths:
        fp.addPath(path)

    return fp

