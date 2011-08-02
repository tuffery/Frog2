"""Fingerprint
 generate fingerprints from molecules.

 fp = Fingerprint(numIntegers)
   Generate a fingerprint object that stores fingerprints
   in an array of numIntegers.

 fp.addPath(path)
   add a path to the fingerprint.  path is any str'able
   value.

 fp in fp2
   returns 1 if the fingerprint fp is a contained in the fingerprint
   fp2.

   This might be a bastardazation of __contains__ but I sort
   of like it.
   
 The method for creating fingerprints is remarkably simple.

 a sequence of non branching paths is extracted from a molecule.
 A path is a string value named atoms and named bonds
 through the traversal.  For instance:
 'C-C-C-C-N' or even 'Carbon through single bond to Carbon'.
 Any string will do as long as the string is always the
 same for the same path through the same or other molecules.

 For each path
  1 convert the string to an integer value and use it to
    seed a random number generator
    (random.seed can use any hashable value as the seed!
     python, again, is cool!)
    random.seed(path)
  2 pull out two random integers from the seeded generator
    index = int(random.random() * NUM_INTS)
    bit   = int(random.random() * INT_SIZE)

  fingerprint[index] = fingerprint[index] | 1<<bit

 we store a fingerprint as an array of integers.  Each integer
 has a certain number of bits that can be flipped.  The process
 of adding a path to a fingerprint is simply choosing the index
 and bit position for a path.  The above procedure does this
 in a deterministic fashion.
"""
import random
import sys
if sys.version_info[0:2] == (2,3):
    random = random.WichmannHill()
del sys

# We also need a good way to cache these to disk.  I suppose
# cPickle will work just fine for now...
class Fingerprint:
    def __init__(self, numIntegers=1, fingerprint=None):
        if fingerprint:
            self.fingerprint = fingerprint
        else:
            self.fingerprint = [0] * numIntegers

    def addPath(self, path):
        random.seed(path)
        fingerprint = self.fingerprint
        index = int(random.random() * len(fingerprint))
        # XXX FIX ME
        # might not be 32 on all platforms!
        bit = int(random.random() * 32)
        fingerprint[index] |= 1 << bit
        assert fingerprint[index] != 0
        return index, bit

    def __contains__(self, other):
        # self contains other
        if len(other.fingerprint) != len(self.fingerprint):
            raise "Fingerprints not the same size!"
        
        for a,b in zip(self.fingerprint, other.fingerprint):
            if(a & b) != b:
                return 0
        return 1

    def intTobitlist(self, i):
        # we'll assume that i is indeed a real integer
        result = []
        for loop in range(32):
            bit = i&1
            result.append(bit)
            i = i >> 1
        result.reverse()
        return result
    
    def to_list(self):
        result = []
        int_to_bitlist = self.intTobitlist
        
        for i in self.fingerprint:
            result += int_to_bitlist(i)
        return result

class SplitFingerprint:
    def __init__(self, maxdepth=7, integersPerPrint=[4]*6):
        assert maxdepth-1 == len(integersPerPrint)
        p = self.paths = []
        for numints in integersPerPrint:
            p.append(Fingerprint(numIntegers=numints))

    def addPath(self, length, path):
        self.paths[length-2].addPath(path)

    def __contains__(self, other):
        # does self contain other
        assert len(self.paths) == len(other.paths)
        for pother, pself in zip(other.paths, self.paths):
            if pother not in pself:
                return 0

        return 1

    def to_list(self):
        result = []
        for p in self.paths:
            result += p.to_list()
        return result
    
def test():
    # XXX FIX ME
    # better tests are needed of course :)
    paths = ["my",
             "dog",
             "has",
             "fleas"]
    
    fp = Fingerprint(32)
    fp2 = Fingerprint(32)
    for path in paths:
        fp.addPath(path)

    paths.reverse()
    for path in paths:
        fp2.addPath(path)
        
    assert fp.fingerprint == fp2.fingerprint
    assert fp in fp2

    fp2.addPath("yuma")
    assert fp in fp2
    assert fp2 not in fp

if __name__ == "__main__":
    test()
