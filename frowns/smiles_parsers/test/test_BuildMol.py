import time
import Smiles, BuildMol
import test_Smiles

def test():
    lines = test_Smiles.lines
    h = BuildMol.BuildMol()
    t1 = time.clock()
    for line in lines:
        #print line
        Smiles.tokenize(line, h)
        #s = h.mol.dump()
        #print s
    t2 = time.clock()

    print len(lines), "patt. in", t2 - t1, "s =>",
    print (t2-t1)/len(lines), "s/patt. ==",
    print len(lines)/(t2-t1), "patt./s"

if __name__ == "__main__":
    test()

