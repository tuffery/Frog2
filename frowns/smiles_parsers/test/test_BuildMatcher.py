import time
import Smarts, BuildMatcher
import test_Smarts

def test():
    lines = test_Smarts.lines
    h = BuildMatcher.BuildMatcher()
    t1 = time.clock()
    for line in lines:
        #print "-->", line
        Smarts.tokenize(line, h)
        s = h.mol.dump()
        #print s
    t2 = time.clock()

    print len(lines), "patt. in", t2 - t1, "s =>",
    print (t2-t1)/len(lines), "s/patt. ==",
    print len(lines)/(t2-t1), "patt./s"

if __name__ == "__main__":
    test()

