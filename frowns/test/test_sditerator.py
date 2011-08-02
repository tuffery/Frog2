from frowns.mdl_parsers import sditerator
import time

def test(filename):
    t1 = time.time()
    c = 0
    for m, text, err in sditerator.reader(open(filename)):
        if not m: print err
        c += 1
        if c == 100:
            break

    t2 = time.time()

    print t2-t1
    print (t2-t1)/c
    
if __name__ == "__main__":
    import profile
    profile.run("test('f:\\Chemical Libraries\\TIC\\ALLMOLECULES.SDF')")

    test('f:\\Chemical Libraries\\TIC\\ALLMOLECULES.SDF')
