closures = {}
for i in range(10):
    closures[i] = str(i+1)

for i in range(9, 99):
    closures[i] = "%" + str(i+1)

def tweak_line(line, N):
    line[0] = line[0] + closures[N]
    line[-1] = line[-1] + closures[N]

def make_lines(N):
    first = []
    middle = []
    last = []
    for i in range(N):
        middle.append("C" + closures[i])
        if i % 2 == 0:
            first.append("C" + closures[i])
            last.append("C")
        else:
            first.append("C")
            last.append("C" + closures[i])

    tweak_line(first, N)
    tweak_line(middle, N)
    tweak_line(last, N)
    return "".join(first), "".join(middle), "".join(last)

def nanotube(width, height):
    N = width*2
    first, middle, last = make_lines(width*2)
    lines = [first]
    for i in range(height-1):
        lines.append(middle)
    lines.append(last)
    return ".".join(lines)

tube = nanotube(4, 100)
print tube

from frowns import Smiles
import time
#try:
#    import psyco
#    psyco.full()
#except:
#    print "nope"
#    pass

Smiles.smilin("c1cccc1")
t1 = time.time()
mol = Smiles.smilin(tube)
t2 = time.time()

print mol.cansmiles()

print t2-t1
