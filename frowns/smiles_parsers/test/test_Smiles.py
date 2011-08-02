import string, time
import Smiles, Handler

patterns = """\
CC
[OH3+]
O=C=O
[2H]O[2H]
C#N
[235U]
CCN(CC)CC
F/C=C/F
CC(=O)O
F/C=C\F
C1CCCCC1
N[C@@H](C)C(=O)O
c1ccccc1
N[C@H](C)C(=O)O
C
P
N
S
O
Cl
[S]
[Au]
[Fe+++]
[Fe+3]
[H+]
[Fe+2]
[OH-]
[Fe++]
[OH3+]
[NH4+]
CC
C=O
C=C
O=C=O
COC
C#N
CCO
[H][H]
C=CCC=CCO
C=C-C-C=C-C-O
OCC=CCC=C
CCN(CC)CC
CC(C)C(=O)O
[O-]N(=O)C=CC(CCC)C(C(C)C)CCC
[O-]N(=O)c1ccccc1
c1ccccc1(N(=O)[O-])
C1CCCCC1
CC1=CC(Br)CCC1
CC1=CC(CCC1)Br
C12C3C4C1C5C4C3C25
O1CCCCC1N1CCCCC1
[Na+].[O-]c1ccccc1
c1cc([O-].[Na+])ccc1
C1.C1
[12C]
[13C]
[C]
[13CH4]
F/C=C/C=C/C
F/C=C/C=CC
NC(C)(F)C(=O)O
NC(F)(C)C(=O)O
N[C@@H](C(=O)O)C
N[C@@H](C)C(=O)O
N[C@@](F)(C)C(=O)O
N[C@@]([H])(C)C(=O)O
N[C@H](C(=O)O)C
N[C@H](C)C(=O)O
N[C@](C)(F)C(=O)O
N[C@]([H])(C)C(=O)O
[C@@H](N)(C)C(=O)O
[C@H](N)(C)C(=O)O
[H][C@@](N)(C)C(=O)O
[H][C@](N)(C)C(=O)O
OC(Cl)=[C@]=C(C)F
OC=[C@]=CF
OC(Cl)=[C@AL1]=C(C)F
OC([H])=[C@AL1]=C([H])F
F[Po@SP1](Cl)(Br)I
F[Po@SP2](Br)(Cl)I
F[Po@SP3](Cl)(I)Br
S[As@@](F)(Cl)(Br)C=O
O=C[As@](F)(Cl)(Br)S
S[Co@@](F)(Cl)(Br)(I)C=O
O=C[Co@](F)(Cl)(Br)(I)S
[2H]O[2H]
[H][H]
c1ccc1
c1ccccc1
n1ccccc1
O=n1ccccc1
[O-][N+]c1ccccc1
Cn1cccc1
[nH]1cccc1
C=[N+]=[N-]
C[N+](=O)[O-]
CN(=O)=O
C=[N]=[N]
O=c1[nH]cccc1
Oc1ncccc1
"""

lines = map(string.strip, string.split(patterns, "\n"))

def test():
    h = Handler.TokenHandler()
    t1 = time.clock()
    #h = Handler.WriteHandler()  # uncomment to print
    for line in lines:
        #print "####", line
        Smiles.tokenize(line, h)
    t2 = time.clock()

    print len(lines), "patt. in", t2 - t1, "s =>",
    print (t2-t1)/len(lines), "s/patt. ==",
    print len(lines)/(t2-t1), "patt./s"

if __name__ == "__main__":
    test()
