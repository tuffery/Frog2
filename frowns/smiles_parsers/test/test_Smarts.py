import string, time
import Smarts, Handler

patterns = """\
C
c
a
[#6]
[Ca]
[++]
[R]
[D3]
[X3]
[v3]
C[C@H](F)O
C[C@?H](F)O
C
cc
c:c
c-c
[CH2]
[!C;R]
[!C;!R0]
[n;H]
[n&H]
[nH]
[c,n&H1]
[c,n;H1]
[Cl]
[35*]
[35Cl]
[F,Cl,Br,I]
*C
*CC
[$(*C);$(*CC)]
CaaO
CaaaN
Caa(O)aN
Ca(aO)aaN
C[$(aaO);$(aaaN)]
O
[O;H1]
[O;D1]
[O;D2]
[C,c]
[N;R]
[!C;R]
[n;H1]
[n&H1]
[c,n&H1]
[c,n;H1]
C.C
(C.C)
(C).(C)
(C).C
(C).(C).C
*!@*
*@;!:*
[C,c]=,#[C,c]
FC=CF
F/C=C/F
F/C=C/?F
F/C=C\F
F/C=C\?F
C[C@](F)(Cl)Br
C[C@?](F)(Cl)Br
C[Pt@SP](F)(Cl)Br
F[P@@](Cl)(Br)(I)c1ccccc1
F[P@TB2](Cl)(Br)(I)c1ccccc1
[#1]
[#6]
CCC(C)C
[X4]
[C;D1&H3,D2&H2,D3&H1,D4]
[CX4]
[Sv4]
[r5]
[r0]
[R0]
[R]
[!R0]
[H1]
[*H]
[H,+]
P(~[OD1])(~[OD1])~[OD1]
n
NC=O
NS(=O)=O
[Nr0]~[C,C+]([#6,#7])~[Nr0]
[OD1]~N~[OD1]
[CH3][N,O,S,P,F,Cl,Br,I]
[C@H](N)(C=O)C
[C@@H](N)(C=O)C
O=[C,N;R0]
[O,N;!H0]
O=[C,N;R0]-*~*-*=[O,N;!H0]
[OH,O-]C(=O)[$(c1o(Cl)ccoc1),$(c1cc(Cl)coc1),$(c1ccc(Cl)oc1)]
[!C]
[C&*]
[C&!*]
[C&C]
C[!*]C
C!~C
[a]
[!a]
[!R]
[C&O]
[a&!R]
[#6A]
[!#6,a]
[C&!O]
C1CC1
C1@C@C@1
[X1]*
[X1]-&!@*
[O&$([-])]
C#[$(N#C)]
C#N
[O-]
[$(C=O)$(CO)]
[$(C(=O)O)]
[$([OX1]C=O)$(OCC)]
[$([OX1]C(=O)C)]
[$(C=O)]
[C$(*=O)]
[O;-]
[O&-]
[C&r]
[Cr]
[a&s]
[as]
[+14]
[$([2H++])]
[$([2H]),$([3H])]
[2H,3H]
[$((C.C(C)C.C))]
[$((C))]
C(C)1CCC1
(C).(R).[U].([$(*)])
C[$([2H])$([3H])]
"""
lines = map(string.strip, string.split(patterns, "\n"))

def test():
    h = Handler.TokenHandler()
    #h = Handler.WriteHandler()  # uncomment to print
    t1 = time.clock()
    for line in lines:
        #print "Parsing", line
        Smarts.tokenize(line, h)
    t2 = time.clock()
    print len(lines), "patt. in", t2 - t1, "s =>",
    print (t2-t1)/len(lines), "s/patt. ==",
    print len(lines)/(t2-t1), "patt./s"

if __name__ == "__main__":
    test()
