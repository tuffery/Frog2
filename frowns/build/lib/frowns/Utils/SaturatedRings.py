from frowns.perception import sssr

def getSaturationScore(molecule):
    if not hasattr(molecule, "cycles"):
        sssr.sssr(molecule)
        
    fusedRings = {}
    lookup = {}

    i = 0
    for cycle in molecule.cycles:
        fusedRings[i] = [cycle]
        lookup[cycle] = i
        i += 1

    if not fusedRings:
        return None
    
    for bond in molecule.bonds:
        for cycle in bond.rings:
            index = lookup[cycle]
            for cycle2 in bond.rings:
                if cycle is not cycle2:
                    rings = fusedRings[index]
                    if cycle2 not in rings:
                        rings.append(cycle2)

    # just keep one of the rings for inspection
    for a, cycles in fusedRings.items():
        for cycle in cycles:
            index = lookup[cycle]
            if index > a:
                if index in fusedRings:
                    del fusedRings[index]
    potential = 0
    actual = 0.0
    
    for a, cycles in fusedRings.items():
        bonds = {}

        shared = 0

        for cycle in cycles:
            for bond in cycle.bonds:
                if bond in bonds: shared += 1
                else: bonds[bond] = 1

        double_bonds = 0
        for bond in bonds:
            if bond.bondorder == 2:
                double_bonds += 1

        #print "number of shared bonds", shared
        #print "number of real bonds", len(bonds)
        addition = (len(bonds) - shared)/2
        #if double_bonds > addition:
            #print double_bonds, addition

            
        potential += addition
        actual += double_bonds


    if potential == 0:
        return None
    
    return actual/potential



if __name__ == "__main__":
    test = """C1CCCCC1.C1CCCCC1       0
CC1=CC(=O)C=CC1=O	1
S(SC1=NC2=CC=CC=C2S1)C3=NC4=C(S3)C=CC=C4	2
OC1=C(Cl)C=C(C=C1[N+]([O-])=O)[N+]([O-])=O	3
[O-][N+](=O)C1=CNC(=N)S1	4
NC1=CC2=C(C=C1)C(=O)C3=C(C=CC=C3)C2=O	5
OC(=O)C1=C(C=CC=C1)C2=C3C=CC(=O)C(=C3OC4=C2C=CC(=C4Br)O)Br	6
CN(C)C1=C(Cl)C(=O)C2=C(C=CC=C2)C1=O	7
CC1=C(C2=C(C=C1)C(=O)C3=CC=CC=C3C2=O)[N+]([O-])=O	8
CC(=NO)C(C)=NO	9
C1=CC=C(C=C1)P(C2=CC=CC=C2)C3=CC=CC=C3	10
CC(C)(C)C1=C(O)C=C(C(=C1)O)C(C)(C)C	11
CC1=NN(C(=O)C1)C2=CC=CC=C2	12
NC1=CC=NC2=C1C=CC(=C2)Cl	13
CCCCCC[CH]1CCCCN1	14
O=CC1=C2C=CC=CC2=CC3=C1C=CC=C3	15
BrN1C(=O)CCC1=O	16
CCCCCCCCCCCCCCCC1=C(N)C=CC(=C1)O	17
C(COC1=C(C=CC=C1)C2=CC=CC=C2)OC3=CC=CC=C3C4=CC=CC=C4	18
CCCCSCC	19
CC(=O)NC1=NC2=C(C=C1)C(=CC=N2)O	20
CC1=C2C=CC(=NC2=NC(=C1)O)N	21
CCOC(=O)C1=CN=C2N=C(N)C=CC2=C1O	22
CC1=CC(=NC=C1)N=CC2=CC=CC=C2	23
C[N+](C)(C)CC1=CC=CC=C1	24
C[N+](C)(C)C(=O)C1=CC=CC=C1	25
ICCC(C1=CC=CC=C1)(C2=CC=CC=C2)C3=CC=CC=C3	26
CC1=CC(=C(C[N+](C)(C)C)C(=C1)C)C	27
C[C](O)(CC(O)=O)C1=CC=C(C=C1)[N+]([O-])=O	28
CC1=CC=C(C=C1)C(=O)C2=CC=C(Cl)C=C2	29
ON=CC1=CC=C(O)C=C1	30
CC1=CC(=C(N)C(=C1)C)C	31
CC1=CC=C(C=C1)C(=O)C2=CC=C(C=C2)[N+]([O-])=O	32
CC(O)(C1=CC=CC=C1)C2=CC=CC=C2	33
ON=CC1=CC(=CC=C1)[N+]([O-])=O	34
OC1=C2C=CC(=CC2=NC=C1[N+]([O-])=O)Cl	35"""

    smiles_strings = [x.split()[0] for x in test.split("\n")]
    
    from frowns import Smiles
    for input in smiles_strings:
        m = Smiles.smilin(input)
        saturation_coef = getSaturationScore(m)
        if saturation_coef is None:
            label = "N/A"
        else:
            label = "%0.3f"%saturation_coef

        print "%s\t%s"%(label, input)

        
