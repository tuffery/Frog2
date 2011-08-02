SP3 = "SP3"
SP2 = "SP2"
SP = "SP"
SD3 = "SD3"
DSP3 = "DSP3"
D2SP3 = "D2SP3"
ANY = None

# symbol, connected atoms, charge

_hybrid = {
    ("C", 4, 0): SP3,
    ("C", 3, -1): SP3,
    ("C", 3, 0): SP2,
    ("C", 4, +1): SP2,
    ("C", 2, ANY): SP,
    ("C", 1, 0): SP2,
    ("C", 1, -1): SP,
    
    ("N", 4, +1): SP3,
    ("N", 3, 0): SP3,
    ("N", 3, +1): SP2,
    ("N", 2, ANY): SP2,
    ("N", 1, -1): SP2,
    ("N", 1, 0): SP,

    ("O", 2, 0): SP3,
    ("O", 1, -1): SP3,
    ("O", 1, 0): SP2,
    ("O", 2, 1): SP2,

    ("P", 3, 0): SP3,
    ("P", 4, 1): SP3,
    ("P", 4, 0): SD3,
    ("P", 5, 0): DSP3,
    ("P", 6, 0): D2SP3,

    ("S", 2, 0): SP3,
    ("S", 1, -1): SP3,
    ("S", 1, 0): SP2,
    ("S", 3, 1): SP3,
    ("S", 3, 0): SD3,
    ("S", 4, 0): SD3,
}


def getHybridAtoms(molecule):
    """(molecule) -> mark all hybridiced atoms in a molecule
    that is set the atom.hybrid property to the correct value"""
    for atom in molecule.atoms:
        
        if atom.aromatic:
            # assume the combined bondorder is 3
            atom.hybrid = SP2
            continue
        symbol = atom.symbol
        hcount = atom.valences[0] - atom.sumBondOrders()
        connections = len(atom.bonds) + hcount
        charge = atom.charge
        if atom.symbol in ["C", "N"] and connections == 2:
            charge = ANY

        atom.hybrid =  _hybrid.get((symbol, connections, charge), "")
        #print atom, `atom`, hcount, atom.hybrid

        
        
            
if __name__ == "__main__":
    from frowns import Smiles
    from frowns.perception import TrustedSmiles, sssr
    
    smiles = "CC(=O)OC1C(N(c2c(CC1c1ccc(OC)cc1)cc(cc2)Oc1ccccc1)CCN(C)C)=O"
    smiles = "O=C1NC(N)=NC2=C1N=CNC2"
    #smiles = "O=C1C=CC(=O)C=C1"

    print smiles
    mol = Smiles.smilin(smiles, [sssr.sssr])
            
    getHybridAtoms(mol)
    for cycle in mol.cycles:
        print [(atom, atom.symbol, atom.hybrid) for atom in cycle.atoms]

