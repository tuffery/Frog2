"""TrustedSmiles

exports
  molecule = trusted_smiles(molecule)
    returns a molecule that trusts the original input for the
    declaration of aromaticity.


 This is a really stupid algorithm to accept standard smiles input
 and basically trust it to be correct.

 Essentially all non specified bonds are converted to either
 aromatic or single bonds depending on two things

 An aromatic bond must:
  1 be between two aromatic atoms
  2 must be in a ring with either unspecified
    bonds or aromatic bonds.

 Otherwise it is considered a single bond.

 BUG
 There is one class of degenerate cases for this stupid
 algorithm for pre Daylight 4.71 output


     ___
    / o \       Imagine three aromatic rings surrounding
    \___/        a non aromatic ring
 ___/   \___    If the non aromatic ring doesn't have it's
/ o \___/ o \    single bonds specified, then the single bonds
\___/   \___/    and hence the inner ring will be considered aromatic.

 This algorithm really is for testing purposes, so don't get your
 knickers in a bunch :)
"""
def atoms_are_aromatic(atoms):
    for atom in atoms:
        if not atom.aromatic:
            return 0
    return 1

def bonds_are_unspecified_or_aromatic(bonds):
    for bond in bonds:
        if bond.fixed and not bond.aromatic:
            return 0
    return 1

def trusted_smiles(molecule):
    """(molecule)->molecule
    Trust the original smiles string in order to generate the proper
    aromaticity.  Will produce odd results with some smiles strings
    produced by early versions of Daylight (pre 4.71)"""



    for cycle in molecule.cycles:
        if atoms_are_aromatic(cycle.atoms) and \
           bonds_are_unspecified_or_aromatic(cycle.bonds):
            cycle.set_aromatic()

    # set the rest to single bonds
    for bond in molecule.bonds:
        if not bond.fixed:
            bond.bondtype = bond.bondorder = 1
            bond.symbol = "-"
            bond.fixed = 1
            
    # check the valences of the atoms
    # set the symmetry order to the input order
    index = 0
    for atom in molecule.atoms:
        atom.symorder = index
        if atom.valences:
            lowestValence = atom.valences[0]
            sumBondOrders = atom.sumBondOrders()
            original_hcount = atom.hcount + atom.explicit_hcount
            atom.hcount = max(original_hcount, lowestValence + atom.charge - sumBondOrders)
            atom.imp_hcount = atom.hcount - original_hcount

##            if atom.hcount + sumBondOrders > lowestValence:
##                if atom.aromatic and atom.symbol == "N" and atom.hcount == 1:
##                    # pyrole nitrogen (maybe)
##                    # skip it
##                    pass
##                #else:
##                #    print "Warning High Valence Atom", atom.index
        index += 1
    

    return molecule

    
