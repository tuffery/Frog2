from frowns.Atom import Atom
from frowns.Bond import Bond

def addHydrogens(mol):
    if not mol.explicitHydrogens:
        for atom in mol.atoms:
            for i in range(atom.hcount):
                hatom = Atom()
                hatom.symbol = 'H'
                
                #is the coordinates of hydrogen stored somewhere in the
                #atom the hydrogen is attached to?
##                hatom.x = 
##                hatom.y = 
##                hatom.z =

                #does the hcount need to be changed?
##                atom.hcount -= 1 
                
                mol.add_atom(hatom)
                bond = Bond()
                mol.add_bond(bond, atom, hatom)

        #reset atom indices
        index = 0
        for atom in mol.atoms:
            atom.index = index
            index += 1

        mol.explicitHydrogens = 1 #change flag to let know that hydrogens are explicit in this mol 
    return mol


from frowns import Smiles
m = Smiles.smilin("CCCCc1ccccc1")
m.explicitHydrogens = 0

addHydrogens(m)

print m.cansmiles()
