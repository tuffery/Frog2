from frowns import Smiles
import random

# Take a molecule
#  read in the molecule, rearrange the atoms and bonds
#  and see if the output is still the same after
#  Canonicalization
#
#  This code won't work for the Editable Molecule class by the
#  way :)
smiles = "CC1=C(OP(=O)(OC2=C(C)C=CC=C2)OC3=C(C)C=CC=C3)C=CC=C1"

mol = Smiles.smilin(smiles)

smilefunc = mol.oldsmiles

template = smilefunc()

def _rearrange(l, numSwaps=10):
    l = l[:]
    length = len(l)
    swaps = []
    for i in range(numSwaps):
        a = random.randrange(0, length)
        b = random.randrange(0, length)
        l[a], l[b] = l[b], l[a]
        swaps.append((a,b))
    return l, swaps
                   
start_atoms = mol.atoms
start_equiv_class = [atom.equiv_class for atom in start_atoms]
start_sym_class = [atom.symclass for atom in start_atoms]
indices = [atom.index for atom in start_atoms]
assert indices == range(len(start_atoms))

for i in range(10000):
    atoms, swaps = _rearrange(start_atoms)
    mol.atoms = atoms
    mol.dirty = 1

    out = smilefunc()
    
    equiv_class = [atom.equiv_class for atom in start_atoms]
    sym_class = [atom.symclass for atom in start_atoms]
    if start_equiv_class != equiv_class:
        print start_equiv_class
        print equiv_class
    if sym_class != start_sym_class:
        print "oops!  Symclasses changed!"
        print start_sym_class
        print sym_class
        asdf

    if out != template:
        print "non canonical!!!"
        print swaps
        print "\t", out
        print "\t", template
        break
