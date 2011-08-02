from frowns import Smarts, Smiles

mol = Smiles.smilin('O=[N+]([O-])c1ccc(NC(C(C)(C)O)=O)cc1')
pattern = Smarts.compile('[cH0]-[!C;!N]')
match = pattern.match(mol)
if match:
    for path in match:
        a1, a2 = path.atoms
        b1 = path.bonds[0]
        print a1, a1.index, a1.aromatic, b1.aromatic, b1.bondtype, b1.bondorder, a2, a2.index, a2.aromatic


