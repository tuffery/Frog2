from Tkinter import *
from frowns.Depict.MoleculeDock import MoleculeDock
from frowns import Smiles

# read in a molecule
smiles = ["c1ccccc1C=O", "c1ccc2cc1cccc2",
          "CCN", "CCCC(=O)NNCC"]


# create the moleculedock widget and place it
# into a tk window
tk = top = Tk()
m = MoleculeDock(top)
m.pack(fill=BOTH, expand=1)

for smile in smiles:
    mol = Smiles.smilin(smile)
    m.addMolecule(mol)


mainloop()
