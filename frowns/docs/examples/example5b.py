from Tkinter import *
from frowns.Depict.MoleculeDock import MoleculeDock
from frowns import Smiles

# read in a molecule
smiles = ["c1ccccc1C=O",
          "c1ccc2cc1cccc2",
          "CCN",
          "CCCC(=O)NNCC",
          "CCC(CC)([CH](OC(N)=O)C1=CC=CC=C1)C2=CC=CC=C2",
          "ON=C(CC1=CC=CC=C1)[CH](C#N)C2=CC=CC=C2",
          "C1=CC=C(C=C1)C(N=C(C2=CC=CC=C2)C3=CC=CC=C3)C4=CC=CC=C4",
          ]


# create the moleculedock widget and place it
# into a tk window
tk = top = Tk()
m = MoleculeDock(top, cols=4)
m.pack(fill=BOTH, expand=1)

for smile in smiles:
    mol = Smiles.smilin(smile)
    mol.name = smile
    m.addMolecule(mol)

mainloop()
