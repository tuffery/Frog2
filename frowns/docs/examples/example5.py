from Tkinter import *
from frowns.Depict.MoleculeDock import MoleculeDock
from frowns import MDL

# read in a molecule
reader = MDL.mdlin(open("../../test/data/bad.sdf"))
mol = reader.next()

# create the moleculedock widget and place it
# into a tk window
tk = top = Tk()
m = MoleculeDock(top)
m.pack(fill=BOTH, expand=1)

# add some molecules
m.addMolecule(mol)
m.addMolecule(mol)
m.addMolecule(mol)
m.addMolecule(mol)

mainloop()

        
