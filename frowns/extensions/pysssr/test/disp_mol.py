from Tkinter import *
from frowns.Depict.TkMoleculeDrawer import MoleculeDrawer
from frowns import MDL

# read in a molecule
reader = MDL.sdin(open("bad.sdf"))
mol, error, text = reader.next()

# create the moleculedock widget and place it
# into a tk window
tk = top = Tk()
m = MoleculeDrawer(top, mol, 0)
m.pack(fill=BOTH, expand=1)

mainloop()

        
