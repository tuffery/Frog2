from frowns import MDL
from frowns.Depict.MoleculeDrawer import *
file = open("data\\bad.sdf")

reader = MDL.mdlin(file, stripHydrogens=0)
mol = reader.next()


from Tkinter import *
root = Tk()
d = DrawMolHarness(root, mol, "bad")
d.pack(expand=1, fill=BOTH)

mainloop()
    
