from frowns import Smiles
import time
from frowns.Depict import TkMoleculeDrawer

from Tkinter import *

class Generator:
    def __init__(self, smiles):
        self.smiles = smiles
        self.tk = Tk()
        button = Button(self.tk, text="advance", command=self.next)
        button.pack()
        self.drawer = None

    def next(self):
        if not self.smiles:
            return None
        smile = self.smiles.pop()
        m = Smiles.smilin(smile)

        #addCoords(m)
        #toplevel = Toplevel(self.tk)
        if self.drawer:
            self.drawer.frame.pack_forget()
        drawer = self.drawer=TkMoleculeDrawer.MoleculeDrawer(self.tk, m, name=smile)
        drawer.pack(expand=1, fill=BOTH)

file = open("../test/data/NCI_small")
lines = file.readlines()
smiles = [x.split()[0] for x in lines]
smiles = ["Oc1ccc(CC(C(=O)[O-])N)cc1", "c12ccccc1CC(=O)N2",     "CC12CCC3C(CCC4(O)CC(O)CCC34\\C=N\\NC(=O)C5=C(N=CN5)C(=O)NC6=NC=CC=C6)C1(O)CCC2C7=CC(=O)OC7"]
gen = Generator(smiles)
gen.next()
mainloop()
