from Tkinter import *
from Tkinter import _cnfmerge
from ScrollingFrame import ScrollingFrame
from frowns.Depict.TkMoleculeDrawer import MoleculeDrawer

class MoleculeDock(ScrollingFrame):
    def __init__(self, master=None, cnf={}, **kw):
        """Construct a new Molecule dock with the parent master.

        Valid resource names: background, bd, bg, borderwidth, class,
        colormap, container, cursor, height, highlightbackground,
        highlightcolor, highlightthickness, relief, takefocus, visual, width."""

        # get the number of columns to use in the viewer
        if kw.has_key("cols"):
            self.cols = kw["cols"]
            del kw["cols"]
        else:
            self.cols = 2

        cnf = _cnfmerge((cnf, kw))
        ScrollingFrame.__init__(self, master)

        # now add the internal Scrolling Frame
        self._int = Frame(self)
        self._int.pack(expand=1, fill=BOTH)

        # current row and column
        self.r = 0
        self.c = -1
        for col in range(self.cols):
            self._int.grid_columnconfigure(col, weight=1)
        self._int.grid_rowconfigure(self.r, weight=1)
        
    def addMolecule(self, mol, name=None):
        if name is None:
            name = mol.name

        
        d = MoleculeDrawer(self._int,
                           mol,
                           name)
        self.c += 1
        if self.c >= self.cols:
            self.c = 0
            self.r += 1
            self._int.grid_rowconfigure(self.r, weight=1)

        d.grid(row=self.r, column=self.c, sticky="NEWS")


if __name__ == "__main__":
    from frowns import MDL
    reader = MDL.mdlin(open("../test/data/bad.sdf"))
    mol = reader.next()
    tk = top = Tk()

    m = MoleculeDock(top)
    m.pack(fill=BOTH, expand=1)
    m.addMolecule(mol)
    m.addMolecule(mol)
    m.addMolecule(mol)
    m.addMolecule(mol)

    mainloop()

        
