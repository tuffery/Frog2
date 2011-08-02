from MoleculeDrawer import DrawMolHarness
from Tkinter import *

class TkMixin:    
    def _init(self, master, height, width):
        frame = self.frame = Frame(master, height=height, width=width)
        canvas = self.canvas = Canvas(frame, bg="white", height=height,
                                      width=width)
        self.width = width or int(master.cget("width"))
        self.height = height or int(master.cget("height"))

        canvas.pack(fill=BOTH, expand=1)
        self.canvas.bind("<Configure>", self._resize)

    def _drawOval(self, x, y, xh, yh):
        self.canvas.create_oval(x, y, xh, yh)
        
    def _drawText(self, text, font, fontsize, x, y, color, bg="white"):
            current = self.canvas.create_text(x,
                                              y,
                                              font=("Helvetica", fontsize),
                                              text=text,
                                              fill=color)
            bbox = self.canvas.bbox(current)
            i = self.canvas.create_rectangle(bbox,
                                             fill=bg,
                                             outline=bg,
                                             tag="marker")
            self.canvas.lower(i, current)

    def _clear(self):
        # returns -1 on failure
        for item in self.canvas.find_all():
            self.canvas.delete(item)

        if not self.molecule:
            self.canvas.create_text(
                self.width/2,
                self.height/2,
                text="UNABLE TO DEPICT", fill="red")
            return -1

        
        self.canvas.create_text(self.width/2,
                                self.height-self.border, 
                                text=self.name, fill="red")

        return 0
        
    def _resize(self, event):
        """(event) -> resive the drawing to event.height, event.width"""
        # the label might dissappear but so what for now...
        self.height = event.height
        self.width = event.width
        self.draw()
        
    def _drawLine(self, x1, y1, x2, y2, color):
        self.canvas.create_line(x1, y1, x2, y2, fill=color)
        
    def pack(self, *a, **kw):
        apply(self.frame.pack, a, kw)#

    def pack_forget(self):
        self.frame.pack_forget()
        
    def grid(self, *a, **kw):
        apply(self.frame.grid, a, kw)

    def postscript(self, *a, **kw):
        """return a postscript image of the current molecule arguments
        are sent to the Tkinter canvas postscript method"""
        return apply(self.canvas.postscript, a, kw)



class TkDrawer(DrawMolHarness, TkMixin): pass

MoleculeDrawer = TkDrawer

if __name__ == "__main__":
    from frowns import MDL, Smiles
    from ScrollingFrame import ScrollingFrame
    reader = MDL.mdlin(open("../test/data/bad.sdf"))
    #reader = MDL.mdlin(open("../apps/SaltFinder/SAVED_MOLECULE0.sdf"))

    mol = reader.next()
    tk = Tk()
    frame = Frame(tk)
    
    frame.pack(expand=1, fill=BOTH)

    drawer = MoleculeDrawer(frame, mol, drawAromatic=1)
    drawer.grid(row=0, column=0, sticky=N+E+W+S)
    frame.rowconfigure(0, weight=1)
    frame.columnconfigure(0, weight=1)
    
    
##    mol = reader.next()
##    components = getComponents(mol)
##    for component in range(len(components)):
##        f = Frame(frame, border=2, relief="sunken")
##        f.pack()
##        b = Button(f, text="keep", command=lambda x=None:x)
##        b.pack(side=LEFT)
##        d = DrawMolHarness(f, mol, name=mol.name,
##                           component=component)
##        d.pack(expand=1, fill=BOTH, side=LEFT)

    mainloop()
        
