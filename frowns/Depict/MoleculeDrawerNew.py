from Drawables import generateDrawables
from LineGenerator import make_double_line
from Tkinter import *
import GraphViz

GeneratorError = GraphViz.Error


####################################################
# Routines to split the components of a molecule
def _recurse(atom, visited):
    for bond in atom.bonds:
        for next in bond.atoms:
            if not visited.has_key(next):
                visited[next] = 1
                _recurse(next, visited)
    
                
def getComponents(molecule):
    atoms = {}
    for atom in molecule.atoms:
        atoms[atom] = 1

    components = []

    while atoms:
        start = atoms.keys()[0]
        visited = {start:1}
        _recurse(start, visited)
        bonds = {}
        for atom in visited.keys():
            del atoms[atom]
            for bond in atom.bonds:
                bonds[bond] = 1
                
        components.append((visited.keys(), bonds.keys()))

    return components

#        if not molecule:
#            self.canvas.create_text(
#                self.width/2,
#                self.height/2,
#                text="UNABLE TO DEPICT", fill="red")
#            return -1
#
#        
#        self.canvas.create_text(self.width/2,
#                                self.height-self.border, 
#                                text=self.name, fill="red")

class Drawer:
    def __init__(self, drawAromatic=1, coordinateGenerator=GraphViz.addCoords):
        self.drawAromatic = drawAromatic
        self.points, self.lines, self.bbox = None, None, None
        self.border = 24
        self.coordinateGenerator = coordinateGenerator
        self.md = None
        self.components = None

    def setMolecule(self, molecule, name=None, atoms=None, bonds=None):
        """(molecule, name, atoms=None, bonds=None) set the molecule
        to be drawn.  Use only atoms and bonds if a portion is desired and drawAromatic
        controls whether aromatic rings are drawn"""
        self.molecule = molecule
        if name is not None:
            self.name = name
        else:
            self.name = molecule.name or ""

        self.components = None
        
        # try to generate coordinates if needed
        if self.molecule.atoms and self.molecule.atoms[0].x == None:
            try:
                self.coordinateGenerator(molecule)
            except GeneratorError:
                self.molecule = None

        if atoms or bonds:
            self.setupComponent(atoms, bonds)
        else:
            self.setupMolecule()

    def setComponentNumber(self, component):
        if self.components is None:
            self.components = getComponents(molecule)
        atoms, bonds = components[component]
        self.setupComponent(atoms, bonds)        


    def setupComponent(self, atoms, bonds):
        """(atoms, bonds) setup to draw a portion of the molecule using
        only atoms and bonds."""
        m = self.molecule
        if m is None:
            return
        points, lines, circles, self.bbox, self.md = generateDrawables(m.atoms, m.bonds, drawAromatic=0)
        self.points, self.lines, self.circles, bbox, self.md = \
                     generateDrawables(atoms, bonds, self.drawAromatic)

        
    def setupMolecule(self):
        """->setup to draw the entire molecule"""

        if self.molecule:
            m = self.molecule
            self.points, self.lines, self.circles, self.bbox, self.md =  \
                         generateDrawables(m.atoms, m.bonds, self.drawAromatic)
        else:
            self.points, self.lines, self.circles, self.bbox, self.md = None, None, None, None

    def draw(self, renderer):
        """draw the molecule to the canvas"""
        renderer.clear()

        border = self.border
        points, lines, bbox = self.points, self.lines, self.bbox 
        height = renderer.height-2*border
        width = renderer.width-2*border
 
        if renderer.height < border or renderer.width < border:
            return
        
        self.bbox = bbox
        if bbox == None or bbox == (None, None, None, None):
            return

        xdist = float(bbox[1] - bbox[0])
        ydist = float(bbox[3] - bbox[2])

        if xdist == 0.0:
            xdist = 0.1
            
        if ydist == 0.0:
            ydist = 0.1
            
        self.xscale = (width-2*border)/xdist
        self.yscale = (height-2*border)/ydist
        if self.xscale < self.yscale:
            self.scale = self.xscale
        else:
            self.scale = self.yscale
            
        self.xmin = bbox[0]
        self.ymin = bbox[2]
        addPoint = self.addPoint
        addLine = self.addLine
        addCircle = self.addCircle
        
        for bondtype, l in lines:
            addLine(renderer, bondtype, l)

        for circle in self.circles:
            x,y,radius=circle
            addCircle(renderer, circle)

        for p in points:
            addPoint(renderer, p)

    def addCircle(self, renderer, circle):
        x, y, radius = circle
        x = (x-self.xmin) * self.xscale + self.border
        y = (y-self.ymin) * self.yscale + self.border
        xradius = radius * self.xscale
        yradius = radius * self.yscale
        xl = x - xradius
        yl = y - yradius
        xh = x + xradius
        yh = y + yradius
        renderer.drawOval(xl, yl, xh, yh)
        
    def addPoint(self, renderer, point):
        border = self.border
        symbol, x, y = point
        x = (x-self.xmin) * self.xscale + border
        y = (y-self.ymin) * self.yscale + border

        try:
            fontsize = int((self.md * self.scale)/3.0)
        except OverflowError:
            fontsize = 0

        if fontsize > 24: fontsize=24
        
        if symbol and fontsize:
            renderer.drawText(symbol, "Helvetica", fontsize, x, y, "black")
            
    
    def addLine(self, renderer, bondtype, line):        
        border = self.border
        x1,y1 = line[0]
        x2,y2 = line[1]
        x1 = (x1-self.xmin) * self.xscale + border
        y1 = (y1-self.ymin) * self.yscale + border

        x2 = (x2-self.xmin) * self.xscale + border
        y2 = (y2-self.ymin) * self.yscale + border

        if bondtype in [2, 3]:
            a, b, c, d = make_double_line((x1, y1), (x2, y2), split=3.0)
            xa, ya = a
            xb, yb = b
            xc, yc = c
            xd, yd = d
            renderer.drawLine(xa, ya, xc, yc)
            renderer.drawLine(xb, yb, xd, yd)
            if bondtype == 3:
                renderer.drawLine(x1,y1,x2,y2)
        else:
            renderer.drawLine(x1,y1,x2,y2)



class TkMoleculeDrawer:    
    def __init__(self, master, height, width, drawAromatic=1):
        frame = self.frame = Frame(master, height=height, width=width)
        canvas = self.canvas = Canvas(frame, bg="white", height=height,
                                      width=width)
        self.width = width or int(master.cget("width"))
        self.height = height or int(master.cget("height"))

        canvas.pack(fill=BOTH, expand=1)
        self.drawer = Drawer(drawAromatic)
        self.canvas.bind("<Configure>", self._resize)

    def setMolecule(self, mol):
        self.drawer.setMolecule(mol)
        self._resize()
        
    def drawOval(self, x, y, xh, yh):
        self.canvas.create_oval(x, y, xh, yh)
        
    def drawText(self, text, font, fontsize, x, y, color, bg="white"):
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

    def clear(self):
        for item in self.canvas.find_all():
            self.canvas.delete(item)
        
    def _resize(self, event=None):
        """(event) -> resive the drawing to event.height, event.width"""
        # the label might dissappear but so what for now...
        if event:
            self.height = event.height
            self.width = event.width
        self.drawer.draw(self)
        
    def drawLine(self, x1, y1, x2, y2):
        self.canvas.create_line(x1, y1, x2, y2)
        
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




MoleculeDrawer = TkMoleculeDrawer

if __name__ == "__main__":
    from frowns import MDL
    from ScrollingFrame import ScrollingFrame
    reader = MDL.mdlin(open("../test/data/bad.sdf"))
    #reader = MDL.mdlin(open("../apps/SaltFinder/SAVED_MOLECULE0.sdf"))
    mol = reader.next()
    tk = Tk()
    frame = Frame(tk)
    
    frame.pack(expand=1, fill=BOTH)

    drawer = TkMoleculeDrawer(frame, 256, 256, drawAromatic=1)
    drawer.setMolecule(mol)
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
        
