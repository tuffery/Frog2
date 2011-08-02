from Drawables import generateDrawables
from LineGenerator import make_double_line

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

class DrawMolHarness:
    def __init__(self,
                 master,
                 molecule=None,
                 name=None,
                 component=None,
                 height=256, width=256,
                 coordinateGenerator=GraphViz.addCoords,
                 drawAromatic=0):
        self.master = master
        self.molecule = molecule

        self._init(master, height, width)
        
        self.points, self.lines, self.bbox = None, None, None
        self.border = 24
        self.coordinateGenerator = coordinateGenerator
        self.md = None
        self.colors = {}
        if molecule:
            if type(component) == type(1):
                components = getComponents(molecule)
                atoms, bonds = components[component]
                self.setMolecule(molecule, name, atoms, bonds,
                                 drawAromatic=drawAromatic)
            elif component:
                atoms, bonds = component
                self.setMolecule(molecule, name, atoms, bonds,
                                 drawAromatic=drawAromatic)
            else:
                self.setMolecule(molecule, name, molecule.atoms, molecule.bonds,
                                 drawAromatic=drawAromatic)

    def setMolecule(self, molecule, name=None, atoms=None, bonds=None, drawAromatic=1):
        """(molecule, name, atoms=None, bonds=None, drawAromatic=1) set the molecule
        to be drawn.  Use only atoms and bonds if a portion is desired and drawAromatic
        controls whether aromatic rings are drawn"""
        self.molecule = molecule
        if name is not None:
            self.name = name
        else:
            self.name = molecule.name or ""
                
        if self.molecule.atoms and self.molecule.atoms[0].x == None:
            try:
                self.coordinateGenerator(molecule)
            except GeneratorError:
                self.molecule = None

        if atoms or bonds:
            self.drawComponent(atoms, bonds, drawAromatic=drawAromatic)
        else:
            self.drawMolecule(drawAromatic=drawAromatic)

    def draw(self):
        """draw the molecule to the canvas"""
        if self._clear() == -1:
            return
        
        border = self.border
        points, lines, bbox = self.points, self.lines, self.bbox
        if self.bbox == (None, None, None, None):
            return
        height = abs(self.bbox[2] - self.bbox[0])
        width = abs(self.bbox[3] - self.bbox[1])
        self._minDoubleLineSplit = min(height*0.02, width*0.02)
        height = self.height-2*border
        width = self.width-2*border
 
        if self.height < border or self.width < border:
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
        addWedge = self.addWedge
        
        for bondtype, bond, l in lines:
            addLine(bondtype, bond, l)

        for circle in self.circles:
            x,y,radius=circle
            addCircle(circle)

        for wedge in self.wedges:
            addWedge(wedge)
            
        for p in points:
            addPoint(p)

    def addCircle(self, circle):
        x, y, radius = circle
        x = (x-self.xmin) * self.xscale + self.border
        y = (y-self.ymin) * self.yscale + self.border
        xradius = radius * self.xscale
        yradius = radius * self.yscale
        xl = x - xradius
        yl = y - yradius
        xh = x + xradius
        yh = y + yradius
        self._drawOval(xl, yl, xh, yh)
        
    def addPoint(self, point):
        border = self.border
        symbol, x, y = point
        x = (x-self.xmin) * self.xscale + border
        y = (y-self.ymin) * self.yscale + border

        try:
            size = int((self.md * self.scale)/3.0)
            if size < 6:
                size = 6
            fontsize = size
##            fontsize = min(14, int((self.md * self.scale)/3.0))
        except OverflowError:
            fontsize = 0

        if symbol and fontsize:
            self._drawText(symbol, "Helvetica", fontsize, x, y, "black")
            
    def addWedge(self, wedge):
        bondtype, bond, c1, c2, c3, filled = wedge
        if not filled:
            self.addLine(1, bond, [c1, c2])
            self.addLine(1, bond, [c2, c3])
            self.addLine(1, bond, [c3, c1])
        else:
            xmin, ymin = self.xmin, self.ymin
            border = self.border
            c1 = (c1[0] - xmin) * self.xscale + border, \
                 (c1[1] - ymin) * self.yscale + border
            c2 = (c2[0] - xmin) * self.xscale + border, \
                 (c2[1] - ymin) * self.yscale + border
            c3 = (c3[0] - xmin) * self.xscale + border, \
                 (c3[1] - ymin) * self.yscale + border
            
            self._drawPoly([c1, c2, c3], color="black")
#            self.addFilledPoly(bond, [c1, c2, c3])
        
    def addLine(self, bondtype, bond, line):        
        border = self.border
        x1,y1 = line[0]
        x2,y2 = line[1]
        x1 = (x1-self.xmin) * self.xscale + border
        y1 = (y1-self.ymin) * self.yscale + border

        x2 = (x2-self.xmin) * self.xscale + border
        y2 = (y2-self.ymin) * self.yscale + border
        drawLine = self._drawLine
        color = self.colors.get(bond, "black")
        
        if bondtype in [2, 3]:
            a, b, c, d = make_double_line((x1, y1), (x2, y2), split=self._minDoubleLineSplit)
            xa, ya = a
            xb, yb = b
            xc, yc = c
            xd, yd = d
            drawLine(xa, ya, xc, yc, color)
            drawLine(xb, yb, xd, yd, color)
            if bondtype == 3:
                drawLine(x1, y1, x2, y2, color)
        else:
            drawLine(x1, y1, x2, y2, color)


    def drawComponent(self, atoms, bonds, drawAromatic=1):
        """(atoms, bonds, drawAromatic=1) draw a portion of the molecule using
        only atoms and bonds.  set drawAromatic to zero if aromatic rings are
        not desired"""
        m = self.molecule
        if m is None:
            return
        points, lines, circles, self.wedges, self.bbox, self.md = generateDrawables(m.atoms, m.bonds, drawAromatic=0)
        self.points, self.lines, self.circles, self.wedges, bbox, self.md = \
                     generateDrawables(atoms, bonds, drawAromatic=drawAromatic)

        
    def drawMolecule(self, drawAromatic=1):
        """(drawAromatic=1) draw the current molecule set drawAromatic to 0
        if drawing of aromatic rings is not desired"""
        if self.molecule:
            m = self.molecule
            self.points, self.lines, self.circles, self.wedges, self.bbox, self.md =  \
                         generateDrawables(m.atoms, m.bonds, drawAromatic=drawAromatic)
        else:
            self.points, self.lines, self.circles, self.wedges, \
                         self.bbox, self.md = None, None, None, None


