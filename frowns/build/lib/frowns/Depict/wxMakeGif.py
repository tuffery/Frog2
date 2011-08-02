"""Draw a picture form a molecule"""

from frowns.Depict.MoleculeDrawer import DrawMolHarness
from wxPython.wx import *
import os, tempfile

class Drawer(DrawMolHarness):
    def __init__(self, x0=0, y0=0, width=256, height=256):
        self.x0 = x0
        self.y0 = y0
        self.width = width
        self.height = height
        self.border = 5
        self._textFont = "ARIAL"
        self.fgBrush = wxBrush(wxWHITE, wxTRANSPARENT)
        self.bgBrush = wxBrush(wxWHITE, wxSOLID)

        self.fgPen = wxBLACK_PEN
        self.bgPen = wxWHITE_PEN
        
    def _clear(self):
        self.dc.DrawRectangle(self.x0, self.y0, self.width, self.height)
    
    def _drawOval(self, x, y, xh, yh):
        self.dc.DrawEllipse(x+self.x0, y+self.y0, xh-x, yh-y)
        
    def _drawLine(self, x1, y1, x2, y2, color):
        color = wxTheColourDatabase.FindColour(color)
        dc = self.dc
        pen = dc.GetPen()
        if color is not None:
            dc.SetPen(wxPen(color, 1, wxSOLID))
            
        dc.DrawLine(self.x0+x1, self.y0+y1, self.x0+x2, self.y0+y2)
        dc.SetPen(pen)

    def _drawPoly(self, coords, color="black"):
        x0,y0 = self.x0, self.y0
        coords = [(x+x0,y+y0) for x,y in coords]
        dc = self.dc
        pen = dc.GetPen()
        brush = dc.GetBrush()
        if color is not None:
            dc.SetPen(wxPen(color, 1, wxSOLID))
            dc.SetBrush(wxBrush(color, wxSOLID))
        self.dc.DrawPolygon(coords)
        dc.SetPen(pen)
        dc.SetBrush(brush)
        
    def _drawText(self, text, font, fontsize, x, y, color, bg="white"):
        if not text:
            # don't try to draw empty text strings.
            # this is just stupid...
            return
        dc = self.dc
        font = wxFont(fontsize, wxSWISS, wxNORMAL, wxNORMAL,
                      0, self._textFont)
        dc.SetFont(font)

        w,h = dc.GetTextExtent(text)
        w+=2
        h+=2
        x -= w/2
        y -= h/2
            

        pen = dc.GetPen()
        dc.SetPen(self.bgPen)
        dc.DrawRectangle(x+self.x0, y+self.y0, w, h)
        dc.SetPen(pen)
        dc.DrawText(text, x+self.x0, y+self.y0)

    def makeGif(self, mol, width=256, height=256, type=wxBITMAP_TYPE_BMP):
        self.height=height
        self.width=width
        self.x0=0
        self.y0=0

        class Rect:
            def __init__(self, **kw):
                self.__dict__.update(kw)

        rect = Rect(
            x = 0,
            y = 0,
            width = width,
            height = height
            )
        
        self.fgBrush = wxBrush(wxWHITE, wxTRANSPARENT)
        self.bgBrush = wxBrush(wxWHITE, wxSOLID)

        self.fgPen = wxBLACK_PEN
        self.bgPen = wxWHITE_PEN
        filename = tempfile.mktemp()
        dc = wxMemoryDC()
        bmp = wxEmptyBitmap(width, height)
        self._draw(mol, dc, rect)
        bmp.SaveFile(filename, type)
        text = open(filename).read()
        os.remove(filename)
        return text

    def _draw(self, mol, dc, rect):
        """actual drawing code for the widget.  This is seperated into a private
        function to allow for drawing on any supplied DC like a postscript or
        printer DC"""
        dc.SetBackgroundMode(wxSOLID)
        dc.SetBrush(self.bgBrush)
        dc.SetPen(self.fgPen)
        dc.DrawRectangle(rect.x, rect.y, rect.width, rect.height)
        dc.SetBrush(self.fgBrush)
        
        if mol.name:
            fontsize = min(16, (rect.width/float(2))/(2*len(mol.name)))
        else:
            fontsize = 0
        if not mol:
            return

        self.dc = dc
                    
        self.setMolecule(mol)
        self.colors = {}

#        # XXX FIX ME This probably shouldn't be here.
#        if hasattr(self.table, "bond_indices"):
#            bonds = [mol.bonds[i] for i in self.table.bond_indices[displayedRow]]
#            for bond in bonds:
#                self.colors[bond] = "red"

        self.border = 5

        self.x0 = rect.x
        self.y0 = rect.y
        self.width = rect.width
        self.height = rect.height - 2*fontsize
        self.draw()

        if mol.atoms:
            self._drawText(mol.name, "Arial", fontsize,
                           rect.width/2, rect.height-fontsize, "red")
        else:
            if mol.name:
                fontsize = min(16, float(rect.width)/len(mol.name))
                self._drawText(mol.name, "Arial", fontsize,
                               rect.width/2, rect.height/2, "red")
                

        # clear the dc when done
        self.dc = None

__test__ = """
  -ISIS-  05260006032D

 10 10  0  0  0  0  0  0  0  0999 V2000
    0.3541   -0.2083    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.3541   -1.0333    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3583    0.2083    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
   -0.3583    1.0333    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3583   -1.4416    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0708    1.4458    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0666   -0.2083    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0666   -1.0333    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3583    1.4458    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0708   -1.4458    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  1  3  1  0  0  0  0
  2  5  1  0  0  0  0
  2 10  1  0  0  0  0
  3  4  1  1  0  0  0
  3  7  1  0  0  0  0
  4  6  2  0  0  0  0
  4  9  1  0  0  0  0
  5  8  1  0  0  0  0
  7  8  1  0  0  0  0
M  END
$$$$
"""
if __name__ == "__main__":
   import StringIO
   from frowns import MDL
   mol, error, text = MDL.sdin(StringIO.StringIO(__test__)).next()

   print mol
#   app = wxPySimpleApp()
   d = Drawer()
   print d.makeGif(mol)
