import MoleculeDrawerNew
from wxPython.wx import *

class wxMoleculeRendererHelper:
    def SetXY(self, x, y):
        self.x0, self.y0 = x, y
        
    def clear(self):
        pass       

    def drawOval(self, x, y, xh, yh):
        self.dc.DrawEllipse(self.x0 + x, self.y0 + y, xh-x, yh-y)
        
    def drawLine(self, x1, y1, x2, y2):
        self.dc.DrawLine(self.x0 + x1, self.y0 +y1,
                         self.x0 + x2, self.y0 + y2)
        
    def drawText(self, text, font, fontsize, x, y, color, bg="white"):
        dc = self.dc
        font = wxFont(fontsize, wxDEFAULT, wxNORMAL, wxNORMAL,
                      0, self._textFont)
        dc.SetFont(font)
        w,h = dc.GetTextExtent(text)
        w+=2
        h+=2
        x = x - w/2 + self.x0
        y = y - h/2 + self.y0

        pen = dc.GetPen()
        dc.SetPen(wxWHITE_PEN)
        dc.DrawRectangle(x, y, w, h)
        dc.SetPen(pen)
        dc.DrawText(text, x, y)

class wxMoleculeDrawer(wxScrolledWindow, wxMoleculeRendererHelper):
    def __init__(self, parent, drawAromatic=1):
        size = parent.GetClientSize()
        self.parent = parent
        wxScrolledWindow.__init__(self, parent, -1, (0,0), size,
                                  wxSUNKEN_BORDER)

        self._textFont = "HELVETICA"
        self.SetBackgroundColour(wxNamedColor("WHITE"))
        self.drawer = MoleculeDrawerNew.Drawer(drawAromatic)
        self.SetXY(0,0)
        
        # flag used to redraw only when necessary
        EVT_PAINT(self, self.OnPaint)
        EVT_SIZE(self, self.OnSize)

    def setMolecule(self, mol):
        self.drawer.setMolecule(mol)
        #self.OnSize(None)
        
    def OnSize(self, event):
        self.Refresh()
    
    def OnPaint(self, event):
        self.width, self.height = self.GetSize()
        dc = self.dc = wxPaintDC(self)
        self.PrepareDC(dc)
        dc.BeginDrawing()
        self.drawer.draw(self)
        dc.EndDrawing()
        # we need to clean up the dc otherwise terror ensues!
        self.dc = None

        

class testApp(wxApp):
    def OnInit(self):
        from frowns import MDL
        reader = MDL.mdlin(open("../test/data/bad.sdf"))
        
        mol = reader.next()
        wxInitAllImageHandlers()
        frame = wxFrame(None, -1, "", size=(350,200))

        
        drawer = self.drawer = wxMoleculeDrawer(frame, drawAromatic=1)
        drawer.setMolecule(mol)
        frame.Show(TRUE)
        return TRUE

 


      

if __name__ == "__main__":
    _app = testApp(0)
    _app.MainLoop()
