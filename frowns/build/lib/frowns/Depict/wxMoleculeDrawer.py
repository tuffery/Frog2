import MoleculeDrawer
from wxPython.wx import *

class wxMixin(wxScrolledWindow):
    def _init(self, parent, height=None, width=None):
        self.x0 = 0
        self.y0 = 0
        self.name = ""
        if height is None and width is None:
            size = parent.GetClientSize()
        else:
            size = (width, height)
        self.parent = parent
        wxScrolledWindow.__init__(self, parent, -1, (0,0), size,
                                  wxSUNKEN_BORDER)

        self._textFont = "HELVETICA"
        self.SetBackgroundColour(wxNamedColor("WHITE"))
        
        # flag used to redraw only when necessary
        EVT_PAINT(self, self.OnPaint)
        EVT_SIZE(self, self.OnSize)

    def _clear(self):        
        return 1
    
    def OnSize(self, event):
        self.Refresh()
        
    def OnPaint(self, event, dc=None):
        self.width, self.height = self.GetSize()

        if dc is None:
            dc = self.dc = wxPaintDC(self)
        else:
            self.dc = dc
            
        self.PrepareDC(dc)
        dc.BeginDrawing()

        if self.name:
            fontsize = min(16, self.width/len(self.name))
            self._drawText(self.name,  "Arial", fontsize, self.width/2, self.height-2*fontsize, "red")
            
        self.draw()

        dc.EndDrawing()
        # we need to clean up the dc otherwise terror ensues!
        self.dc = None

    def _drawPoly(self, coords, color="black"):
        dc = self.dc
        pen = dc.GetPen()
        brush = dc.GetBrush()
        if color is not None:
            dc.SetPen(wxPen(color, 1, wxSOLID))
            dc.SetBrush(wxBrush(color, wxSOLID))
        self.dc.DrawPolygon(coords)
        dc.SetPen(pen)
        dc.SetBrush(brush)
        
    def _drawOval(self, x, y, xh, yh):
        self.dc.DrawEllipse(self.x0+x, self.y0+y, xh-x, yh-y)
        
    def _drawLine(self, x1, y1, x2, y2, color):
        color = wxTheColourDatabase.FindColour(color)
        dc = self.dc
        pen = dc.GetPen()
        if color is not None:
            dc.SetPen(wxPen(color, 1, wxSOLID))
            
        self.dc.DrawLine(self.x0+x1, self.y0+y1, self.x0+x2, self.y0+y2)
        dc.SetPen(pen)
        
    def _drawText(self, text, font, fontsize, x, y, color, bg="white"):
        dc = self.dc
        font = wxFont(fontsize, wxDEFAULT, wxNORMAL, wxNORMAL,
                      0, self._textFont)
        dc.SetFont(font)
        w,h = dc.GetTextExtent(text)
        w+=2
        h+=2
        x -= w/2
        y -= h/2

        pen = dc.GetPen()
        # hack        
        dc.SetPen(wxWHITE_PEN)
        dc.DrawRectangle(self.x0+x, self.y0+y, w, h)
        
        dc.DrawText(text, self.x0+x, self.y0+y)
        dc.SetPen(pen)

class MoleculeDrawer(MoleculeDrawer.DrawMolHarness, wxMixin):
    pass
