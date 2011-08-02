from frowns import MDL
from frowns.Depict.wxMoleculeDrawer import MoleculeDrawer

for mol, error, text in MDL.sdin(open("data\chiral.sdf")):
    for bond in mol.bonds:
        print bond.stereo
    print mol, error
    print mol.cansmiles()

from wxPython.wx import *
class testApp(wxApp):
    def OnInit(self):
        from frowns import MDL
        reader = MDL.sdin(open("data\chiral3.sdf"))
        
        mol, error, text = reader.next()

        wxInitAllImageHandlers()
        frame = wxFrame(None, -1, "", size=(350,200))

        
        drawer = self.drawer = MoleculeDrawer(frame, drawAromatic=1)
        drawer.setMolecule(mol)
        frame.Show(TRUE)
        return TRUE

 


      

if __name__ == "__main__":
    _app = testApp(0)
    _app.MainLoop()
