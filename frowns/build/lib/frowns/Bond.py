"""Bond class
Bonds """
from IdGenerator import defaultGenerator
import weakref

#DAYLIGHT CIS/TRANS support
DX_CHI_CIS=1
DX_CHI_TRANS=2
DX_CHI_NO_DBO=3

#Chirality Support
UP = "UP"
DOWN = "DOWN"
EITHER = "EITHER"

class Bond(object):
    __slots__ = ["adjunct", "atoms", "index", "rings", "_closure",
                 "dbo", "stereo", "parent", "handle",
                 "symbol", "aromatic", "bondorder", "equiv_class",
                 "bondtype", "fixed", "_traverseOrder", "_line"]
    
    def __init__(self,
                 symbol="-",
                 bondorder=1,
                 bondtype=1,
                 fixed=0,
                 stereo=None,                
                 generator=defaultGenerator):
        self.reset(symbol, bondorder, bondtype, fixed, stereo)

        self.adjunct = 0
        self.atoms = []

        self.index = -1
        self.rings = []
        self._closure = 0
        self.dbo = None
        self.stereo = None
        self.parent = None
        
        self.handle = id(self) #generator() # generate a new unique
                               #   # handle
           
    def setdbo(self, bond1, bond2, dboval):
        """Set the double bond orientation for bond1 and bond2
        based on this bond"""
        # this bond must be a double bond
        if self.bondtype != 2:
            raise FrownsError("To set double bond order, center bond must be double!")
        assert dboval in [DX_CHI_CIS, DX_CHI_TRANS, DX_CHI_NO_DBO], "bad dboval value"
        
        self.dbo.append(bond1, bond2, dboval)
        
            
    def reset(self, symbol, bondorder, bondtype, fixed, stereo):
        self.symbol = symbol
        self.bondorder = bondorder
        self.equiv_class = self.bondtype = bondtype
        self.fixed = fixed
        self.stereo = stereo
        # set the aromaticity correctly
        if self.bondtype == 4:
            self.aromatic = 1
        else:
            assert self.bondtype == self.bondorder
            self.aromatic = 0
        
    def set_symbol(self, symbol):
        """(symbol, bondorder) -> set the bondsymbol
        of the molecule"""
        raise "Deprecated"
        self.symbol, self.bondtype, bondorder, self.equiv_class = \
                     BONDLOOKUP[symbol]
        if self.bondtype == 4:
            self.aromatic = 1
        else:
            self.aromatic = 0


    def set_parent(self, parent):
        self.parent = weakref.ref(parent)

    def destroy(self):
        self.parent = None
        if self.atoms:
            a1, a2 = self.atoms
            a1.bonds.remove(self)
            a2.bonds.remove(self)
            a1.oatoms.remove(a2)
            a2.oatoms.remove(a1)
        self.rings = []
                
    def xatom(self, atom):
        """(atom)->return the atom at the other end of this bond
        or None if atom is not part of this bond"""
        handle = atom.handle
        
        if handle == self.atoms[0].handle:
            return self.atoms[1]
        elif handle == self.atoms[1].handle:
            return self.atoms[0]
        return None
        

    def __repr__(self):
        return "%s(%s)"%(self.__class__.__name__,
                         self.index)

