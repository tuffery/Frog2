from frowns.Objects import FrownTypes
from IdGenerator import defaultGenerator
from Chirality import TetraHedral
from FrownsError import FrownsError

defaultAtomTypes = FrownTypes.defaultAtomTypes
#
# Things might change a little bit, but here are some
# implementation specific details
#
# chiral ordering is not defined yet...
#
# self.rings -> what rings is this atom in...
#               I don't know if I need yet
#
# self._findbonds -> lookup table for the findbond function
#                    one of the ring detectors uses this
#                    function frequently, so an optimization
#                    was in order.

CLOCKWISE = "@@"
ANTICLOCKWISE = "@"

class Atom(object):
    __slots__ = ["symbol", "hcount", "explicit_hcount", "imp_hcount",
                 "charge", "weight", "aromatic", "chiral_order", "adjunct",
                 "equiv_class", "symclass", "symorder", "index",
                 "rings", "bonds", "oatoms", "_closure", "chirality",
                 "chiral_class", "x", "y", "z", "parent", "handle",
                 "valences", "number", "mass", "negativity", "name",
                 "_chirality", "_line", "has_explicit_hcount"
                ]
    def __init__(self, generator=defaultGenerator):
        self.symbol = None
        self.hcount = 0
        self.explicit_hcount = 0
        self.has_explicit_hcount = False
        self.imp_hcount = 0
        self.charge = 0
        self.weight = 0
        self.aromatic = 0
        self.chiral_order = None
        self.adjunct = 0
        self.equiv_class  = -1
        self.symclass = -1
        self.symorder = -1
        self.index = -1
        
        self.rings = []
        self.bonds = []
        self.oatoms = []
        self._closure = 0
        self.chirality = None
        self._chirality = None
        self.chiral_class = None

        # coordinates for drawing and the like
        self.x = None
        self.y = None
        self.z = None

        self._line = ""
        
        self.parent = None
        self.handle = id(self)#generator()
        
    def set_symbol(self, symbol):
        self.symbol, self.number, self.name, self.mass, self.negativity, \
                     self.valences, self.equiv_class = defaultAtomTypes[symbol]

## XXX implicit_hcount is the responsibility of the modifing code
##  (actually, everything is...)  maybe I should try
## properties

    def findbond(self, otherAtom):
        handle = otherAtom.handle
        for atom, bond in zip(self.oatoms, self.bonds):
            if handle == atom.handle:
                return bond
        return None

    def sumBondOrders(self):
        result = 0
        for x in self.bonds:
            if x.bondtype == 4:
                # XXX FIX ME
                # this is a hack to fix bad conjugation
                # this will be fixed soon
                result += 1.5
            else:
                result += x.bondorder
        return result

    def destroy(self):
        self.rings = []
        self.bonds = []
        self.oatoms = []

    def chival(self, bonds):
        """compute the chiral value around an atom given a list of bonds"""
        # XXX I'm not sure how this works?
        order = [bond.xatom(self) for bond in bonds]
        return self._chirality(order)
    
    def setchival(self, bondorder, rotation):
        """compute chiral ordering of surrounding atoms"""
        rotation = [None, "@", "@@"][(rotation % 2)]
        # check to see if the bonds are attached
        if not bondorder: # use the default xatoms
            if len(self.oatoms) < 3 and self.explicit_hcount != 1:
                raise FrownsError("Need to have an explicit hydrogen when specifying "\
                                  "chirality with less than three bonds")
                

            self._chirality = TetraHedral.T(self.oatoms,
                                            rotation)
            return
        if len(bondorder) != len(self.bonds):
            raise AtomError("The order of all bonds must be specified")
        
        for bond in bondorder:
            if bond not in self.bonds:
                raise AtomError("Specified bonds to assign chirality are not attatched to atom")

        order = [bond.xatom(self) for bond in bonds]
        self._chirality = TetraHedral.T(order, rotation)
        
    def __repr__(self):
        return "%s(%s)"%(self.__class__.__name__,
                         self.index)

