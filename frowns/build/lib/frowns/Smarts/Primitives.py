"""
SMARTS Primitives taken from daylight's web site
http://www.daylight.com

SMARTS Atomic primitives

Symbol   Symbol name            Atomic property            Default
                                requirements
======   ===========            ============               ===========                                
*        wildcard               any atom                   (no default) 
a        aromatic               aromatic                   (no default) 
A        aliphatic              aliphatic                  (no default) 
D<n>     degree                 <n> explicit connections   (no default) 
H<n>     total-H-count          <n> attached hydrogens     exactly one 
h<n>     implicit-H-count       <n> implicit hydrogens     exactly one 
R<n>     ring membership in     <n> SSSR rings             any ring atom 
r<n>     ring size              in smallest SSSR ring      any ring atom 
                                of size <n>   
v<n>     valence                total bond order <n>       (no default) 
X<n>     connectivity <n>       total connections          (no default) 
-<n>     negative charge -<n>   charge                     -1 charge (-- is -2, etc) 
+<n>     positive charge +<n>   formal charge              +1 charge (++ is +2, etc) 
#n       atomic number          atomic number <n>          (no default) 
@        chirality              anticlockwise              anticlockwise, default class 
@@       chirality              clockwise                  clockwise, default class 
@<c><n>  chirality              chiral class <c>           chirality <n> (nodefault) 
@<c><n>? chiral or              <c><n> orunspecified       (no default) 
         unspec chirality 
<n>      atomic mass            explicit atomic mass       unspecified mass

SMARTS Bond Primitives
Symbol   Atomic property           requirements 
======   ===============           ============
-        single bond               (aliphatic) 
/        directional single bond   "up" 
\        directional single bond   "down" 
/?       directional bond          "up or unspecified" 
\?       directional bond          "down or unspecified" 
=        double bond 
#        triple bond 
:        aromatic bond 
~        any bond (wildcard) 
@        any ring bond
"""
import string
#
# The goal of a matcher is to build an object that can be
# used to indicate wheter an atom or bond in a particular
# context is equivelent to a compiled expression.
#
# A compiled expression is an instance that has an overloaded
# equivalence operator (__eq__) so that the following operations
# can be performed
#
# expression == atom
# expression == bond
#

class PropertyMatch:
    """Matches a named property of an object exactly"""
    property = None
    def __init__(self, value):
        self.value = value
    def __eq__(self, atom):
        return self.value == getattr(atom, self.property)
    def __str__(self):
        return "%s == %s" % (self.property, self.value)

class BooleanMatch(PropertyMatch):
    def __init__(self, value = 1):
        assert value in (0, 1)
        PropertyMatch.__init__(self, value)
    def __str__(self):
        if self.value :
            return "%s is true" % (self.property, )
        else:
            return "%s is false" % (self.property, )

class SymbolMatch(PropertyMatch):
    property = "symbol"
    def __eq__(self, atom):
        if self.value == '*': return 1
        return self.value == getattr(atom, self.property)

class AromaticMatch(PropertyMatch):
    property = "aromatic"

class AtomicNumberMatch(PropertyMatch):
    property = "number"

class WeightMatch(PropertyMatch):
    property = "weight"

class ChargeMatch(PropertyMatch):
    property = "charge"

class TotalHMatch(PropertyMatch):
    property = "hcount"

class ImplicitHMatch(PropertyMatch):
    property = "imp_hcount"

class DegreeMatch(PropertyMatch):
    def __eq__(self, atom):
        # we have to perform a computation here...
        return self.value == len(atom.bonds)

class RingMembershipMatch(PropertyMatch):
    def __eq__(self, atom):
        # are we in self.value SSSR rings?
        return self.value == len(atom.rings)
    def __str__(self):
        return "atom in %s rings"%(self.value)
    
class BooleanRingMembershipMatch(BooleanMatch):
    def __eq__(self, atom):
        return len(atom.rings) == 1
    def __str__(self):
        return "atom in 1 ring"

class RingSizeMatch(PropertyMatch):
    # XXX FIX ME
    # I don't understand this one too well.  The description says
    # atom in smallest SSSR ring of size <n>
    # although that is a bit ambiguous isn't it?
    # let's just say for now that the atom is in a ring of
    # size n and leave it at that...
    def __eq__(self, atom):
        for cycle in atom.rings:
            if len(cycle) == self.value:
                return 1
        return 0

class ValenceMatch(PropertyMatch):
    property = "valence"
    def __eq__(self, atom):
        order = 0
        for bond in atom.bonds:
            order += bond.bondorder
        return self.value == order

class ConnectivityMatch(PropertyMatch):
    property = "connectivity"
    def __eq__(self, atom):
        return self.value == atom.hcount + len(atom.bonds)

####################################################
# We don't support chirality yet
# This will be kind of bizzare...
# we will do a small internal match inside
# of this object
class ChiralClassMatch(PropertyMatch):
    property = "XXX ChiralClass"

class ChiralCountMatch(PropertyMatch):
    property = "XXX ChiralCount"

####################################################
# The recursive matcher is not implemented at this
# point.  I probably will need to modified the
# graph searcher a little bit to anchor the search
# at the starting atom of the recursion
class RecursiveMatcher:
    def __init__(self, mol):
        self.mol = mol
    def __eq__(self, atom):
        raise NotImplementedError
    def __str__(self):
        s = self.mol.dump()
        lines = string.split(s, "\n")
        new_lines = []
        for line in lines:
            new_lines.append("recursive> " + line)
        s = string.join(new_lines, "\n")
        return "atom in\n%s\n" % (s,)

class BondSymbolMatch(PropertyMatch):
    property = "symbol"

class UnspecifiedBondMatch(PropertyMatch):
    def __init__(self):
        self.value = "-:"
    
    def __eq__(self, bond):
        if bond.symbol in ["-", ":"]:
            return 1
        return 0

class AnyBond(PropertyMatch):
    property = "anybond"
    def __eq__(self, bond):
        return 1
    
class AnyRingBond(PropertyMatch):
    property = "anybond"
    def __eq__(self, bond):
        a1, a2 = bond.atoms
        if a1.rings and a2.rings:
            return 1
        return 0

