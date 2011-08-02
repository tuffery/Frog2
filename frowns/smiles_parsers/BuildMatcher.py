"""
SMARTS Primitives taken from daylight's web site
http://www.daylight.com

SMARTS Atomic primitives

Symbol   Symbol name            Atomic property            Default
                                requirements
======   ===========            ============                ===========                                
*        wildcard               any atom                   (no default) 
a        aromatic               aromatic                   (no default) 
A        aliphatic              aliphatic                  (no default) 
D<n>     degree <n>             explicit connections       (no default) 
H<n>     total-H-count <n>      attached hydrogens         exactly one 
h<n>     implicit-H-count <n>   implicit hydrogens         exactly one 
R<n>     ring membership in <n> SSSR rings                 any ring atom 
r<n>     ring size              in smallest SSSR ring      any ring atom 
                                of size <n>   
v<n>     valence                total bond order <n>       (no default) 
X<n>     connectivity <n>       total connections          (no default) 
- <n>    negative charge -<n>   charge                     -1 charge (-- is -2, etc) 
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

# Build a simple Molecule Matcher object given the events from the
# Smarts tokenizer.

# Notes:
#  - A lot of optimization could be done on the graph
#  - The "match" functions don't work and shouldn't be done this way
#  - I think the recursive SMARTS need to be done first, once, and
#      the results used for the higher-level matches
#  - I'm not sure about the zero-level SMARTS, but what I have does
#      match the Daylight toolkit behaviour
#  - This is a first-draft, proof of principle parser.

import string
import Handler

class PropertyMatch:
    property = None
    def __init__(self, value):
        self.value = value
    def match(self, atom):
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

class AromaticMatch(PropertyMatch):
    property = "aromatic"

class AtomicNumberMatch(PropertyMatch):
    property = "number"

class WeightMatch(PropertyMatch):
    property = "weight"

class ChargeMatch(PropertyMatch):
    property = "charge"

class TotalHMatch(PropertyMatch):
    property = "explicit_hcount"

class ImplicitHMatch(PropertyMatch):
    property = "imp_hcount"

class DegreeMatch(PropertyMatch):
    property = "degree"

class RingMembershipMatch(PropertyMatch):
    property = "XXX RingMembership"

class BooleanRingMembershipMatch(BooleanMatch):
    property = "XXX2 BooleanRingMembership"

class RingSizeMatch(PropertyMatch):
    property = "XXX RingSize"

class ValenceMatch(PropertyMatch):
    property = "XXX Valence"

class ConnectivityMatch(PropertyMatch):
    property = "XXX Connectivity"

class ChiralClassMatch(PropertyMatch):
    property = "XXX ChiralClass"

class ChiralCountMatch(PropertyMatch):
    property = "XXX ChiralCount"

class RecursiveMatcher:
    def __init__(self, mol):
        self.mol = mol
    def match(self, atom):
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

class NotMatch:
    def __init__(self, child):
        self.child = child
    def match(self, obj):
        return self.child.match(obj)
    def __str__(self):
        return "not (%s)" % (self.child,)

class AndMatch:
    def __init__(self, left, right):
        self.left = left
        self.right = right
    def match(self, obj):
        return self.left.match(obj) and self.right.match(obj)
    def __str__(self):
        return "AND(%s, %s)" % (self.left, self.right)
    
class OrMatch:
    def __init__(self, left, right):
        self.left = left
        self.right = right
    def match(self, obj):
        return self.left.match(obj) or self.right.match(obj)
    def __str__(self):
        return "OR(%s, %s)" % (self.left, self.right)

bool_unary_not = 76
bool_strong_and = 77
bool_or = 78
bool_weak_and = 79
binary_operators = [bool_strong_and, bool_or, bool_weak_and]
boolean_operators = binary_operators + [bool_unary_not]
text_to_bool = {
    "&": bool_strong_and,
    ",": bool_or,
    ";": bool_weak_and,
    "!": bool_unary_not,
    }

class ExpressionList:
    def __init__(self):
        self.matchers = []
    def __nonzero__(self):
        return len(self.matchers) != 0
    def add_matcher(self, obj):
        if self.matchers and self.matchers[-1] not in boolean_operators:
            self.matchers.append(bool_strong_and)
        self.matchers.append(obj)
    def add_operator(self, op):
        assert op in binary_operators or op == bool_unary_not
        if __debug__:
            if self.matchers:
                if op in binary_operators:
                    assert self.matchers[-1] not in binary_operators
            else:
                assert op not in binary_operators
        self.matchers.append(op)
    def make_matcher(self):
        matchers = self.matchers[:]
        i = 0
        while i < len(matchers):
            if matchers[i] == bool_unary_not:
                matchers[i:i+2] = [NotMatch(matchers[i+1])]
            else:
                i = i + 1
        i = 1
        while i < len(matchers):
            if matchers[i] == bool_strong_and:
                matchers[i-1:i+2] = [AndMatch(matchers[i-1], matchers[i+1])]
            else:
                i = i + 1
        i = 1
        while i < len(matchers):
            if matchers[i] == bool_or:
                matchers[i-1:i+2] = [OrMatch(matchers[i-1], matchers[i+1])]
            else:
                i = i + 1
        i = 1
        while i < len(matchers):
            if matchers[i] == bool_weak_and:
                matchers[i-1:i+2] = [AndMatch(matchers[i-1], matchers[i+1])]
            else:
                i = i + 1
        assert len(matchers) == 1, matchers
        return matchers[0]

class AtomExpression(ExpressionList):
    pass

class BondExpression(ExpressionList):
    pass
        

class Atom:
    def __init__(self, matcher, component_number):
        self.matcher = matcher
        self.component_number = component_number
        self.bonds = []
    def match(self, atom):
        return self.matcher(atom)
    def dump_info(self, bonds):
        text = "[%d] <%s>" % (self.component_number, str(self.matcher))
        text = text + " bonds = ["
        bond_ids = []
        for bond in self.bonds:
            bond_ids.append(str(bonds.index(bond)))
        text = text + string.join(bond_ids, ", ") + "]"
        return text


class Bond:
    def __init__(self, matcher, component_number):
        self.matcher = matcher
        self.component_number = component_number
        self.atoms = []
    def match(self, bond):
        return self.matcher(bond)
    def dump_info(self, atoms):
        return "[%d] %d <%s> %d" % (self.component_number,
                                    atoms.index(self.atoms[0]),
                                    self.matcher,
                                    atoms.index(self.atoms[1]))
        

class Molecule:
    def __init__(self):
        self.atoms = []
        self.bonds = []
    def dump(self):
        lines = ["atom_id [component_number] <matcher> "
                 "bonds = [list of bond ids]"]
        i = 0
        for atom in self.atoms:
            lines.append(str(i) + " " + atom.dump_info(self.bonds))
            i = i + 1
        lines.append("bond_id [component_number> "
                     "atom1_id <matcher> atom2_id")
        i = 0
        for bond in self.bonds:
            lines.append(str(i) + " " + bond.dump_info(self.atoms))
            i = i + 1
        return string.join(lines, "\n")
    def __del__(self):
        for atom in self.atoms:
            del atom.bonds
        for bond in self.bonds:
            del bond.atoms

def get_symbol_aromatic(text):
    if text[0] in "cnosp":
        return string.upper(text), 1
    return text, 0

def normalize_closure(text):
    if text[:1] == "%":
        return int(text[1:])
    return int(text)

class BuildMatcher(Handler.TokenHandler):
    _save_state = ("closures", "mol", "_atom_expr", "_prev_atoms",
                   "_pending_bond_expr")
    def begin(self):
        self._reset()
        self.component_number = 0
        self.recursive_smarts_stack = []
        
    def _reset(self):
        self.closures = {}
        self.mol = Molecule()
        self._atom_expr = None
        self._prev_atoms = []

        # None occurs after a '.'
        self._pending_bond_expr = None

    def _verify(self):
        if len(self._prev_atoms) >= 2:
            raise AssertionError("Missing ')'")
        if self._pending_bond_expr:
            raise AssertionError("Missing an atom expression after the bond")
        if self.closures:
            raise AssertionError("Missing closures for %s" %
                                 (self.closures.keys(),))
    def end(self):
        self._verify()
        if self.recursive_smarts_stack:
            raise AssertionError("Inside of a recursive SMARTS")


    def add_token(self, field, pos, text):
        getattr(self, "do_" + field)(text)

    def do_open_recursive_smarts(self, text):
        d = {}
        for k in self._save_state:
            d[k] = getattr(self, k)
        self.recursive_smarts_stack.append(d)
        self._reset()
        
    def do_close_recursive_smarts(self, text):
        self._verify()
        mol = self.mol
        self._reset()
        d = self.recursive_smarts_stack.pop()
        for k, v in d.items():
            setattr(self, k, v)
        self._atom_expr.add_matcher(RecursiveMatcher(mol))

    def add_atom(self, atom):
        if self._pending_bond_expr is not None and not self._pending_bond_expr:
            # Implicit single or aromatic bond
            self._pending_bond_expr = BondExpression()
            self._pending_bond_expr.add_matcher(BondSymbolMatch("-"))

        if self._pending_bond_expr is not None:
            bond = Bond(self._pending_bond_expr.make_matcher(),
                        self.component_number)
            prev_atom = self._prev_atoms[-1]
            bond.atoms[:] = [prev_atom, atom]
            prev_atom.bonds.append(bond)
            atom.bonds.append(bond)
            self.mol.bonds.append(bond)
        self._pending_bond_expr = BondExpression()
        if not self._prev_atoms:
            self._prev_atoms.append(atom)
        else:
            self._prev_atoms[-1] = atom
        self.mol.atoms.append(atom)
        
    def do_raw_atom(self, text):
        symbol, aromatic = get_symbol_aromatic(text)
        atom_match = AndMatch(SymbolMatch(symbol), AromaticMatch(aromatic))
        self.add_atom(Atom(atom_match,
                           self.component_number))

    def do_raw_aromatic(self, text):
        self.add_atom(Atom(AromaticMatch(1),
                           self.component_number))
    def do_raw_aliphatic(self, text):
        self.add_atom(Atom(AliphaticMatch(0),
                           self.component_number))
    def do_raw_b_unknown(self, text):
        1/0
    def do_raw_f_unknown(self, text):
        1/0
    def do_raw_h_unknown(self, text):
        1/0
    def do_raw_i_unknown(self, text):
        1/0
    def do_raw_r_unknown(self, text):
        # I think this is right
        self.add_atom(Atom(BooleanRingMembershipMatch(),
                           self.component_number))
    def do_raw_R_unknown(self, text):
        # I think this is right
        self.add_atom(Atom(BooleanRingMembershipMatch(),
                           self.component_number))

    def do_open_bracket(self, text):
        self._atom_expr = AtomExpression()

    def do_atom_not(self, text):
        self._atom_expr.add_operator(text_to_bool[text])
    def do_atom_binary(self, text):
        self._atom_expr.add_operator(text_to_bool[text])

    def do_atomic_number(self, text):
        self._atom_expr.add_matcher(AtomicNumberMatch(text[1:]))
    def do_weight(self, text):
        self._atom_expr.add_matcher(WeightMatch(text))
    def do_element(self, text):
        symbol, aromatic = get_symbol_aromatic(text)
        atom_matcher = AndMatch(SymbolMatch(symbol),
                                AromaticMatch(aromatic))
        self._atom_expr.add_matcher(atom_matcher)
    def do_chiral_count(self, text):
        self._atom_expr.add_matcher(ChiralCountMatch(text[1:]))
    def do_chiral_named(self, text):
        self._atom_expr.add_matcher(ChiralClassMatch(text[1:3]))
        self._atom_expr.add_matcher(ChiralCountMatch(int(text[3:])))
    def do_chiral_symbols(self, text):
        self._atom_expr.add_matcher(ChiralCountMatch(len(text)))

    def do_aromatic(self, text):
        self._atom_expr.add_matcher(AromaticMatch(1))
    def do_aliphatic(self, text):
        self._atom_expr.add_matcher(AromaticMatch(0))
    def do_total_hcount(self, text):
        if text == "H":
            count = 1
        else:
            count = int(text[1:])
        self._atom_expr.add_matcher(TotalHMatch(count))
    def do_imp_hcount(self, text):
        self._atom_expr.add_matcher(ImplicitHMatch(text[1:]))
    def do_degree(self, text):
        self._atom_expr.add_matcher(DegreeMatch(text[1:]))
    def do_ring_membership(self, text):
        if text == "R":
            self._atom_expr.add_matcher(BooleanRingMembershipMatch())
        else:
            self._atom_expr.add_matcher(RingMembershipMatch(text[1:]))
    def do_ring_size(self, text):
        if text == "r":
            self._atom_expr.add_matcher(BooleanRingMembershipMatch())
        else:
            self._atom_expr.add_matcher(RingSizeMatch(text[1:]))
    def do_valence(self, text):
        self._atom_expr.add_matcher(ValenceMatch(text[1:]))
    def do_connectivity(self, text):
        self._atom_expr.add_matcher(ConnectivityMatch(text[1:]))

    def do_positive_count(self, text):
        self._atom_expr.add_matcher(ChargeMatch(int(text)))
    def do_positive_symbols(self, text):
        self._atom_expr.add_matcher(ChargeMatch(len(text)))
    def do_negative_count(self, text):
        self._atom_expr.add_matcher(ChargeMatch(int(text)))
    def do_negative_symbols(self, text):
        self._atom_expr.add_matcher(ChargeMatch(-len(text)))

    def do_close_bracket(self, text):
        self.add_atom(Atom(self._atom_expr.make_matcher(),
                           self.component_number))
        self._atom_expr = None


    def do_bond(self, text):
        self._pending_bond_expr.add_matcher(BondSymbolMatch(text))
    def do_bond_not(self, text):
        self._pending_bond_expr.add_operator(text_to_bool[text])
    def do_bond_binary(self, text):
        self._pending_bond_expr.add_operator(text_to_bool[text])

    def do_dot(self, text):
        assert not self._pending_bond_expr, "not possible"
        self._pending_bond_expr = None

    def do_closure(self, text):
        num = normalize_closure(text)
        if self.closures.has_key(num):
            prev_atom, prev_bond_expr = self.closures[num]
            del self.closures[num]

            # Because things like ".1" or "1C" can't occur
            assert self._pending_bond_expr is not None, "Can't happen"
            assert prev_bond_expr is not None, "Can't happen either"
            
            if not prev_bond_expr:
                if not self._pending_bond_expr:
                    bond_matcher = BondSymbolMatch("-")
                else:
                    bond_matcher = self._pending_bond_expr.make_matcher()
            else:
                if not self._pending_bond_expr:
                    bond_matcher = prev_bond_expr.make_matcher()
                else:
                    raise NotImplementedError("Need to check if they match")

            bond = Bond(bond_matcher, self.component_number)

            atom = self._prev_atoms[-1]
            if prev_atom is atom:
                raise AssertionError("cannot close a ring with itself")
            bond.atoms[:] = [prev_atom, atom]
            prev_atom.bonds.append(bond)
            atom.bonds.append(bond)
            self.mol.bonds.append(bond)
        else:
            self.closures[num] = (self._prev_atoms[-1], self._pending_bond_expr)
        self._pending_bond_expr = BondExpression()

    def do_open_branch(self, text):
        self._prev_atoms.append(self._prev_atoms[-1])
    
    def do_close_branch(self, text):
        self._prev_atoms.pop()

    def do_open_zero(self, text):
        # I think this component thing is right
        pass
    
    def do_close_zero(self, text):
        self._verify()
        self.component_number = self.component_number + 1
    
def test():
    import Smarts
    h = BuildMatcher()
    for smi in ["C",
                "CC",
                "c",
                "C(N)O",
                "[O]",
                "c1ccccc1",
                ]:
        print "*"*44
        print "-->", smi        
        Smarts.tokenize(smi, h)
        print h.mol.dump()

if __name__ == "__main__":
    test()
