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

from frowns.smiles_parsers import Smarts
from frowns.smiles_parsers import Handler
from SmartsMatcher import Matcher
from Pattern import *
from Expressions import *
from Primitives import *
import string

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
        self.mol = SmartsPattern()
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
            self._pending_bond_expr.add_matcher(UnspecifiedBondMatch())
                #BondSymbolMatch("-"))

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
        if text == "*":
            atom_match = SymbolMatch(text)
        else:
            symbol, aromatic = get_symbol_aromatic(text)
            atom_match = AndMatch(SymbolMatch(symbol), AromaticMatch(aromatic))

        self.add_atom(Atom(atom_match,
                           self.component_number))

    def do_raw_aromatic(self, text):
        self.add_atom(Atom(AromaticMatch(1),
                           self.component_number))

    def do_raw_aliphatic(self, text):
        self.add_atom(Atom(AromaticMatch(0),
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
        self._atom_expr.add_matcher(AtomicNumberMatch(int(text[1:])))

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
        self._atom_expr.add_matcher(ImplicitHMatch(int(text[1:])))

    def do_degree(self, text):
        self._atom_expr.add_matcher(DegreeMatch(int(text[1:])))

    def do_ring_membership(self, text):
        if text == "R":
            self._atom_expr.add_matcher(BooleanRingMembershipMatch())
        else:
            self._atom_expr.add_matcher(RingMembershipMatch(int(text[1:])))

    def do_ring_size(self, text):
        if text == "r":
            self._atom_expr.add_matcher(BooleanRingMembershipMatch())
        else:
            self._atom_expr.add_matcher(RingSizeMatch(int(text[1:])))

    def do_valence(self, text):
        self._atom_expr.add_matcher(ValenceMatch(text[1:]))

    def do_connectivity(self, text):
        self._atom_expr.add_matcher(ConnectivityMatch(int(text[1:])))


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
        if text == "~":
            self._pending_bond_expr.add_matcher(AnyBond(text))
        elif text == "@":
            self._pending_bond_expr.add_matcher(AnyRingBond(text))
        elif text == "-":
            # explicit single bonds can't be aromatic
            match = AndMatch(AromaticMatch(0), BondSymbolMatch(text))
            self._pending_bond_expr.add_matcher(match)
        else:
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
                    bond_matcher = UnspecifiedBondMatch()#BondSymbolMatch("-")
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

def compile(smarts_pattern):
    h = BuildMatcher()
    Smarts.tokenize(smarts_pattern, h)
    return Matcher(h.mol)

def test():
    from frowns import Smiles
    import time
    h = BuildMatcher()
    for smi in ["C",
                "CC",
                "c",
                "C(N)O",
                "[O]",
                "c1ccccc1",
                ]:
        matcher = compile(smi)
        print matcher.dump()
        smiles = Smiles.smilin(smi)
        pathset = matcher.match(smiles)


if __name__ == "__main__":
    test()
