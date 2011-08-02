"""
atom := [cnosp]
atom := 'B' | 'C' | 'N' | 'O' | 'F' | 'P' | 'S' | 'Cl' |  'Br' |  'I'
atom := '[' mass? symbol chiral? hcount? charge? ']'

symbol := a lot of names ...

hcount := 'H' \d+
mass := \d+

chiral := '@'  (chiral_repeat | chiral_count | chiral_name)
chiral_repeat := '@'+
chiral_count := \d+
chiral_name := TH[12] | AL[12] | SP[123] | TB(1[0-9]?|20?|[3-9]) |
                OH(1[0-9]?|2[0-9]?|30?|[4-9])

charge := "+" ("+"+ | \d+)?
charge := "-" ("-"+ | \d+)?

bond := [=#/\\:~-]

dot := \.

closure = \d | '%'\d\d?

start -> atom
#start -> dot     # not allowed by Daylight

atom -> atom     # CC
atom -> bond     # C=C
atom -> branch_out  # CCC(C)C   (the 'C(')
atom -> branch_in   # CCC(C)C   (the 'C)')
atom -> dot      # C.C
atom -> closure  # c1ccccc1

bond -> atom     # C=C
bond -> closure  # C1CCCCC-1
bond -> branch_out   # CCC=(O)CC  # but no branch_in

branch_out -> atom   # CC(C)C
branch_in -> atom
branch_out -> bond   # CC(=O)C(C)=O
branch_in -> bond
branch_out -> dot    # C(.C).C
branch_in -> dot
branch_in -> branch_out  # CCC(C)(C)CC
branch_in -> branch_in   # CC(C(C))C

dot -> atom      # C.C
#dot -> dot       # not allowed by Daylight

closure -> atom   # c1ccccc1
closure -> bond   # C1.C1#C
closure -> closure # C12C3C4C1C5C4C3C25
closure -> branch # C1.C1(C)C
closure -> dot    # c1ccccc1.O

"""

import re, string

element_symbols = [
    '*', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',   # 10
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',     # 20
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',   # 30
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr',   # 40
    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',  # 50
    'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',   # 60
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',  # 70
    'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',   # 80
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',  # 90
    'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',   # 100
    'Md', 'No', 'Lr', 'Rf', 'Ha', 'Sg', 'Ns', 'Hs', 'Mt',
    ]
# Add in the aromatics
element_symbols.extend(list("cnosp"))

_special_order = "CONSPHcnosp" + string.letters + "*"

def _special_order_cmp(a, b):
    return cmp(_special_order.index(a), _special_order.index(b))

def _make_element_symbols_pattern():
    # Turn the list of element symbols into something that can
    # be used in a regular expression match.  Do some premature
    # optimization so the most common first letters are checked
    # first.
    name_table = {}
    allow_single = {}
    for name in element_symbols:
        if len(name) == 1:
            # Single letter symbol
            name_table.setdefault(name, [])
            allow_single[name] = 1
        else:
            # Two letter symbol
            c1, c2 = name
            name_table.setdefault(c1, []).append(c2)

    # Order so things like "C" are tested early
    keys = name_table.keys()
    keys.sort(_special_order_cmp)

    terms = []
    for key in keys:
        seconds = name_table[key]
        if len(seconds) == 0:
            # Only a single element starts with this letter and
            # the symbol is only one letter long. (eg, "U")
            s = re.escape(key)
        elif len(seconds) == 1:
            assert re.escape(seconds[0]) == seconds[0]
            # Either only one element exists with this letter and
            # the symbol is two letters long (eg, "Xe")
            # -or-
            # Two elements exist, one with a one letter long name
            # and the other with two letters (eg, "O" and "Os")
            s = re.escape(key) + seconds[0]
            if allow_single.get(key):
                s = s + "?"
        else:
            # Several elements start with the same first letter
            s = "%s[%s]" % (re.escape(key), string.join(seconds, ""))
            if allow_single.get(key):
                s = s + "?"
        terms.append(s)
    return string.join(terms, "|")

element_symbols_pattern = _make_element_symbols_pattern()

atom_fields = [
    "raw_atom",
    "open_bracket",
    "weight",
    "element",
    "chiral_count",
    "chiral_named",
    "chiral_simple",
    "hcount",
    "charge_count",
    "charge_simple",
    "error_1",
    "error_2",
    "close_bracket",
    ]

atom = re.compile(r"""
(?P<raw_atom>Cl|Br|[cnospBCNOFPSI]) |  # Don't need brackets
(
  (?P<open_bracket>\[)                 # Start bracket
  (?P<weight>\d+)?                     # Atomic weight (optional)
  (                                    # valid term or error
   (                                   #   valid term
    (?P<element>""" + element_symbols_pattern + r""")  # element
    (                                  # Chirality can be
     (?P<chiral_count>@\d+) |          #   @1 @2 @3 ...
     (?P<chiral_named>                 # or
       @TH[12] |                       #   @TA1 @TA2
       @AL[12] |                       #   @AL1 @AL2
       @SP[123] |                      #   @SP1 @SP2 @SP3
       @TB(1[0-9]?|20?|[3-9]) |        #   @TB{1-20}
       @OH(1[0-9]?|2[0-9]?|30?|[4-9])) | # @OH{1-30}
     (?P<chiral_simple>@+)             # or @@@@@@@...
    )?                                 # and chirality is optional
    (?P<hcount>H\d*)?                  # Optional hydrogen count
    (                                  # Charges can be
     (?P<charge_count>[+-]\d+) |       #   +<number> or -<number>
     (?P<charge_simple>\++|\-+)        #   +++...  or ---
    )?                                 # and are optional
    (?P<error_1>[^\]]+)?               # If there's anything left, it's an error
  ) | (                                # End of parsing stuff in []s, except
    (?P<error_2>[^\]]*)                # If there was an error, we get here
  ))
  (?P<close_bracket>\])                # End bracket
)
""", re.X)

bond_fields = ["bond"]
bond = re.compile(r"(?P<bond>[=#/\\:~-])")

dot_fields = ["dot"]
dot = re.compile(r"(?P<dot>\.)")

closure_fields = ["closure"]
closure = re.compile(r"(?P<closure>\d|%\d\d?)")

branch_in_fields = ["branch_in"]
branch_in = re.compile(r"(?P<branch_in>\))")

branch_out_fields = ["branch_out"]
branch_out = re.compile(r"(?P<branch_out>\()")

# Make from the state name to the regular expession to
# match and the list of field names in the regexp.
# (There's no way to use sre to get that as an ordered list.)
info = {
    "atom": (atom, atom_fields),
    "bond": (bond, bond_fields),
    "dot": (dot, dot_fields),
    "closure": (closure, closure_fields),
    "branch_in": (branch_in, branch_in_fields),
    "branch_out": (branch_out, branch_out_fields),
    }

# Mapping from current state to allowed states
table = {
    "start": ("atom",),
    "atom": ("atom", "bond", "branch_in", "branch_out", "dot", "closure"),
    "bond": ("atom", "closure", "branch_out"),
    "branch_in": ("atom", "bond", "dot", "branch_in", "branch_out"),
    "branch_out": ("atom", "bond", "dot"),
    "dot": ("atom","dot"),
    "closure": ("atom", "bond", "closure", "branch_in", "branch_out", "dot"),
}

# Parse a SMILES string and print the events found
import Generator

class ParseError(Exception):
    pass

def tokenizer(s):
    tokens = []
    expected = table["start"]
    while s:
        for state in expected:
            pat, fields = info[state]
            m = pat.match(s)
            if m:
                break
        else:
            raise ParseError("Could not parse starting with %s" % (repr(s),))
        d = m.groupdict()

        for field in ["error_1", "error_2"]:
            if d.get(field, None) is not None:
                    raise ParseError(
                        "Caught error starting at %s" % (repr(s),)
                        )
        tokens.append((state, d))
        expected = table[state]
        s = s[m.end(0):]

    return tokens

def molecule_generator(tokens, generator=Generator.Generator):
    # quick scans for closure tokens
    # it's easier for us to do the bond checks here
    # than in the molecule generator.  Hopefully this
    # will save us some time.

    mol = generator()
    for state, dict in tokens:
        mol.addState(state, dict)
    return mol.mol()
