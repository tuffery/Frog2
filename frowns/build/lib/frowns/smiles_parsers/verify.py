
# Used when building the element matching regular expression string.
# Not to be imported directly.

import string, re

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

aromatic_symbols = list("cnosp")

_special_order = "CONSPHcnosp" + string.letters + "*"

def _special_order_cmp(a, b):
    return cmp(_special_order.index(a), _special_order.index(b))

def make_re_pattern(symbols):
    # Turn the list of element symbols into something that can
    # be used in a regular expression match.  Do some premature
    # optimization so the most common first letters are checked
    # first.
    name_table = {}
    allow_single = {}
    for name in symbols:
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

if __name__ == "__main__":
    import Smiles, Smarts
    err = 0
    expected = make_re_pattern(element_symbols + aromatic_symbols)
    if Smiles.element_symbols_pattern != expected:
        print "In Smiles, expected:"
        print repr(expected)
        print "but the Smiles module has:"
        print repr(Smiles.element_symbols_pattern)
        err = 1

    symbols = element_symbols[:]
    symbols.remove("H")
    expected = make_re_pattern(symbols + aromatic_symbols)
    if Smarts.element_symbols_pattern != expected:
        print "In Smarts, expected:"
        print repr(expected)
        print "but the Smarts module has:"
        print repr(Smarts.element_symbols_pattern)
        err = 1

    assert not err, "Please fix things"
                                                     
else:
    raise ImportError("This module should only be used to verify constants.")
