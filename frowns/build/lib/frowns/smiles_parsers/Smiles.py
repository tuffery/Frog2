import re, string

import Handler

# To verify this is correct, run
#  support.make_re_pattern(support.element_symbols + support.aromatic_symbols)
# Yes, I included aromatics in this list.
element_symbols_pattern = \
  r"C[laroudsemf]?|Os?|N[eaibdpos]?|S[icernbmg]?|P[drmtboau]?|"  \
  r"H[eofgas]?|c|n|o|s|p|A[lrsgutcm]|B[eraik]?|Dy|E[urs]|F[erm]?|"  \
  r"G[aed]|I[nr]?|Kr?|L[iaur]|M[gnodt]|R[buhenaf]|T[icebmalh]|" \
  r"U|V|W|Xe|Yb?|Z[nr]|\*"

atom_fields = [
    "raw_atom",
    "open_bracket",
    "weight",
    "element",
    "chiral_count",
    "chiral_named",
    "chiral_symbols",
    "hcount",
    "positive_count",
    "positive_symbols",
    "negative_count",
    "negative_symbols",
    "error_1",
    "error_2",
    "close_bracket",
    "error_3",
    ]

atom = re.compile(r"""
(?P<raw_atom>Cl|Br|[cnospBCNOFPSI]) |  # "raw" means outside of brackets
(
  (?P<open_bracket>\[)                 # Start bracket
  (?P<weight>\d+)?                     # Atomic weight (optional)
  (                                    # valid term or error
   (                                   #   valid term
    (?P<element>""" + element_symbols_pattern + r""")  # element or aromatic
    (                                  # Chirality can be
     (?P<chiral_count>@\d+) |          #   @1 @2 @3 ...
     (?P<chiral_named>                 # or
       @TH[12] |                       #   @TA1 @TA2
       @AL[12] |                       #   @AL1 @AL2
       @SP[123] |                      #   @SP1 @SP2 @SP3
       @TB(1[0-9]?|20?|[3-9]) |        #   @TB{1-20}
       @OH(1[0-9]?|2[0-9]?|30?|[4-9])) | # @OH{1-30}
     (?P<chiral_symbols>@+)            # or @@@@@@@...
    )?                                 # and chirality is optional
    (?P<hcount>H\d*)?                  # Optional hydrogen count
    (                                  # Charges can be
     (?P<positive_count>\+\d+) |       #   +<number>
     (?P<positive_symbols>\++) |       #   +++...  This includes the single '+'
     (?P<negative_count>-\d+)  |       #   -<number>
     (?P<negative_symbols>-+)          #   ---...  including a single '-'
    )?                                 # and are optional
    (?P<error_1>[^\]]+)?               # If there's anything left, it's an error
  ) | (                                # End of parsing stuff in []s, except
    (?P<error_2>[^\]]*)                # If there was an error, we get here
  ))
  ((?P<close_bracket>\])|              # End bracket
   (?P<error_3>$))                     # unexpectedly reached end of string
)
""", re.X)

bond_fields = ["bond"]
bond = re.compile(r"(?P<bond>[=#/\\:~-])")

dot_fields = ["dot"]
dot = re.compile(r"(?P<dot>\.)")

closure_fields = ["closure"]
closure = re.compile(r"(?P<closure>\d|%\d\d?)")

close_branch_fields = ["close_branch"]
close_branch = re.compile(r"(?P<close_branch>\))")

open_branch_fields = ["open_branch"]
open_branch = re.compile(r"(?P<open_branch>\()")

# Make from the state name to
#   1. the regular expession to try to match and
#   2. the list of field names in that regexp.
# (There's no way to use sre to get #2 as an ordered list so
#  it needs to be done manually.)
info = {
    "atom": (atom, atom_fields),
    "bond": (bond, bond_fields),
    "dot": (dot, dot_fields),
    "closure": (closure, closure_fields),
    "close_branch": (close_branch, close_branch_fields),
    "open_branch": (open_branch, open_branch_fields),
    }

# Mapping from current state to allowed states
table = {
    # Could allow a dot
    "start": ("atom",),

    # CC, C=C, C(C)C, C(C)C, C.C, C1CCC1
    "atom": ("atom", "bond", "close_branch", "open_branch", "dot", "closure"),

    # C=C, C=1CCC=1
    "bond": ("atom", "closure"),

    # C(C)C, C(C)=C, C(C).C, C(C(C))C, C(C)(C)C
    "close_branch": ("atom", "bond", "dot", "close_branch", "open_branch"),

    # C(C), C(=C), C(.C) (really!)
    "open_branch": ("atom", "bond", "dot"),

    # C.C  -- allow a dot? as in C..C
    "dot": ("atom",),

    # C1CCC1, C1=CCC1, C12CC1C2, C1C(CC)1, C1(CC)CC1, c1ccccc1.[NH4+]
    "closure": ("atom", "bond", "closure", "close_branch", "open_branch", "dot"),
}

# Parse a SMILES string and print the events found
def tokenize(s, handler = Handler.TokenHandler()):
    expected = table["start"]
    n = len(s)
    i = 0
    handler.begin()
    while i < n:
        # Of the expected states, find one that matches the
        # text at the current position.
        for state in expected:
            pat, fields = info[state]
            m = pat.match(s, i)
            if m:
                break
        else:
            # No matches found, so this was an error
            handler.error("Unknown character", i, s[i:])
            # The handler is allowed to not throw an
            # exception, but we are done, so return.
            return
        #print "New state:", state,

        # Get the dictionary of matched name groups
        d = m.groupdict()

        # Go through the list of fields that could have matched.
        # Needs to go in a order so the token text can be converted
        # back into the original string.
        for field in fields:
            # See if there was a match for the given named field
            if d[field] is not None:
                # Was it an error match?
                if field[:5] == "error":
                    pos = m.start(field)
                    if field == "error_3":
                        handler.error("Missing ']'", pos, s[pos:])
                    else:
                        handler.error("Unknown character", pos, s[pos:])
                    return
                # Success, so send the token to the callback
                #print "--> ", m.group(field), field
                handler.add_token(field, i, m.group(field))

        # Get the new list of expected states, and move to
        # the end of the previous match.
        expected = table[state]
        i = m.end(0)

    handler.end()
