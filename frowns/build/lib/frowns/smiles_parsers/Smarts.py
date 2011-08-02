#!/usr/bin/enb python

import string, re
import Handler

#######################

# Define some regular expressions inside a quoted string
# then turn the string into the actual data structure.
# (I found it was easiest to understand when done this way.)

definitions = r"""
# These are the atomic symbols Daylight allows outside of []s
# See "atom_class" for names like "a" and "A"
raw_atom   Cl|Br|[cnospBCNOFPSI*]

# For atoms inside of []s
open_bracket    \[
close_bracket   \]

# See "element_modifiers" for the patterns for element names
# charges, chiralities, H count, etc.

# [235U]
weight  \d+

# [#6]
atomic_number   #\d+

# [!C]
atom_not  !

# & is highest (an "and")
# , is next (an "or")
# ; is lowest (an "and")
#  [n&H]   [n,H]  [c,h;H1]
atom_binary   [&,;]

# C.C
dot  \.

# -   single bond (aliphatic)
# /   directional single bond "up"
# \   directional single bond "down"
# /?  directional bond "up or unspecified"
# \?  directional bond "down or unspecified"
# =   double bond
# #   triple bond
# :   aromatic bond
# ~   any bond (wildcard)
# @   any ring bond
bond  [/\\]\??|[=#:~@-]

# *!:*  -- not aromatic
bond_not  !

# *@;!:*  -- same as !:
bond_binary   [&;,]

# (C).(C)
open_zero   \(
# C(C)
open_branch  \(

# [$(*C);$(*CC)]
open_recursive_smarts  \$\(

# special cased because it closes open_zero, open_branch, and
# recursive_smarts
close_parens  \)

# Ring closures, 1, %5 %99 (and even %00 for what it's worth)
closure    \d|%\d\d?
"""
#######################

# Turn the above string into key/value pairs where the
# values are the compiled regular expressions.
info = {}
for line in string.split(definitions, "\n"):
    line = string.strip(line)
    if not line or line[:1] == "#":
        continue
    name, pattern = string.split(line)
    info[name] = re.compile(pattern)
del line, name, pattern

info["atom_class"] = re.compile(r"""
(?P<raw_aromatic>a)|    # Not really sure what these mean
(?P<raw_b_unknown>b)|
(?P<raw_f_unknown>f)|
(?P<raw_h_unknown>h)|
(?P<raw_i_unknown>i)|
(?P<raw_r_unknown>r)|
(?P<raw_aliphatic>A)|
(?P<raw_R_unknown>R)
""", re.X)

# 'H' is used for the hydrogen count, so those searches require a
# special recursive SMARTS definition.  Eg, for deuterium or tritium
#  [$([2H]),$([3H])]
# This is implemented as a special-case hack.  Note: if there's
# an error in the parse string in this section then the error
# location will point to the start of this term, not at the
# character that really caused the error.  Can be fixed with an
# 'error_'  like I did for the SMILES -- not needed for now.  XXX

hydrogen_term_fields = [
    "open_recursive_smarts",
    "open_bracket",
    "weight",
    "element",
    "positive_count",
    "positive_symbols",
    "negative_count",
    "negative_symbols",
    "close_bracket",
    "close_recursive_smarts",
    ]

info["hydrogen_term"] = re.compile(r"""
(?P<open_recursive_smarts>\$\()
(?P<open_bracket>\[)
(?P<weight>\d+)?                # optional molecular weight [2H]
(?P<element>H)                  # Must be a hydrogen
(                               # optional charge
 (?P<positive_count>\+\d+)|     # +3
 (?P<positive_symbols>\++)|     # ++
 (?P<negative_count>\-\d+)|     # -2
 (?P<negative_symbols>\-+)|     # ---
)?
(?P<close_bracket>\])
(?P<close_recursive_smarts>\))
""", re.X)

element_symbols_pattern = \
   r"C[laroudsemf]?|Os?|N[eaibdpos]?|S[icernbmg]?|P[drmtboau]?|" \
   r"H[eofgas]|c|n|o|s|p|A[lrsgutcm]|B[eraik]?|Dy|E[urs]|F[erm]?|" \
   r"G[aed]|I[nr]?|Kr?|L[iaur]|M[gnodt]|R[buhenaf]|T[icebmalh]|" \
   r"U|V|W|Xe|Yb?|Z[nr]|\*"

info["element_modifier"] = re.compile(r"""
 (?P<element>
      # This does *not* contain H.  Hydrogen searches must be done
      # with a special recursive SMARTS.  On the other hand, it does
      # include the lower case aromatic names.
""" + element_symbols_pattern + r"""
 )|
 (?P<aromatic>a)|            # aromatic
 (?P<aliphatic>A)|           # Aliphatic
 (?P<degree>D\d+)|           # Degree<n>
 (?P<total_hcount>H\d*)|     # total Hydrogen count<n> (defaults to 1)
 (?P<imp_hcount>h\d*)|       # implicit hydrogen count<n> (defaults to 1)
 (?P<ring_membership>R\d*)|  # in <n> Rings (no n means any rings)
 (?P<ring_size>r\d*)|        # in a ring of size <n> (no n means any rings)
 (?P<valence>v\d+)|          # total bond order of <n>
 (?P<connectivity>X\d+)|     # <n> total connections

 (?P<positive_count>\+\d+)|  # +2  +3
 (?P<positive_symbols>\++)|  # +   ++   +++
 (?P<negative_count>\-\d+)|  # -1  -4
 (?P<negative_symbols>\-+)|  # --  -    -------

  # XXX What about chiral_count?
 (?P<chiral_named>           # The optional '?' means "or unspecified"
   @TH[12]\??|               # @TH1  @TH2?
   @AL[12]\??|               # @AL2?
   @SP[123]\??|              # @SP3  @SP1?
   @TB(1[0-9]?|20?|[3-9])\??|           # @TH{1 through 20}
   @OH(1[0-9]?|2[0-9]?|30?|[4-9])\??    # @OH{1 through 30}
 )|
 (?P<chiral_symbols>@@?\??)  # @ (anticlockwise) or @@ (clockwise)
""", re.X)

# The ')' closes three different open parens.  This maps from the
# previous open state to the appropriate close state.
close_parens_states = {
    "open_branch": "close_branch",
    "open_recursive_smarts": "close_recursive_smarts",
    "open_zero": "close_zero",
}

#### Some helpful definitions to reduce clutter and complication

# Possible transitions from the start node.  Also visited after
# a '.' disconnect or in a recursive SMARTS.
expecting_start = ("raw_atom", "atom_class", "open_bracket", "open_zero")

# Looking for node definition, like "C" or "a" or "["
expecting_atom = ("raw_atom", "atom_class", "open_bracket")

# Inside of []s: 235U, #6, R, $([2H]), $(*=C), !
expecting_element_start = ("weight", "atomic_number",
                           "element_modifier", "hydrogen_term",
                           "open_recursive_smarts", "atom_not")

# the ';' in [n;H1] or the ']' at the end
expecting_element_end = ("atom_binary", "close_bracket")

# All bonds start with a '!' or one of the bond symbols
expecting_bond_start = ("bond", "bond_not")

expecting_raw_term = expecting_atom + expecting_bond_start + \
                     ("close_parens", "open_branch", "dot", "closure")
expecting_modifier = ("element_modifier", "open_recursive_smarts")

table = {
    "start": expecting_start,

    # (C).(R).[U].([$(*)])
    "open_zero": ("raw_atom", "atom_class", "open_bracket"),
    # as well as (CC(C))
    "close_zero": ("dot", "close_parens"),

    # A raw term are the things like 'C', '[U]', '%10', '.', '(', '!#'
    "raw_atom": expecting_raw_term,
    # An atom_class is a non-specific atom term, like 'A' or 'r'
    "atom_class": expecting_raw_term,

    # the []s
    "open_bracket": expecting_element_start,
    "close_bracket": expecting_raw_term,

    # Yes, '[!!!!C]' is legal, according to the docs, but it isn't
    # supported by the parser, unless you optimze it.
    "atom_not": expecting_element_start,
    "atom_binary": expecting_element_start,

    # "14N", "14a", ...
    # Note that weight can only be set once so it isn't a modifier
    # Also, "14#6" isn't legal (tested against the toolkit)
    "weight": expecting_modifier,

    # "#6R2" or "#8," or "#7]"
    # The atomic_number can only be set once so it isn't a modifier
    "atomic_number": expecting_modifier + expecting_element_end,

    # All of these are type of modifiers
    "element_modifier": expecting_modifier + expecting_element_end,
    "hydrogen_term": expecting_modifier + expecting_element_end,
    "close_recursive_smarts": expecting_modifier + expecting_element_end,

    # This it the recursive part -- goes back to the beginning
    "open_recursive_smarts": expecting_start,

    # C=C, C1CCC=1, C~-C, C=(C)C, C=,-C
    "bond": expecting_atom + ("closure", "bond", "open_branch",
                              "bond_binary"),

    # C!!=C
    "bond_not": expecting_bond_start,
    # C=,-C
    "bond_binary": expecting_bond_start,

    "closure": expecting_raw_term,

    "close_branch": expecting_raw_term,
    "open_branch": expecting_atom + expecting_bond_start + ("dot",),

    # After a "." we can start all over again
    "dot": expecting_start,
}

def tokenize(s, handler = Handler.TokenHandler()):
    expected = table["start"]
    parens_stack = []
    n = len(s)
    i = 0
    handler.begin()
    while i < n:
        for state in expected:
            m = info[state].match(s, i)
            if m:
                break
        else:
            handler.error("Unknown character", i, s[i:])
            return
        if close_parens_states.has_key(state):
            parens_stack.append(state)
        elif state == "close_parens":
            try:
                state = close_parens_states[parens_stack.pop()]
            except IndexError:
                # Too many close parens
                handler.error("Too many ')'", i, s[i:])
                return

        d = m.groupdict()
        if d and state == "hydrogen_term":
            # Special case the hydrogen term
            for field in hydrogen_term_fields:
                if d[field] is not None:
                    handler.add_token(field, i, d[field])
            #print " --> New state:", state
        else:
            name = state
            if d:
                # There should only be one match
                for name, v in d.items():
                    if v is not None:
                        break

            handler.add_token(name, i, m.group(0))
            
        expected = table[state]
        i = m.end(0)

    handler.end()
