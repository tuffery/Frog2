"""perception

Potentially perceptable people perceive partial properties.

A lot of computation chemistry involves the perception of various
chemical properties.  The term perception is used because some
chemists have different beliefs about what defines chemical
properties.

The list of properties that frowns perceives is the following:

Smallest Set of Smallest Rings
Aromaticity (which depends on the smallest set of smallest rings)

These routines are designed to be flexible, that is they can be
relaplaced by other routines easily.

Each of the derived perception routines has the following
signature:

molecule = perception(molecule)

This format is retained because in the future a copy of the molecule
with the appropriate annotations might be returned.

Currently for the sake of speed the molecule is modified (not
a copy of the molecule)  This may change in the future.
"""

