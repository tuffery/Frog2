"""Render a molecule using AT&T's GraphViz.

This is a poor man's rendering package and should be used
accordingly.  It is intended merely as a debugging tool
to see if molecules are formed somewhat correctly as
debugging smiles strings tends to make my brain wilt.

That being said, it's pretty cool :)  Be aware that while
GraphViz is opensourced it is not free for commercial
applications and any such use must be cleared with AT&T.

You may find GraphViz here
http://www.research.att.com/sw/tools/graphviz/
but please read their FAQ and license before downloading.

This package must be able to find the 'neato' executable
in the search path.
"""
import os, re
try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO

NEATO = "neato"
NEATO = "C:\\Progra~1\\ATT\\Graphviz\\bin\\neato.exe"

DOTPAT = re.compile(r'(atom[0-9]+) \[pos="([^,]+),([^"]+)')

class Error(Exception): pass

def addCoords(mol):
    """(mol) Assign x,y coordinates to the atoms of a molecule
    using AT&T's graph layout algorithm.  This technique is
    very fast but doesn't look very good so use with
    caution.

    The executable 'neato' must exist in the search path"""
    
    ofile = StringIO()
    atoms = mol.atoms
    bonds = mol.bonds
    print >> ofile, "graph TEST {"
    
    d = {}
    for atom in mol.atoms:
        d["atom%d"%atom.handle] = atom
        
    for bond in mol.bonds:
        atom1, atom2 = bond.atoms
        print >> ofile, "\tatom%s -- atom%s;"%(atom1.handle, atom2.handle)
        
    print >> ofile, "}"

    try:
        r, w = os.popen2(NEATO)
    except:
        raise Error("Unable to locate neato")

    r.write(ofile.getvalue())
    r.close()

    added = 0
    for line in w.readlines():
        line = line.strip()
        groups = DOTPAT.match(line)
        if groups:
            name, x,y = groups.groups()
            d[name].x = int(x)
            d[name].y = int(y)
            added += 1
    if added != len(atoms):
        raise Error("neato failed depiction")
