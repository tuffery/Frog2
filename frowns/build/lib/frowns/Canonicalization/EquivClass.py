"""computes equivalence classes for atoms

equiv_class = compute_equiv_class(atom)
"""

def compute_equiv_class(atom):
    """(atom)->Computes a unique integer for an atom"""    
    try:
        equiv_class = atom.number + \
                      1000*(atom.charge+10) + \
                      100000*(atom.hcount) + \
                      1000000*(atom.weight)
    except TypeError:
        raise ValueError, \
              "Can't compute number from atom.number %s atom.charge %s atom.hcount %s"\
              " atom.weight %s"%(atom.number, atom.charge, atom.hcount, atom.weight)
    return equiv_class
