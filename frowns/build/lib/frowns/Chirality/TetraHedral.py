"""TetraHedral.T(order, chirality)
order = atoms in the input order of the smiles string.
chirality = @ anti-clockwise @@ clockwise

TetraHedral.T.getChirality(neworder) -> return the chirality
of the new order of the atoms

order -> atom order of the molecules
         around the center.
         If the number of atoms is 3 then the
         a hydrogen atom is assumed to complete
         the chirality.

Lookup table
chirality = '@'
input order = (1,2,3,4)
neworder, chirality
1,2,3,4 : '@'
"""
REVERSE = 0
SAME = 1

chiral_table = {
    (0,1,2,3): SAME,
    (0,2,3,1): SAME,
    (0,3,2,1): SAME,
    (1,0,2,3): REVERSE,
    (1,2,3,0): REVERSE,
    (1,2,3,0): REVERSE,
    (2,0,1,3): SAME,
    (2,1,3,0): SAME,
    (2,3,0,1): SAME,
    (3,0,1,2): REVERSE,
    (3,1,2,0): REVERSE,
    (3,2,0,1): REVERSE,
    (0,1,3,2): REVERSE,
    (0,3,2,1): REVERSE,
    (0,2,1,3): REVERSE,
    (1,0,3,2): SAME,
    (1,3,2,0): SAME,
    (1,2,0,3): SAME,
    (2,0,3,1): REVERSE,
    (2,3,1,0): REVERSE,
    (2,1,0,3): REVERSE,
    (3,0,2,1): SAME,
    (3,2,1,0): SAME,
    (3,1,0,2): SAME,
    (0,1,2): SAME,
    (1,2,0): SAME,
    (2,0,1): SAME,
    (0,2,1): REVERSE,
    (2,1,0): REVERSE,
    (1,0,2): REVERSE
    }

class T:
    def __init__(self, order, chirality):
        # store the initial order
        self.order = order
        self._initialOrder = [x.handle for x in order]
        self.chirality = "@"
        # normalize to anti-clockwise ordering
        if chirality == "@@":
            order[-1], order[-0] = order[-0], order[-1]

    def __str__(self):
        text = "Chirality\n"\
               " %s -> %s"%(self.order, self.chirality)
        return text
        
    def getChirality(self, order):
        """(order)->what is the chirality of a given order of
        atoms?"""
        indices = tuple([self._initialOrder.index(atom.handle)
                         for atom in order])
        same = chiral_table[indices]
        if same:
            return self.chirality
        else:
            if self.chirality == "@": return "@@"
            else: return "@"
                
if __name__ == "__main__":
    class Atom:
        def __init__(self, id):
            self.handle = id
        def __repr__(self):
            return "Atom(%s)"%self.handle
        
    a = Atom(0)
    b = Atom(1)
    c = Atom(2)
    d = Atom(3)

    c0 = T([a,b,c,d], "@")

    print c0.getChirality([a,b,c,d])
    print c0.getChirality([a,b,d,c])
    
    
    
