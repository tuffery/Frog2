class T:
    def __init__(self, order, chirality):
        order = self.order = [atom.handle for atom in order]
        self.chirality = "@"
        # normalize to anti-clockwise ordering
        if chirality == "@@":
            order[-2], order[-1] = order[-1], order[-2]

    def getChirality(self, order):
        """(order)->what is the chirality of a given order of
        atoms?"""
        indices = [self.order.index(atom.handle) for atom in order]

        # if less than four atoms are shown, then the first one
        # is assumed to be a hydrogen
        if len(self.order) == 4:            
            initial, rest = indices[0], indices[1:]
        else:
            initial, rest = None, indices

        minimum = min(rest)
        while rest[0] != minimum:
            first = rest.pop(0)
            rest.append(first)

        print rest

if __name__ == "__main__":
    from frowns.Atom import Atom
    a = Atom()
    b = Atom()
    c = Atom()
    d = Atom()
    center = Atom()

    T([a,b,c,d], "@")
    
    
