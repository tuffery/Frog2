def bfsCheckAtoms(atom):
    """(atom) -> See if the bfs atoms can form a chiral center"""
    q = {}
    for oatom in atom.oatoms:
        q[oatom] = [{atom:1, oatom:1}, oatom.mass, [oatom]]


    weights = {}
    if len(q) == 3:
        if oatom.hcount != 1:
            return 0

        weights['H'] = 1


    level = 1
    while q:
        molecularWeights = {}
        for start, boundary in q.items():
            next_layer = []

            visited, mw, queue = boundary
            for oatom in queue:
                for next in oatom.oatoms:
                    if next not in visited:
                        mw += level * next.mass + next.hcount
                        visited[next] = 1
                        next_layer.append(next)

            weights[start] = mw            

            if not next_layer:
                del q[start]
            else:
                q[start] = visited, mw, next_layer

        level += 1

        d = {}
        for value in weights.values():
            d[value] = 1
            
        if len(d) == 4:
            return 1

    return 0

def getStereoCenters(molecule):
    """(molecule) -> return a list of the atom stereo centers"""
    atomsToCheck = []
    
    for atom in molecule.atoms:
        if atom.symbol == "C" and (
            (len(atom.bonds) == 4 and atom.hcount == 0) or \
            (len(atom.bonds) == 3 and atom.hcount == 1) ):
            atomsToCheck.append(atom)

    chiralAtoms = []
    for atom in atomsToCheck:
        types = {}
        for attachedAtom in atom.oatoms:
            symbol = attachedAtom.symbol
            
            if symbol in types:
                if bfsCheckAtoms(atom):
                    chiralAtoms.append(atom)
                break
                
            types[symbol] = 1
            
        else:
            chiralAtoms.append(atom)

    return chiralAtoms
    

if __name__ == "__main__":
    from frowns import Smiles
    for smiles in ["C(N)(O)(S)",
                   "C(N)(O)(S)(C)",
                   "C(N)(C)(O)[Fe]"]:
                   
        m = Smiles.smilin(smiles)
        print getStereoCenters(m)    
