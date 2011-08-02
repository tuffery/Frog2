def isChiral(atom):
    """(atom) -> See if the bfs atoms can form a chiral center"""
    q = {}
    weights = {}
    for oatom in atom.oatoms:
        q[oatom] = [{atom:1, oatom:1}, oatom.number, [oatom]]
        weights[oatom] = oatom.number + oatom.hcount


    if len(q) == 3:
        if atom.hcount != 1:
            return 0

        weights['H'] = 1

    level = 1
    while q:
        molecularWeights = {}
        for start, boundary in q.items():
            next_layer = []

            visited, mw, queue = boundary

            oldmw = mw
            for oatom in queue:
                for next in oatom.oatoms:
                    if next not in visited:
                        mw += level * (next.number)
                        visited[next] = 1
                        next_layer.append(next)

            if mw != oldmw:
                weights[start] = level, mw            

            if not next_layer:
                del q[start]
            else:
                q[start] = visited, mw, next_layer

        level += 1

        d = {}
        for start, value in weights.items():
            d[value] = 1
        if len(d) == 4:
            return 1

    return 0

def getStereoCenters(molecule):
    """(molecule) -> return a list of the atom stereo centers"""
    atomsToCheck = []
    
    for atom in molecule.atoms:
        if atom.symbol == "C" and \
            ( (len(atom.bonds) == 4 and atom.hcount == 0) or \
              (len(atom.bonds) == 3 and atom.hcount == 1) ):
            atomsToCheck.append(atom)

    chiralAtoms = []

    #print atomsToCheck
    for atom in atomsToCheck:
        types = {}
	#print atom.oatoms
        for attachedAtom in atom.oatoms:
            symbol = attachedAtom.symbol
            
            if symbol in types:
                if isChiral(atom):
                    chiralAtoms.append(atom)
                break
            
            types[symbol] = 1
            
        else:
            chiralAtoms.append(atom)
	    print "else"
            
    return chiralAtoms
        

##if __name__ == "__main__":
##    from frowns import Smiles
##    from frowns import MDL
##    from frowns.Depict import wxMoleculeDrawer
    
##    #mdl = MDL.mdlin(open("..\\..\\DrugProto1\\Widgets\\1.sdf"))
##    mdl = MDL.mdlin(open("test_6.sdf"))
##    mol = mdl.next()

##    i = 0

##    atoms = getStereoCenters(mol)

##    colors = {}
##    for atom in atoms:
##        for bond in atom.bonds:
##            colors[bond] = "red"
            
##    _app = wxMoleculeDrawer.testApp(0)
##    _app.draw(mol, colors)
##    _app.MainLoop()

####    for smiles in ["C(N)(O)(S)",
####                   "C(N)(O)(S)(C)",
####                   "C(N)(C)(O)[Fe]"]:
                   
####        m = Smiles.smilin(smiles)
####        print getStereoCenters(m)    
