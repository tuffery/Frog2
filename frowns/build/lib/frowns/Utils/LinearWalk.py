def _dfswalk(atom, visitedAtoms, visitedBonds, nonPathBonds, path, paths, depth, mindepth, maxdepth):
    if depth >= maxdepth:
        return

    mark = len(path)
    markNonPathBonds = len(nonPathBonds)
    
    for bond in atom.bonds:
        oatom = bond.xatom(atom)

            
        if not visitedAtoms.has_key(oatom):
            visitedAtoms[oatom] = 1
            visitedBonds[bond] = 1
            
            path.append(bond)
            path.append(oatom)

            # keep a record of all the bonds that are not in the linear path
            # but branch back to it
            for bond2 in oatom.bonds:
                if bond2 not in visitedBonds and bond2 not in nonPathBonds:
                    a1, a2 = bond2.atoms
                    if a1 is oatom:
                        check = a2
                    else:
                        check = a1
                    if check in visitedAtoms:
                        print bond2, nonPathBonds
                        assert bond2 not in nonPathBonds
                        nonPathBonds.append(bond2)
            
            if (len(path)+1)/2 >= mindepth:
                paths[tuple(path), tuple(nonPathBonds)] = 1


            _dfswalk(oatom, visitedAtoms, visitedBonds, nonPathBonds, path, paths, depth+1, mindepth, maxdepth)

            while len(path) != mark:
                path.pop()

            while len(nonPathBonds) != markNonPathBonds:
                nonPathBonds.pop()

            del visitedAtoms[oatom]
            del visitedBonds[bond]

def generatePaths(molecule, mindepth=6,  maxdepth=7):
    """(molecule, mindepth maxdepth) -> linear paths
    Generate all linear paths through a molecule from mindepth atoms
    to maxdepth atoms
    """
    paths = {}
    for atom in molecule.atoms:
        _dfswalk(atom, {atom:1}, {}, [], [atom], paths, 1, mindepth, maxdepth)

    return paths.keys()


if __name__ == "__main__":
    from frowns import Smiles
    m = Smiles.smilin("CCCNc1ccccc1")
    for paths in generatePaths(m):
        print paths
