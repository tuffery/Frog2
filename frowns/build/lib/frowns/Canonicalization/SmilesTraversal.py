"""Traversals - a traversal (or path) through a graph

Traversals are stored as a list of tokens representing each object in a graph.
Each token object has __str__ defined so once a list has been formed
it can be turned into a string object by

 "".join(traversal_list)

Tokens are defined by the representation desired.  See SmilesTokens.py for
a list of SmilesTokens.

One thing to remember is that closures are considered atom properties for
the purposes of traversals.  That is the AtomToken is reponsible for storing
information about closures.  This information is built up during the
depth first search of the graph.
"""

import SmilesTokens

class IDGenerator:
    def __init__(self):
        self.index = 0

    def next(self):
        self.index += 1
        return self.index
    
class SmilesTraversal:
    # override the following to create different
    # tokens for traversals
    AtomToken = SmilesTokens.Atom
    BondToken = SmilesTokens.Bond
    BranchToken = SmilesTokens.Branch
    BranchEndToken = SmilesTokens.BranchEnd

    def __init__(self, parent=None):
        if parent is None:
            atomsDone = {}
            idGenerator = IDGenerator()
            closureIdGenerator = IDGenerator()
            closures = {}
        else:
            atomsDone = parent.atomsDone
            idGenerator = parent.idGenerator
            closureIdGenerator = parent.closureIdGenerator
            closures = parent.closures

        self.atoms = []
        self.bonds = []
        self.data = []
        self.atomsDone = atomsDone
        self.idGenerator = idGenerator
        self.closureIdGenerator = closureIdGenerator
        self.closures = closures

    def addAtom(self, atom):
        #print "atom", atom,
        atomToken = self.AtomToken(atom,
                                   self.closures,
                                   self.closureIdGenerator)
        self.atomsDone[atom] = atomToken
        self.data.append(atomToken)
        self.atoms.append(atom)
        #print self.atomsDone

    def addBond(self, bond):
        #print "bond", bond
        self.data.append(self.BondToken(bond))
        # XXX HACK -> give the bond a traversal order
        self.bonds.append(bond)

    def addClosure(self, atom1, atom2, bond):
        #print "closure", atom1, atom2, bond
        atomsDone = self.atomsDone
        closures1 = atomsDone[atom1].closures
        closures2 = atomsDone[atom2].closures
        assert closures1 is not closures2
        id = self.idGenerator.next()
        
        bondToken = self.BondToken(bond)
        # only pass the bond info to the first atom
        # the other gets a None (or don't print bond info)
        closures1.append((id, bondToken))
        closures2.append((id, None))
        self.bonds.append(bond)

    def addBranch(self):
        self.data.append(self.BranchToken())

    def addBranchEnd(self):
        self.data.append(self.BranchEndToken())

    def append(self, traverse):
        """(traverse)->append the traverse to the current traverse"""
        self.data.extend(traverse.data)
        self.atoms.extend(traverse.atoms)
        self.bonds.extend(traverse.bonds)

    def __str__(self):
        # XXX HACK -> find the bond traversal order
        # and make a note of it.  This is pretty ugly
        index = 1
        for x in self.data:
            if x.__class__ == SmilesTokens.Bond:
                x.bond._traverseOrder = index
                index += 1
        return "".join(map(str,self.data))


class SmartsTraversal(SmilesTraversal):
    AtomToken = SmilesTokens.SmartsAtom
    BondToken = SmilesTokens.SmartsBond

class IsomericSmilesTraversal(SmilesTraversal):
    AtomToken = SmilesTokens.IsomericAtom
