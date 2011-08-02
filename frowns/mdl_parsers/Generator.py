"""Generator

 Generator is a class that is used for reading linearized graph notations
 examples of linearized graph notations are smiles (http://www.daylight.com)
 although Tripos has one as well (http://www.tripos.com)

 A linearized notation consists of these parts:

  1 atoms
  2 bonds
  3 branches - specified by ( and ) in smiles
  4 closures - specifed by numbers or %12 in smiles

 A generator conforms to the following interface

 class Generator:
   def __init__(self):
      ...
   def addEvent(self, event, data):
       '''(event, data) -> add event with event data 'data' to
       the generator.

       the event and forms of the data are given from the
       description of the parsing routine'''

   def mol(self):
       '''-> return the generated molecule'''
"""
from frowns.Atom import Atom
from frowns.Bond import Bond
from frowns.Molecule import Molecule
from frowns.Objects.ObjectGenerator import parse_atom, parse_bond

def clearDict(dict):
    for k,v in dict.items():
        if v == None:
            del dict[k]
    return dict

# this simply allows a dictionary to be
#  hased
class IndexGenerator:
    def __init__(self, start=0):
        self.index = start -1

    def next(self):
        self.index += 1
        return self.index
    
class DictWrapper:
    def __init__(self, data, idGenerator=IndexGenerator()):
        global INDEX
        self.data = data
        self.index = idGenerator.next()

    def __cmp__(self, other):
        return cmp(self.data, other.data)

    def __str__(self):
        return "%s.%s"%(self.index, str(self.data))

    def __repr__(self):
        return "%s.%s"%(self.index, str(self.data))
    
    def __hash__(self):
        return self.index

#
# The generator has callbacks for the following events
#  'atom' -> adds an atom
#  'bonds' -> adds a bond
#  'branch_out' -> pushes a branch
#  'branch_in' -> ends the last branch
#  'closure' -> adds a closure or completes a closure
#  'dot' -> a special case of a blank (empty) bond
#
# Implementation Notes
#  lastAtom = None
#  lastBond = None
#  closures = {}
#  when an atom is added
#    if there is no lastAtom then do nothing.
#    else create a blank lastBond if none exists
#         use lastBond to connect atom to lastAtom
#    update lastAtom to be the new atom
#    clear lastBond because it has been used
#
#  when a bond is added
#    make sure that lastBond is None otherwise
#    there are two bonds in a row
#
#  branch_out
#    there is a branch stack which stores the current
#    lastAtom state
#
#  branch_in
#    go back to the last branch stack state thus
#    closing the branch
#
#  closure
#    a closure maps an atom to a closure index
#    when the closure index is seen again the current
#    atom is bonded to the closure atom.
#
#    XXX FIX ME
#     Daylight adds the bond to the molecule when the
#     closure index first appears.
#
#     This code adds the bond when the closure is completed.
#     So bonds will not be in the same order as daylight.
class Generator:
    def __init__(self, Atom=Atom, Bond=Bond, Molecule=Molecule):
        # This dispatches the parser events
        self.dispatch = {'atom': self.addAtom,
                         'bond': self.addBond,
                         'branch_out': self.pushBranch,
                         'branch_in': self.popBranch,
                         'closure': self.addClosure,
                         'dot': self.clearLastBond
                         }
        self.atoms = []
        self.bonds = []
        
        self.bondedAtoms = {}
        self.branches = []
        self.closures = {}
	self.lastBond = None
        self.lastAtom = None
        self.Atom = Atom
        self.Bond = Bond
        self.Molecule = Molecule

    def addState(self, state, data):
        self.dispatch[state](data)

    def addAtom(self, atom):
        """(atom)->add atom to the atom stack
        combine with lastBond if one exists and set lastBond
        to None since it has been used"""
        atom = parse_atom(atom, self.Atom)
        self.atoms.append(atom)
        
        #print "adding atom", atom
        if self.lastAtom and not self.lastBond:
            self.addBond({"bond":""})

        if self.lastBond:
            self.bondedAtoms[self.lastBond] = (self.lastAtom, atom, 0)
            self.bonds.append(self.lastBond)
            #print '\t', self.lastBond, "->", (self.lastAtom, atom)
            self.lastBond = None
        self.lastAtom = atom
        
    def addBond(self, bond):
        bond = DictWrapper(bond)
        #print "adding bond", bond
        self.lastBond = bond

    def addClosure(self, closure):
        """(closure) -> add a closure"""	
        closureIndex = closure['closure']
        closures = self.closures

        # if we don't have the closure index
        #  then create one
        if not closures.has_key(closureIndex):
            # need to create an unspecified bond if
            # it doesn't exist
            if self.lastBond is None:
                self.lastBond = DictWrapper({'bond':''})
                
            closures[closureIndex] = self.lastAtom, self.lastBond
        else:
            closureAtom, closureBond = closures[closureIndex]
            # consume the closure (we may see the same used
            #                      in different closures more than once)
            del closures[closureIndex]
            
            # we need to determine which of the closure bonds
            # specify the bond type
            if self.lastBond is None:   # lastBond is not specified
                bondToUse = closureBond
            elif not closureBond.data['bond']:       # closureBond is not specified
                bondToUse = self.lastBond
            elif not self.lastBond.data['bond']:     # lastBond is not specified
                bondToUse = closureBond
            else:
                # we need to make sure that both bonds are the same
                bondToUse = closureBond
                if closureBond != self.lastBond:
                    raise "Closure specified with two different bond types %s!=%s"%(closureBond, self.lastBond)

            self.bondedAtoms[bondToUse] = (self.lastAtom, closureAtom, 1)
            self.bonds.append(bondToUse)
                
        self.lastBond = None        
                               
    def clearLastBond(self, data):
        if self.lastBond is not None:
            raise "Syntax Error . follows bond"
        self.lastAtom = None
    
    def pushBranch(self, data):
        self.branches.append(self.lastAtom)

    def popBranch(self, data):
        if not self.branches:
            raise "Unbalanced branch"
        if self.lastBond:
            raise "Illegal bond before branch end"
        self.lastAtom = self.branches.pop()

    def mol(self):
        if self.closures:
            raise "Closure not completed %s"%(self.closures.keys()[0])
        mol = self.Molecule()
        index = 0
        for atom in self.atoms:
            atom.symorder = index
            mol.add_atom(atom)
            index += 1

        bondedAtoms = self.bondedAtoms
        Bond = self.Bond
        for bond_data in self.bonds:
            bond = parse_bond(bond_data.data, Bond)
            a1, a2, closure = bondedAtoms[bond_data]
            if closure:
                bond._closure = 1
            assert a1.parent and a2.parent
            mol.add_bond(bond, a1, a2)
        
        return mol
    
