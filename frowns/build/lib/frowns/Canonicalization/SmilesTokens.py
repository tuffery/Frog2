"""SmilesTokens

This is a list of smiles tokens that define the smiles language.
There are for types of tokens

 Atoms - C or [1000C]
 Bonds - single, double, sacrifice fly
 closures - things like C1CCC1 (the ones are closures)
 branches (Branch Start, Branch End)
          - these are ( and ) respectively

 Each Token Class below can convert the appropriate object into
 a string through overloading the __str__ mechanism.
"""
#
# Canonical labeling rules for daylight atoms
# From the Daylight Theory Manual
# Organic subset atoms may be written without brackets if the number of
# attached hydrogens conforms to the lowest normal valence consistent
# with explicit bonds. "Lowest normal valences" are B (3), C (4),
# N (3,5), O (2), P (3,5), S (2,4,6), and 1 for the halogens.
# Atoms in aromatic rings are specified by lower case letters, e.g.,
# aliphatic carbon is represented by the capital letter C, aromatic
# carbon by lower case c.
ORGANIC_SUBSET = ['B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I']

class Atom:
    def __init__(self, atom, closureMappings, idGenerator):
        self.atom = atom
        self.closures = []
        self.closureMappings = closureMappings
        self.idGenerator = idGenerator

    def name_atom(self):
        atom = self.atom
        symbol = "%s"%(atom.symbol,)
        weight = atom.weight
        charge = atom.charge
        hcount = atom.hcount
        explicit_hcount = atom.explicit_hcount
        
        if atom.aromatic:
            out_symbol = symbol[0].lower() + symbol[1:]
        else:
            out_symbol = symbol

        # pyrole like nitrogens
        if atom.aromatic and symbol == "N" and charge == 0 and \
           weight == 0 and explicit_hcount == 1:
            # XXX Fix Me
            # There should only be one of these per five membered
            # aromatic ring
            return "[nH]"

        if symbol in ORGANIC_SUBSET and atom.valences and \
           not weight and not charge:
            sumOrders = atom.sumBondOrders()
            hcount = atom.hcount
            for valence in atom.valences:
                if sumOrders + hcount == valence:
                    return out_symbol#+"<%s>"%atom.handle

        if not weight: weight = ""
        
        if charge == -1: charge = "-"
        elif charge == 1: charge = "+"
        elif charge > 1: charge = "+%s"%charge
        else: charge = ""

        if hcount == 1: hcount = "H"
        elif hcount == 0: hcount = ""
        elif hcount > 1: hcount = "H%s"%hcount
        else:
            raise "Negative hcount!!!"

        return "[%s%s%s%s]"%(weight, out_symbol, hcount, charge)
                    
    def __str__(self):
        # easy for now
        atom = self.atom
        label = [self.name_atom()]
        
        idGenerator = self.idGenerator
        orderOfClosures = []
        closureMappings = self.closureMappings
        for id, bond in self.closures:
            if not closureMappings.has_key(id):
                closureMappings[id] = idGenerator.next()
            id = closureMappings[id]
            orderOfClosures.append((id,bond))
        orderOfClosures.sort()

        for id, bond in orderOfClosures:
            if bond: label.append(str(bond))
            if id > 9:
                label.append("%%%s"%id)
            else:
                label.append("%s"%id)
        return "".join(label)
#
# XXX FIX ME -> a lot of shared code here, it would be nicer
# if we can get rid of them
class IsomericAtom(Atom):
    def name_atom(self):
        atom = self.atom
        symbol = "%s"%(atom.symbol,)
        weight = atom.weight
        charge = atom.charge
        hcount = atom.hcount
        explicit_hcount = atom.explicit_hcount
        chirality = atom._chirality
        
        if atom.aromatic:
            out_symbol = symbol[0].lower() + symbol[1:]
        else:
            out_symbol = symbol

        # pyrole like nitrogens
        if atom.aromatic and symbol == "N" and charge == 0 and \
           weight == 0 and explicit_hcount == 1:
            # XXX Fix Me
            # There should only be one of these per five membered
            # aromatic ring
            return "[nH]"

        if symbol in ORGANIC_SUBSET and atom.valences and \
           not weight and not charge and not chirality:
            sumOrders = atom.sumBondOrders()
            hcount = atom.hcount
            for valence in atom.valences:
                if sumOrders + hcount == valence:
                    return out_symbol#+"<%s>"%atom.handle

        if not weight: weight = ""
        
        if charge == -1: charge = "-"
        elif charge == 1: charge = "+"
        elif charge > 1: charge = "+%s"%charge
        else: charge = ""

        if hcount == 1: hcount = "H"
        elif hcount == 0: hcount = ""
        elif hcount > 1: hcount = "H%s"%hcount
        else:
            raise "Negative hcount!!!"

        if chirality:
            bonds = [(bond._traverseOrder, bond) for bond in atom.bonds]
            bonds.sort()
            bonds = [bond[1].xatom(atom) for bond in bonds]
            chiralstr = chirality.getChirality(bonds)
        else:
            chiralstr = ""
        return "[%s%s%s%s%s]"%(weight, out_symbol, chiralstr, hcount, charge)
        return "[%s%s%s%s<%s>]"%(weight, out_symbol, hcount, charge, atom.handle)
    
class SmartsAtom(Atom):
    def name_atom(self):
        atom = self.atom
        symbol = atom.symbol

        if atom.aromatic:
            symbol = symbol.lower()
            
        if atom.symbol in ORGANIC_SUBSET:
            return symbol
        else:
            return "[%s]"%symbol
    

class Bond:
    def __init__(self, bond):
        self.bond = bond
        self.lookup = {
            1:"",
            2:"=",
            3:"#",
            4:"",
            5:"\\",
            6:"/"}
        
    def __str__(self):
        return self.lookup[self.bond.bondtype]    

class SmartsBond(Bond):
    def __init__(self, bond):
        self.bond = bond
        self.lookup = {
            1:"-",
            2:"=",
            3:"#",
            4:":",
            5:"\\",
            6:"/"}
        
class Branch:       
    def __str__(self):
        return "("

class BranchEnd:
    def __str__(self):
        return ")"
