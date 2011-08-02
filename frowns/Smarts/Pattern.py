import string
from frowns.IdGenerator import defaultGenerator

class Atom:
    def __init__(self, matcher, component_number,
                 generator = defaultGenerator):
        self.matcher = matcher
        self.component_number = component_number
        self.bonds = []
        self.handle = generator()
        
    def __eq__(self, atom):
        return self.matcher == atom

    def dump_info(self, bonds):
        text = "[%d] <%s>" % (self.component_number, str(self.matcher))
        text = text + " bonds = ["
        bond_ids = []

        for bond in self.bonds:
            index = 0
            for b in bonds:
                if bond is b: break
                index += 1
            bond_ids.append(str(index))
            
        text = text + string.join(bond_ids, ", ") + "]"
        return text
    
class Bond:
    def __init__(self, matcher, component_number,
                 generator = defaultGenerator):
        self.matcher = matcher
        self.component_number = component_number
        self.atoms = []
        self.handle = generator()
        
    def __eq__(self, bond):
        return self.matcher == bond

    def dump_info(self, atoms):
        index1 = index2 = index = 0
        for atom in atoms:
            if self.atoms[0] is atom: index1 = index
            if self.atoms[1] is atom: index2 = index
            index += 1
            
        return "[%d] %d <%s> %d" % (self.component_number,
                                    index1,
                                    self.matcher,
                                    index2)
    
class SmartsPattern:
    def __init__(self):
        self.atoms = []
        self.bonds = []

    def dump(self):
        lines = [self.__class__.__name__,
                 "atom_id [component_number] <matcher> "
                 "bonds = [list of bond ids]"]
        i = 0
        for atom in self.atoms:
            lines.append(str(i) + " " + atom.dump_info(self.bonds))
            i = i + 1
        lines.append("bond_id [component_number> "
                     "atom1_id <matcher> atom2_id")
        i = 0
        for bond in self.bonds:
            lines.append(str(i) + " " + bond.dump_info(self.atoms))
            i = i + 1
        return string.join(lines, "\n\t")

    def __del__(self):
        for atom in self.atoms:
            del atom.bonds
        for bond in self.bonds:
            del bond.atoms
