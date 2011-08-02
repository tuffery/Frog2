from frowns import Smiles, Smarts, MDL, Fingerprint
from frowns.perception import RingDetection


def removeTerminalAtoms(m):
    while 1:
        to_remove = []
        for atom in m.atoms:
            if atom.symbol == "C" and len(atom.bonds) == 1:
                to_remove.append(atom)

        if not to_remove:
            break

        for atom in to_remove:
            m.remove_atom(atom)

    return m
            
m = Smiles.smilin("c1ccccc1CCNc1cccc1CCC")

print m.cansmiles()
removeTerminalAtoms(m)
smarts = m.arbsmarts()
print smarts

pat = Smarts.compile(smarts)
assert pat.match(m)

smiles = []

file = open("F:\\Chemical Libraries\\Libraries to Purchase\\Nat_425.sdf")

reader = MDL.mdlin(file)
mol = reader.next()

print "reading natural products"
i = 0
while mol:
    fp = Fingerprint.generateFingerprint(mol)
    smile = mol.cansmiles()
    mol = removeTerminalAtoms(mol)
    smiles.append((len(mol.atoms), len(mol.bonds), smile, mol.arbsmarts(), i, fp))
    i += 1
    mol = reader.next()

search = smiles
search.sort()

smiles = []


print "processing"

dict = {}
indicesUsed = {}

# trim the search list
final_search_list = []

while search:
    atoms, bonds, smile, smarts, index, fp = search.pop(0)
    final_search_list.append((atoms, bonds, smile, smarts, index, fp))
    
    pat = Smarts.compile(smarts)
    next = []
    key = (smarts, index)
    for atoms, bonds, smile, smarts2, index, fp2 in search:
        if smarts == smarts2:
            dict[key] = dict.get(key, []) + [(smile, index)]
        elif fp in fp2 and pat.match(Smiles.smilin(smile), 1):
            dict[key] = dict.get(key, []) + [(smile, index)]
        else:
            next.append((atoms, bonds, smile, smarts2, index, fp2))

    search = next

print "reduced search list to", len(final_search_list)

file = open("F:\\Chemical Libraries\\Libraries to Purchase\\Semi-Nat_9575.sdf")
reader = MDL.mdlin(file)
mol = reader.next()

print "reading semi-natural products"
while mol:
    fp = Fingerprint.generateFingerprint(mol)
    smile = mol.cansmiles()
    mol = removeTerminalAtoms(mol)
    smiles.append((len(mol.atoms), len(mol.bonds), smile, mol.arbsmarts(), i, fp))
    i += 1
    mol = reader.next()
    
smiles.sort()

while final_search_list:
    atoms, bonds, smile, smarts, index, fp = final_search_list.pop(0)
    
    pat = Smarts.compile(smarts)
    next = []
    key = (smarts, index)
    for atoms, bonds, smile, smarts2, index, fp2 in smiles:
        if smarts == smarts2:
            dict[key] = dict.get(key, []) + [(smile, index)]
        elif fp in fp2 and pat.match(Smiles.smilin(smile), 1):
            dict[key] = dict.get(key, []) + [(smile, index)]
        else:
            next.append((atoms, bonds, smile, smarts2, index, fp2))

    smiles = next  
