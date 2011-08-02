"""Generate the drawables for a molecule

generateDrawables(atoms, bonds) -> (points, lines, boundingbox)

points are in the form ((x,y,symbol)...)
lines are in the form ((bondType, (x1,y1), (x2,y2))...
bounding box is in the form (minx, maxx, miny, maxy)

"""
import math
from LineGenerator import make_double_line
from frowns.Canonicalization.SmilesTokens import Atom
ORGANIC_SUBSET = ['B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I']

def name_atom(atom):
    """->symbol
    Return the atom symbol for a depiction.  Carbons atoms
    in general are not returned"""

    symbol = "%s"%(atom.symbol,)
    weight = atom.weight
    charge = atom.charge
    hcount = atom.hcount
    explicit_hcount = atom.explicit_hcount

    out_symbol = symbol

    # pyrole like nitrogens
    if atom.aromatic and symbol == "N" and charge == 0 and \
       weight == 0 and explicit_hcount == 1 and hcount == 0:
        # XXX Fix Me
        # There should only be one of these per five membered
        # aromatic ring
        return "nH"


    if symbol in ORGANIC_SUBSET and atom.valences and \
       not weight and charge == 0:
        sumOrders = atom.sumBondOrders()
        hcount = atom.hcount
        sum = int(hcount + sumOrders)
        for valence in atom.valences:
            if valence == sum:
                if symbol == "C":
                    return ""
                else:
                    return symbol

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
    output = "%s%s%s%s"%(weight, out_symbol, hcount, charge)

    return output

def _recurse(atom, visited):
    """_recursion routine for detecting overlapping
    fragments"""
    
    for bond in atom.bonds:
        for next in bond.atoms:
            if not visited.has_key(next):
                visited[next] = 1
                _recurse(next, visited)    
                
def getOverlappingFragments(atoms):
    """given a list of atoms, return the overlapping
    components in a dictionary
    key = tuple of coordinates
    value = number of overlaps
    """
    _atoms = {}
    for atom in atoms:
        _atoms[atom] = 1

    components = {}

    while _atoms:
        start = _atoms.keys()[0]
        visited = {start:1}
        _recurse(start, visited)
        coords = []
        for atom in visited.keys():
            coords.append((atom.x, atom.y))
            del _atoms[atom]
        coords.sort()
        coords = tuple(coords)
        components[tuple(coords)] = components.get(coords,0)+1

    for key, val in components.items():
        if val == 1:
            del components[key]

    return components

def generateDrawables(atoms, bonds, drawAromatic=1, name_atom=name_atom):
    """(atoms, bonds, drawAromatic=1)->generate drawables from a collection of
    atoms and bonds.  The reason we aren't using a full molecule
    is that we would like to draw a portion of a molecule
    if necessary"""
    atomCoords = {}
    points = []
    lines = []
    wedges = []
    filledWedges = []
    minx = maxx = miny = maxy = None
    atomIndices = {}

    # find the min and max for the x and y coordinates
    for atom in atoms:
        x,y = atom.x, atom.y
        if minx is None or x<minx:
            minx = x
        if maxx is None or x>maxx:
            maxx = x
        if miny is None or y<miny:
            miny=y
        if maxy is None or y>maxy:
            maxy=y

        
    for atom in atoms:
        x,y = atomCoords[atom.handle] = atom.x, maxy-atom.y
        symbol = name_atom(atom)
        
        points.append((symbol, x, y))

    if maxy != None:
        maxy = maxy-miny
        miny = 0
    
    aromaticRings = {}

    minBondDist = 10e99
    
    for bond in bonds:
        atom1, atom2 = bond.atoms
        coord1, coord2 = atomCoords[atom1.handle], atomCoords[atom2.handle]

        if bond.stereo == "DOWN":
            x2, y2 = coord2
            x1, y1 = coord1
        else:
            x1, y1 = coord1
            x2, y2 = coord2
            
        dx = x1-x2
        dy = y1-y2
        d = math.sqrt(dx*dx+dy*dy)
        if d < minBondDist:
            minBondDist = d

        if bond.stereo:
            assert bond.bondtype == 1
            theta = -math.atan2(x2-x1, y2-y1)
            cosTHETA = math.cos(theta)
            sinTHETA = math.sin(theta)
            # XXX FIX ME, we should really know the
            #  font sizes here for proper drawing.
            if bond.atoms[1].symbol != "C":
                # shrink away from the endpoint
                offset = d * 0.8
            else:
                offset = d * 0.95

            if bond.atoms[0].symbol != "C":
                # shrink away from the start point
                dp = d * 0.1
                ox1, oy1 = x1 + dp*cosTHETA, y1 + dp*sinTHETA
            else:
                ox1, oy1 = x1, y1
                
            if bond.stereo == "UP":
                # make a filled polygon
                # coord1, coord2, coord3, filled 1/0
                # translate to "origin"

                length = d * 0.1

                wx1 = (length*cosTHETA - offset*sinTHETA) + x1
                wy1 = (length*sinTHETA + offset*cosTHETA) + y1

                wx2 = (-length*cosTHETA - offset*sinTHETA) + x1
                wy2 = (-length*sinTHETA + offset*cosTHETA) + y1
                wedges.append( (bond.bondorder, bond, (ox1,oy1), (wx1,wy1), (wx2,wy2), bond.stereo=="UP") )
            else:
                # assume bond stereo is down
                # we will make a cross hatch perpindicular to the bond
                for widthFactor in [0.2, 0.4, 0.6, 0.8, 1.0]:
                    newoffset = offset * widthFactor
                    length = 0.2 * newoffset
                    wx1 = (length*cosTHETA - newoffset*sinTHETA) + x1
                    wy1 = (length*sinTHETA + newoffset*cosTHETA) + y1
                    
                    wx2 = (-length*cosTHETA - newoffset*sinTHETA) + x1
                    wy2 = (-length*sinTHETA + newoffset*cosTHETA) + y1
                     
                    lines.append((bond.bondorder, bond, ((wx1, wy1), (wx2, wy2))))
            #lines.append((bond.bondorder, bond, (coord1, coord2)))
        elif not drawAromatic:
            lines.append((bond.bondorder, bond, (coord1, coord2)))
        else:
            if bond.bondtype == 4:
                atom1, atom2 = bond.atoms
                for ring in atom1.rings:
                    if not ring.aromatic:
                        continue
                    coords = aromaticRings.get(ring, {})
                    coords.update({coord1:ring})
                    aromaticRings[ring] = coords
                                        
                for ring in atom2.rings:
                    if not ring.aromatic:
                        continue
                    coords = aromaticRings.get(ring, {})
                    coords.update({coord2:ring})
                    aromaticRings[ring] = coords

            lines.append((bond.bondtype, bond, (coord1, coord2)))

    if minBondDist == 10e99:
        minBondDist = 1
    
    # XXX FIX Me, the circles can go outside
    # the lines occasionally
    circles = []
    for ring, ringAtoms in aromaticRings.items():
        coords = ringAtoms.keys()
        mx = 0.0
        my = 0.0
        for x, y in coords:
            mx += x
            my += y
        mx /= len(coords)
        my /= len(coords)

        # now get 90% of the smallest distance
        mindist = 10e99
        for x, y in coords:
            dx = x-mx
            dy = y-my
            dist = math.sqrt(dx*dx + dy*dy)
            if dist < mindist: mindist = dist

        circles.append((mx, my, mindist*.7))

    return points, lines, circles, wedges, (minx, maxx, miny, maxy), minBondDist
