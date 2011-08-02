"""
   This software is part of Frog, a chemo informatics class able to build 
   3D coordinates for small compounds
    Copyright (C) 2006-2007 P. Tuffery, B.O. Villoutreix, Th. Bohme Leite, D. Gomes, M. Miteva, J. Chomilier

    Frog2 (C) 2009-2010 by P. Tuffery, M. Miteva, F. Guyon

    Using this software, please cite:
        Frog2: Efficient 3D conformation ensemble generator for small compounds.
        Miteva MA, Guyon F, Tuffery P.
        Nucleic Acids Res. 2010 Jul;38(Web Server issue):W622-7. Epub 2010 May 5.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import sys
import math, time
from math import sqrt, fabs
import string
import copy
import iMolecule
#import iMolecule
from frowns import Atom, Bond, Molecule, Smiles
from frowns.Objects import FrownTypes
from frowns.Utils import Components
from frowns.Canonicalization import Disambiguate, Traverse
from frowns.Canonicalization.EquivClass import compute_equiv_class
from parametres import *
from Exceptions import MissingParameter


def getNbOfHInOatoms(atom):
	"""
	Returns the number of H attached to atom
	"""
	nbH = 0
	for at in atom.oatoms:
		if at.symbol == 'H':
			nbH += 1
	return nbH


def getAtom(handle, molecule):
	"""
	Returns an atom in molecule from his handle
	"""
	for at in molecule.atoms:
		if at.handle == handle:
			return at
	return None


def getIndexAtom(atom, molecule):
	"""
	Returns the index of an atom in molecule from his handle
	"""
	try:
		return int(getAtom(atom.handle,molecule).index)
	except:
		return None
	
	
def getAtomsVisitable(atomOri,bondToAvoid):
	"""
	Returns the atoms visitable from atomOri if 
	bondToAvoid is deleted from the molecule
	"""	
	def getAtomsVisitableIntern(atom,atomsVisited):
		atomsVisited.append(atom)
		for oatom in atom.oatoms:
			if oatom not in atomsVisited:
				getAtomsVisitableIntern(oatom,atomsVisited)
			
	atomsVisited = [atomOri]
	for oatom in atomOri.oatoms:
		if (oatom not in atomsVisited) and (atomOri.findbond(oatom) != bondToAvoid):
			getAtomsVisitableIntern(oatom,atomsVisited)
	return atomsVisited


def getDividedMoleculeAtoms(bondToAvoid):
	"""
	Returns the atoms visitable from bondToAvoid.atoms[0] and bondToAvoid.atoms[1] if 
	bondToAvoid is deleted from the molecule
	"""
	if bondToAvoid.rings:
		return bondToAvoid.parent.atoms
	
	atomsLeft = getAtomsVisitable(bondToAvoid.atoms[0],bondToAvoid)
	atomsRight = bondToAvoid.parent.atoms[:]
	for atom in atomsLeft:
		atomsRight.remove(atom)
	return [atomsLeft,atomsRight]


def toCanonicalChain(atomOri,bondToDel):
	mol = copy.deepcopy(atomOri.parent)
	ori = mol.atoms[atomOri.index]
	mol.remove_bond(mol.bonds[bondToDel.index])
	if len(ori.oatoms) == 0:
		return atomOri.symbol
	comp = Components.components(mol)

	if len(comp) == 1:
		atomList = comp[0][0]
	elif len(comp) == 2:
		if ori in comp[0][0]:
			atomList = comp[0][0]
		else:
			atomList = comp[1][0]
	
	mol.pruneToAtoms(atomList)
	Disambiguate.FreedDisambiguate(mol)
	
	max = -1
	for at in mol.atoms:
		if at.symorder > max: 
			max = at.symorder
	while ori.symorder != 0:
		for at in mol.atoms:
			if at.symorder == max: 
				at.symorder = 0
				at.symclass = 0
			else: 
				at.symorder += 1
				at.symclass += 1
	canonical, mol.canonical_list = Traverse.draw(mol)
	
	canonical = str(bondToDel.bondtype) + canonical
	mol.dirty = 0
	
	return canonical

def isRingedRing(atom, verbose = 0):
	"""
	Identifies atoms belonging to macro rings for which we do not
	manage chirality
	"""
	if verbose:
		sys.stderr.write("isRingedRing: considering %s %d\n" % (atom.symbol, atom.index))

	if (atom.aromatic) \
	or not (((atom.symbol in ['C','S','P']) and ((len(atom.oatoms) >= 4) and (getNbOfHInOatoms(atom) <= 1))) \
	               or ((atom.symbol == 'N') and ((len(atom.oatoms) >= 3) and (getNbOfHInOatoms(atom) <= 1)))):
		return False


	# Some special ring case
	if (len(atom.oatoms) == 4) and (atom.rings):
		nbOatomCycl = 0
		for oatom in atom.oatoms:
			if oatom.rings:
				nbOatomCycl += 1
		if nbOatomCycl >= 3:
			if verbose:
				sys.stderr.write("isRingedRing: Found unchiral atom\n")
			return True

	return False
	
def isChiralAtom(atom, verbose = 0):
	"""
	Determines is an atom is chiral
	"""

	if verbose:
		sys.stderr.write("isChiral: considering %s %d\n" % (atom.symbol, atom.index))
	# 1. Avoid some obvious cases not compatible with stereo isomery
	# Not proper atom type, etc
	if (atom.aromatic) \
	or not (((atom.symbol in ['C','S','P']) and ((len(atom.oatoms) >= 4) and (getNbOfHInOatoms(atom) <= 1))) \
	               or ((atom.symbol == 'N') and ((len(atom.oatoms) >= 3) and (getNbOfHInOatoms(atom) <= 1)))):
		return False

	# Some special ring case
	if (len(atom.oatoms) == 4) and (atom.rings):
		nbOatomCycl = 0
		for oatom in atom.oatoms:
			if oatom.rings:
				nbOatomCycl += 1
		if nbOatomCycl >= 3:
			return False	

	# N not in ring with 3 neighbours
	if (atom.number == 7) and (len(atom.oatoms) == 3) and (not atom.rings):
		for oat in atom.oatoms:
			if oat.symbol == 'C' and not oat.rings:
				for bond in oat.bonds:
					if bond.xatom(oat).symbol == 'O' and bond.bondtype == 2:
						return False
	# N in ring with less than 4 neighbours
	elif (atom.number == 7) and (len(atom.oatoms) <= 3) and (atom.rings):
		return False
	# N bound to aromatic C
	elif atom.number == 7:
		for oat in atom.oatoms:
			if oat.symbol == 'C' and oat.aromatic:
				return False

	# 2. Check stereo isomery possible
	neighbours = []
	neighbours[:] = atom.oatoms[:]
	neighbours.sort(key = lambda x: x.symbol)
	toDistinguish = []
	marked = []
	
	# Are all the immediate numbers different ?
	if verbose:
		sys.stderr.write("ischiral: Will consider neighbours:\n")
		for oatom in neighbours:
			sys.stderr.write("%s %d " % (oatom.symbol, oatom.index))
		sys.stderr.write("\n")
	
	for i in range(0,len(neighbours)-1):

		# If current not same atom type than next: OK
		if (neighbours[i] not in marked) and (neighbours[i].symbol != neighbours[i+1].symbol):
			pass

		# Same type as next one: we should investigate further
		elif neighbours[i] not in marked:
			marked.append(neighbours[i+1])
			toDistinguish.append((neighbours[i], neighbours[i].findbond(atom)))
			toDistinguish.append((neighbours[i+1], neighbours[i+1].findbond(atom)))
		elif (neighbours[i] in marked) and (neighbours[i].symbol == neighbours[i+1].symbol):
			marked.append(neighbours[i+1])
			toDistinguish.append((neighbours[i+1], neighbours[i+1].findbond(atom)))
			
	# All neigbours different: stereo site.
	if not toDistinguish:
		if verbose:
			sys.stderr.write("ischiral: True\n")
		return True

	# Are the cansmiles for the neighbors different ?
	chainsToDistinguish = []
	for atomToDistinguish in toDistinguish:
		chainsToDistinguish.append(toCanonicalChain(atomToDistinguish[0],atomToDistinguish[1]))
	
	chainsToDistinguish.sort()
	if verbose:
		sys.stderr.write("ischiral: will distinghuish\n")
		for aChain in chainsToDistinguish:
			sys.stderr.write("    %s\n" % aChain)
	
	for i in range(len(chainsToDistinguish)-1):
		if chainsToDistinguish[i] == chainsToDistinguish[i+1]:		
			return False	
	
	if verbose:
		sys.stderr.write("ischiral: True\n")
	return True

def isChiralBond(bond):
	if (bond.bondtype != 2) or (bond.rings):
		return False
	for atom in bond.atoms:
		if not ((atom.symbol in ['C','N','P']) and (((len(atom.oatoms) == 2) and (getNbOfHInOatoms(atom) == 0)) \
				                                    or ((len(atom.oatoms) == 3) and (getNbOfHInOatoms(atom) <= 1)))):
			return False
	for atom in bond.atoms:
		nDBonds = 0
		for bondd in atom.bonds:
			if bondd.bondtype in [2,4]:
				nDBonds += 1
		if nDBonds >= 2 and len(atom.oatoms) <=2:
			return False
	for atom in bond.atoms:
		toDistinguish = []
		for b in atom.bonds:
			if b != bond:
				toDistinguish.append(toCanonicalChain(atom,b))
		toDistinguish.sort()
		for i in range(len(toDistinguish)-1):
			if toDistinguish[i] == toDistinguish[i+1]:
				return False		
	bond.stereo = True
	return True

def neighboursAtomsOutsideCycle(atom):
	if not atom.rings:
		return []

	neighboursAtomsOutside = atom.oatoms[:]
	for cycle in atom.rings:
		for atom in cycle.atoms:
			try: 
				neighboursAtomsOutside.remove(atom)
			except:
				pass
	return neighboursAtomsOutside

def isAxEqValidPosition(oriCycle,atomsCycle,crds,AEdescript):
	#the first coordinates in crds must be for the atom that the axial/equatorial position stands for
	ori = oriCycle
	descript = AEdescript
							
	#Newell's algorithm for computing n, the normal for the best fit plane of the ori's cycle
	nx = ny = nz = 0.
	n = [0., 0., 0.]
	for k in range(len(atomsCycle)):
		if k == len(atomsCycle) - 1:
			kPlus1 = 0
		else:
			kPlus1 = k + 1
		n[0] += (atomsCycle[k].z + atomsCycle[kPlus1].z) * \
	            (atomsCycle[kPlus1].y - atomsCycle[k].y)
		n[1] += (atomsCycle[k].x + atomsCycle[kPlus1].x) * \
	            (atomsCycle[kPlus1].z - atomsCycle[k].z)
		n[2] += (atomsCycle[k].y + atomsCycle[kPlus1].y) * \
		        (atomsCycle[kPlus1].x - atomsCycle[k].x)
	for k in range(3):
		n[k] *= -1/2
	distN = sqrt(n[0]**2 + n[1]**2 + n[2]**2)
	for k in range(3):
		n[k] /= distN
	
	#check for axial/equatorial position
	c_nb = []
	if len(crds) < 2:
		return True
	for k in range(2):
		c_nb.append([crds[k][0] - ori.x, crds[k][1] - ori.y, crds[k][2] - ori.z])
		d = sqrt(c_nb[k][0]**2 + c_nb[k][1]**2 + c_nb[k][2]**2)
		for kk in range(3):
			c_nb[k][kk] /= d
	scalar = []
	for k in range(2):
		scalar.append(fabs(n[0]*c_nb[k][0] + n[1]*c_nb[k][1] + n[2]*c_nb[k][2]))
	#print "\n\nscalar", scalar,"\n\n"
	if ((scalar[0] > scalar[1]) and (descript == 'E')) or ((scalar[0] < scalar[1]) and (descript == 'A')):
		#print "\n\nFALSE isAxEqValidPosition\n\n"
		return False
	#print "\n\nTRUE isAxEqValidPosition\n\n"
	return True


def conc(a,b, nMin):
    """
    concatener 2 listes (qui ont au moins nMin elements en commun) en enlevant les doublons
    """
    c=[]
    for x in b :
        c.append(x)
    commun=0
    for atom in a :
        for ab in b :
            if atom==ab:
                commun+=1
                break                
        else:
            c.append(atom)
    if commun>=nMin :
        return c
    else :
        return None

def sortIndex(brickList) :
    """
    a partir d'une liste de liste d'atome -> liste de liste des indices tries de ces atomes
    """
    index = []
    res = []
    for fr in brickList :
        for atom in fr :
            index.append(atom.index)
        index.sort()
        res.append(index)
        index = []
    res.sort()
    
    return res


def coordsToSdf(iMol, mCoords, nrjs, fileName):
    ltime = time.localtime()
    _2ndLine = " Frog-RPBS"
    if ltime[1] < 10:
	_2ndLine += "0"
    _2ndLine += str(ltime[1])
    if ltime[2] < 10:
        _2ndLine += "0"
    _2ndLine += str(ltime[2]) + str(ltime[0])[2:5]
    if ltime[3] < 10:
        _2ndLine += "0"
    _2ndLine += str(ltime[3])
    if ltime[4] < 10:
        _2ndLine += "0"
    _2ndLine += str(ltime[4]) + "3D\n"
    _4thLine = "%3i%3i  0  0  0  0  0  0  0  0  1\n"%(len(iMol.molecule.atoms), len(iMol.molecule.bonds))
    bondLines = ""
    for bond in iMol.molecule.bonds:
        bondLines += "%3i%3i%3i%3i%3i%3i\n"%(bond.atoms[0].index+1, bond.atoms[1].index+1, bond.bondtype, 0, 0, 0)
    bondLines += "M  END\n>  <MMFF94 VDW and Electrostatic Energy>\n"
    end = "\n\n>  <SMILES>\n" + iMol.smiles + "\n\n"
    if iMol.axialEquatorial and (iMol.axialEquatorial.find("A") != -1 or iMol.axialEquatorial.find("E") != -1):
  	end += ">  <AX-EQ POS>\n" + iMol.axialEquatorial + "\n\n"
    end += "$$$$\n" 
    finLignesAtoms = []
    for atom in iMol.molecule.atoms:
        finLignesAtoms.append("%3s%3i%3i%3i%3i%3i\n"%(atom.symbol.title(), 0, 0, 0, 0, 0))

    file = open(fileName, "a")
    for index, coords in enumerate(mCoords):
        file.write(iMol.id + "_" + str(index+1) + "\n" + _2ndLine + "\n" + _4thLine)
        for index2, atCoords in enumerate(coords):
            file.write("%10.4f%10.4f%10.4f"%(atCoords[0], atCoords[1], atCoords[2]) + finLignesAtoms[index2])
        file.write(bondLines)
        file.write(nrjs[index])
	file.write(end)
    file.close()


def coordsToPdb(iMol, mCoords, nrjs, fileName):
    nbAtomsBySymbol = {}
    debutLignesAtoms = []
    finLignesAtoms = []
    for atom in iMol.molecule.atoms:
        if nbAtomsBySymbol.has_key(atom.symbol.title()):
            nbAtomsBySymbol[atom.symbol.title()] += 1
        else:
            nbAtomsBySymbol[atom.symbol.title()] = 1
        symbol = atom.symbol.title() + str(nbAtomsBySymbol[atom.symbol.title()])
	if atom.symbol != "H":
            debutLignesAtoms.append("ATOM" + str(atom.index + 1).rjust(7) + "  " + symbol.ljust(3) + " <0>     1    ")
	else:
            debutLignesAtoms.append("ATOM" + str(atom.index + 1).rjust(7) + " " + symbol.ljust(4) + " <0>     1    ")
        finLignesAtoms.append("  1.00  0.00" + atom.symbol.title().rjust(12) + "\n")

    bondLines = ""
    for atom in iMol.molecule.atoms:
        bondLines += "CONECT" + str(atom.index+1).rjust(5)
        for oatom in atom.oatoms:
           bondLines += str(oatom.index + 1).rjust(5)
        bondLines += "\n"
    bondLines += "MASTER        1    0    0    0    0    0    0    0" + str(len(iMol.molecule.atoms)).rjust(5) + \
                 "    0" + str(len(iMol.molecule.atoms)).rjust(5) + "    0\n"

    file = open(fileName, "a")
    for index, coords in enumerate(mCoords):
        file.write("HEADER    " + iMol.id + "_" + str(index+1) + "\nCOMPND    MOL_ID: 1;\nCOMPND   2 MOLECULE: " + iMol.smiles)
	if iMol.axialEquatorial and (iMol.axialEquatorial.find("A") != -1 or iMol.axialEquatorial.find("E") != -1):
	    file.write(" " + iMol.axialEquatorial)
	file.write("\nAUTHOR    generated by FROG v1.01, RPBS-TEAM\nREMARK   6 MMFF94 VDW and Electrostatic Energy = " + nrjs[index] + "\n")
        for index2, atCoords in enumerate(coords):
            file.write(debutLignesAtoms[index2] + str("%8.3f"%(atCoords[0])).rjust(8) + \
                       str("%8.3f"%(atCoords[1])).rjust(8) + str("%8.3f"%(atCoords[2])).rjust(8) + finLignesAtoms[index2])
        file.write(bondLines)
    file.close()
  

def coordsToMol2(iMol, mCoords, nrjs, fileName):
    nbAtomsBySymbol = {}
    debutLignesAtoms = []
    atomstypes = []
    for atom in iMol.molecule.atoms:
        if nbAtomsBySymbol.has_key(atom.symbol.title()):
            nbAtomsBySymbol[atom.symbol.title()] += 1
        else:
            nbAtomsBySymbol[atom.symbol.title()] = 1
	atomstypes.append(getMol2AtomType(atom))
	debutLignesAtoms.append(str(atom.index + 1).rjust(7) + " " + str(atom.symbol.title() + str(nbAtomsBySymbol[atom.symbol.title()])).ljust(8) + "")

    lines3to6 = str(len(iMol.molecule.atoms)).rjust(5) + str(len(iMol.molecule.bonds)).rjust(6) + "     0     0     0\nSMALL\nNO_CHARGES\n\n"
    bondLines = "@<TRIPOS>BOND\n"
    for index, bond in enumerate(iMol.molecule.bonds):
        bondLines += str(index+1).rjust(6)
	if bond.atoms[0].index < bond.atoms[1].index:
	    bondLines += str(bond.atoms[0].index+1).rjust(5) + str(bond.atoms[1].index+1).rjust(5) + " "
	else:
            bondLines += str(bond.atoms[1].index+1).rjust(5) + str(bond.atoms[0].index+1).rjust(5) + " "
        if bond.bondtype == 4:
	    bondLines += "ar\n"
	else:
	    bondLines += str(bond.bondtype) + "\n"

    file = open(fileName, "a")
    for index, coords in enumerate(mCoords):
        file.write("@<TRIPOS>MOLECULE\n" + iMol.id + "_" + str(index+1) + "\n" + lines3to6)
        file.write("VDW energy = " + nrjs[index] + "\n@<TRIPOS>ATOM\n")
	for index2, atCoords in enumerate(coords):
 	    file.write(debutLignesAtoms[index2] + str("%.4f"%(atCoords[0])).rjust(10) + \
                       str("%.4f"%(atCoords[1])).rjust(10) + str("%.4f"%(atCoords[2])).rjust(10) + " " + \
		       atomstypes[index2].ljust(5) + "     1 <O>         0.0000\n")
	file.write(bondLines + "\n")
    file.close()

def getSaveIndex(saveindex, handle):
	if saveindex == None:
		return None
	try:
		return saveindex[handle]
	except:
		return None

def circular(alist):
	"""
	Tool for cirucular permutation of a list (1 step)
	"""
	rs = [alist[-1]]
	rs =  rs + alist[:-1]
	return rs

chiral_map = {0: [[2,3,4], [3,4,2], [4,2,3]],
	      1: ["Ring case!"],
	      2: [[0,4,3], [4,3,0], [3,0,4]],
	      3: [[0,2,4], [2,4,0], [4,0,2]],
	      4: [[0,3,2], [3,2,0], [2,0,3]]}

def isringsamesens(l1,l2, verbose = 0):
	sys.stderr.write("isringsamesens %s %s\n" % (l1,l2) )
	sys.exit(0)

def issamesens(l1,l2, verbose = 0):
	"""
	Check if the way we screen atoms is compatible with smiles chiral order
	as wanted.
	l1 is smiles order (the reference)
	l2 is present order (the proposed order)
	Will return False only if sure order are not compatible.
eval $cmd 2>&1 | grep issam

	"""
	if len(l1) < 5:
		return True
	a = copy.deepcopy(l1)
	b = copy.deepcopy(l2)
	a.sort()
	b.sort()
	# Not same atoms to compare, we pass
	if a != b:
		return True

	amap = []
	for l in l2:
		amap.append(l1.index(l))
	# print l1, l2, amap
	
	pos = l1.index(l2[0]) # 1er atome a placer presentement
	if verbose:
		sys.stderr.write("issamesens %s %s. %s pos is %d. %s vs %s\n" % (str( l1), str(l2), amap, pos, str(amap[2:]), str(chiral_map[amap[0]])))

	# This sounds strange. To check.
	if amap[0] == 1:
		return True
	if amap[2:] in chiral_map[amap[0]]:
		if verbose:
			sys.stderr.write("issamesens: True\n")
		return True
	if verbose:
		sys.stderr.write("issamesens: False\n")
	return False

def getChiralNeighbours(ori, neighboursPlaced, neighboursToPlace, saveindex, verbose = 0):
    """
    ori: chiral center (Atom)
    neighboursPlaced: neighbours with 3D assigned (list of Atom)
    neighboursToPlace: neighbours with 3D unassigned (list of Atom)
    
    returns:
    smilesorder: the index off the atoms in the smiles
    currentorder: the present order to generate 3D
    """
    smilesorder  = None
    currentorder = None
    hsaveindex   = None
    smilesorder = []
    currentorder = []
    for a in neighboursPlaced:
	    smilesorder.append(getSaveIndex(saveindex, a.handle))
	    currentorder.append(getSaveIndex(saveindex, a.handle))
	    if a.symbol == "H":
		    hsaveindex = getSaveIndex(saveindex, a.handle)
    smilesorder.append(getSaveIndex(saveindex, ori.handle))
    currentorder.append(getSaveIndex(saveindex, ori.handle))
    for a in neighboursToPlace:
	    smilesorder.append(getSaveIndex(saveindex, a.handle))
	    currentorder.append(getSaveIndex(saveindex, a.handle))
	    if a.symbol == "H":
		    hsaveindex = getSaveIndex(saveindex, a.handle)
    smilesorder.sort()
    if hsaveindex != None:
	    pos = smilesorder.index(hsaveindex)
	    l2 = smilesorder[:2]
	    l2.append(hsaveindex)
	    for l in smilesorder[2:]:
		    if l not in l2:
			    l2.append(l)
	    smilesorder = l2
    if verbose:
	    sys.stderr.write("setCoord: chiral smiles order is %s\n" % str(smilesorder))
	    sys.stderr.write("setCoord: current build order is %s\n" % str(currentorder))
	    sys.stderr.write("setCoord: is order compatible: %s\n" % str(issamesens(smilesorder, currentorder, verbose = verbose)))
    
    return smilesorder, currentorder

def isRingChiralSens(ori, neighboursPlaced, neighboursToPlace, saveindex, verbose = 1):
    """
    ori: chiral center (Atom)
    neighboursPlaced: neighbours with 3D assigned (list of Atom)
    neighboursToPlace: neighbours with 3D unassigned (list of Atom)
    
    returns:
    smilesorder: the index of the atoms in the smiles
    currentorder: the present order to generate 3D

    This supplementary function comes since the way smiles manages stereochemistry in rings differs from that of aliphatic centers.
    Convention about stereo centers in rings is different from linear ones:
    It is related to the ring closure, as in smiles.
    Convention unclear from daylight' doc.
    Start from H (or atom) towards ring chiral atom, then non ring atom, 
    then ring atoms

    """
    smilesplacedorder   = []
    currentplacedorder  = []
    currentplacedlabel  = []
    currenttoplaceorder = []
    currenttoplacelabel = []
    currentlabel = []
    hsaveindex   = None

    # Identify ring
    for cycle in ori.rings:
	    OK = True
	    for a in neighboursPlaced:
		    if a not in cycle.atoms:
			    OK = False
			    break
    ringindex = []
    for a in cycle.atoms:
	    ringindex.append(getSaveIndex(saveindex, a.handle))
    ringindex.sort()
    isRingClosure = getSaveIndex(saveindex,ori.handle) in (ringindex[0], ringindex[-1])
    isEnd = getSaveIndex(saveindex,ori.handle) == ringindex[0]

    if not isRingClosure:
	    # We follow standard rule
	    if verbose:
		    sys.stderr.write("isRingChiralSens:  Not a ring closure\n")
	    smilesorder  = [getSaveIndex(saveindex,neighboursPlaced[0].handle), getSaveIndex(saveindex,ori.handle), getSaveIndex(saveindex,neighboursPlaced[1].handle)]
	    if smilesorder[0] > smilesorder[2]:
		    b = smilesorder[2]
		    smilesorder[2] = smilesorder[0]
		    smilesorder[0] = b
	    currentorder = [getSaveIndex(saveindex,neighboursPlaced[0].handle), getSaveIndex(saveindex,ori.handle), getSaveIndex(saveindex,neighboursPlaced[1].handle)]
	    placeorder = []
	    for a in neighboursToPlace:
		    if a.symbol == "H":
			    hsaveindex = getSaveIndex(saveindex, a.handle)
		    placeorder.append(getSaveIndex(saveindex,a.handle))
		    currentorder.append(getSaveIndex(saveindex,a.handle))
	    placeorder.sort()
	    smilesorder = smilesorder + placeorder
	    if hsaveindex != None:
		    pos = smilesorder.index(hsaveindex)
		    l2 = smilesorder[:3]
		    l2.append(hsaveindex)
		    for l in smilesorder[3:]:
			    if l not in l2:
				    l2.append(l)
		    smilesorder = l2
 
	    rs  = issamesens(smilesorder,currentorder, verbose = verbose)

	    if verbose:
		    sys.stderr.write("isRingChiralSens: chiral smiles order is %s\n" % str(smilesorder))
		    sys.stderr.write("isRingChiralSens: chiral current order is %s\n" % str(currentorder))
		    sys.stderr.write("isRingChiralSens: found %s\n" % str(rs))
	    # Why not rs is unclear. To check on more cases
	    return not rs

    else: # ring closure, more tricky
	    if verbose:
		    sys.stderr.write("isRingChiralSens:  Ring %s closure Ring atoms: %s ori: %s\n" % ( ("end", "start")[isEnd], str(ringindex), str(getSaveIndex(saveindex,ori.handle)) ) )
	    placeorder = []
	    if verbose:
		    sys.stderr.write("isRingChiralSens: smilesorder label0\n" )

	    for a in neighboursPlaced:
		    placeorder.append(getSaveIndex(saveindex,a.handle))
	    placeorder.sort()
	    if verbose:
		    sys.stderr.write("isRingChiralSens: smilesorder label1 %d %d\n" % (placeorder[0],placeorder[1]) )
	    smilesorder  = [placeorder[0], getSaveIndex(saveindex, ori.handle), placeorder[1]]
	    if verbose:
		    sys.stderr.write("isRingChiralSens: smilesorder label2\n" )
	    currentorder = [getSaveIndex(saveindex,neighboursPlaced[0].handle), getSaveIndex(saveindex,ori.handle), getSaveIndex(saveindex,neighboursPlaced[1].handle)]
	    placeorder = []
	    for a in neighboursToPlace:
		    if a.symbol == "H":
			    hsaveindex = getSaveIndex(saveindex, a.handle)
		    placeorder.append(getSaveIndex(saveindex,a.handle))
		    currentorder.append(getSaveIndex(saveindex,a.handle))
	    placeorder.sort()
	    if verbose:
		    sys.stderr.write("isRingChiralSens: smilesorder label3\n" )
	    smilesorder = smilesorder + placeorder
	    if hsaveindex != None:
		    pos = smilesorder.index(hsaveindex)
		    l2 = smilesorder[:3]
		    l2.append(hsaveindex)
		    for l in smilesorder[3:]:
			    if l not in l2:
				    l2.append(l)
		    smilesorder = l2
	    if verbose:
		    sys.stderr.write("isRingChiralSens:  will issamesens\n" )
	    rs  = issamesens(smilesorder,currentorder, verbose = verbose)
	    if verbose:
		    sys.stderr.write("isRingChiralSens: chiral smiles order is %s\n" % str(smilesorder))
		    sys.stderr.write("isRingChiralSens: chiral current order is %s\n" % str(currentorder))
		    sys.stderr.write("isRingChiralSens: found %s\n" % str(rs))

	    return not rs
    
    # sys.exit(0)
    # Look for atoms outside ring
    outCycle = []
    outCyleAtoms = neighboursAtomsOutsideCycle(ori)
    for a in outCyleAtoms:
	    # outCycle.append(getSaveIndex(saveindex, a.handle))
	    idex = getSaveIndex(saveindex, a.handle)
	    if a.symbol == "H":
	    	    outCycle.insert(0, idex)
	    else:
	    	    outCycle.append(idex)

    # Necessarily, Frog knows about ring atom coordinates
    for a in neighboursPlaced:
	    smilesplacedorder.append(getSaveIndex(saveindex, a.handle))
	    currentplacedorder.append(getSaveIndex(saveindex, a.handle))
	    currentplacedlabel.append(a.symbol)

    smilesplacedorder.sort()

    for a in neighboursToPlace:
	    currenttoplaceorder.append(getSaveIndex(saveindex, a.handle))
	    currenttoplacelabel.append(a.symbol)

    if verbose:
	    sys.stderr.write("isRingChiralSens: outside ring %s\n" % str(outCycle))
	    sys.stderr.write("isRingChiralSens: chiral sens %s\n" % str(ori.chiral_class))
	    sys.stderr.write("isRingChiralSens: chiral smiles order is %s\n" % str(smilesplacedorder))
	    sys.stderr.write("isRingChiralSens: current build order is %s %s\n" % (str(currentplacedorder), str(currenttoplaceorder)))
	    sys.stderr.write("isRingChiralSens: current symbol order is %s %s\n" % (str(currentplacedlabel), str(currenttoplacelabel)))

    # smilesplacedorder is ring closure order
    # Do we see it this way ? if yes, then @ stands for outCycle == currenttoplaceorder
    # if no then @ stands for outCycle != currenttoplaceorder
    if smilesplacedorder == currentplacedorder:  # OK on daylight small test C[C@@H]1CCCCO1 test
	    if ori.chiral_class == 1: 
		    if outCycle == currenttoplaceorder:
			    rs = False
		    rs = True
	    else:
		    if outCycle == currenttoplaceorder:
			    rs = False
		    rs = True
    else: 
	    if ori.chiral_class == 1: # OK on O1CCCC[C@H]1C test
		    if outCycle != currenttoplaceorder:
			    rs = False
		    rs = True
	    else: # OK on O1CCCC[C@@H]1C test
		    if outCycle != currenttoplaceorder:
			    rs = False
		    rs = True
    if verbose:
	    sys.stderr.write("isRingChiralSens: found %s\n" % str(rs))
	    
    return rs
	    

def setCoord(ori, neighboursPlaced, neighboursToPlace, nbOfGhostsToPlace=0, sens=1, saveindex = None, verbose = 0) :
    """
    setCoord: will setup coordinates for one atom given its neighgours

    @param ori              : atom instance
    @param neighboursPlaced : atom list
    @param neighboursToPlace: atom list
    @param sens             : chiral sense, one of 1, -1

    @return: a list of coordinates for the neighboursToPlace
    """

    if verbose:
	    sys.stderr.write("setCoord: ori is %d-%d-%s %s (%d to place) %s %d ghosts\n" % (int(ori.index), int(ori.handle), str(getSaveIndex(saveindex, ori.handle)), str(ori.symbol), len(neighboursToPlace), str(neighboursToPlace), nbOfGhostsToPlace) )
	    sys.stderr.write("setCoord Placed: ")
	    for a in neighboursPlaced:
		    # print a.index, a.symbol
		    sys.stderr.write("%s-%s-%s %s " % (str(a.index), str(a.handle), str(getSaveIndex(saveindex, a.handle)),str(a.symbol)) )
	    sys.stderr.write("\n")
	    sys.stderr.write("setCoord toPlace: ")
	    for a in neighboursToPlace:
		    # print a.index, a.symbol
		    sys.stderr.write("%s-%s-%s %s " % (str(a.index), str(a.handle), str(getSaveIndex(saveindex, a.handle)), str(a.symbol)) )
	    sys.stderr.write("\n")
	    if len(neighboursToPlace) > 1:
		    sys.stderr.write("setCoord: ... (%d to place) %s\n" % (len(neighboursToPlace), str(neighboursToPlace)) )
		    sys.stderr.write("setCoord: issamesens chiral_class %s\n" % ori.chiral_class)
		    sys.stderr.write("setCoord: ori %d type %d chiral %s toplace %d (type %d %d ...) chiral sense %d\n" % (ori.index, ori.number, ori.chiral_class, len(neighboursToPlace), neighboursToPlace[0].number, neighboursToPlace[1].number, sens))
	    elif len(neighboursToPlace):
		    sys.stderr.write("setCoord: ori type %d toplace %d (type %d) chiral sense %d\n" % (ori.number, len(neighboursToPlace), neighboursToPlace[0].number, sens))
	    else:
		    sys.stderr.write("setCoord: Nothing to do\n")
		    
    smilesorder  = None
    currentorder = None
    hsaveindex   = None
#     if ori.chiral_class:
# 	    smilesorder = []
# 	    currentorder = []
# 	    for a in neighboursPlaced:
# 		    smilesorder.append(getSaveIndex(saveindex, a.handle))
# 		    currentorder.append(getSaveIndex(saveindex, a.handle))
# 		    if a.symbol == "H":
# 			    hsaveindex = getSaveIndex(saveindex, a.handle)
# 	    smilesorder.append(getSaveIndex(saveindex, ori.handle))
# 	    currentorder.append(getSaveIndex(saveindex, ori.handle))
# 	    for a in neighboursToPlace:
# 		    smilesorder.append(getSaveIndex(saveindex, a.handle))
# 		    currentorder.append(getSaveIndex(saveindex, a.handle))
# 		    if a.symbol == "H":
# 			    hsaveindex = getSaveIndex(saveindex, a.handle)
# 	    smilesorder.sort()
# 	    if hsaveindex != None:
# 		    pos = smilesorder.index(hsaveindex)
# 		    l2 = smilesorder[:2]
# 		    l2.append(hsaveindex)
# 		    for l in smilesorder[2:]:
# 			    if l not in l2:
# 				    l2.append(l)
# 		    smilesorder = l2
# 	    if verbose:
# 		    sys.stderr.write("setCoord: chiral smiles order is %s\n" % str(smilesorder))
# 		    sys.stderr.write("setCoord: current build order is %s\n" % str(currentorder))
# 		    sys.stderr.write("setCoord: is order compatible: %s\n" % str(issamesens(smilesorder, currentorder, verbose = verbose)))
# 	    # sys.exit(0)
# # setCoord: chiral smiles order is [6, 7, 42, 8, 14]
# # setCoord: current build order is [8, 7, 6, 42, 14]
    
    if ori.chiral_class:
	    if ori.rings:
		    # P. TUFFERY 2011 BUGFIX: added  (len(neighboursPlaced) > 1).
		    if (len(neighboursPlaced) > 1) and (not isRingChiralSens(ori, neighboursPlaced, neighboursToPlace, saveindex, verbose = verbose)):
			    sens *= -1
	    else:
		    smilesorder, currentorder = getChiralNeighbours(ori, neighboursPlaced, neighboursToPlace, saveindex, verbose = verbose)
		    if not issamesens(smilesorder, currentorder):
			    sens *= -1
	    

    A=[0.,0.,0.]
    B=[0.,0.,0.]
    C=[0.,0.,0.]
    crds = []

    if verbose:
	    sys.stderr.write("setCoords: will proceed\n")

    if not neighboursToPlace:
	    # sys.stderr.write("setCoords: nothing to do\n")
	    return crds
    #elif len(v)+len(nb)+ori.hcount==3 and (ori.number==6 or ori.number==7) : # Trigonal (Plan) cxx ou nxx
    #elif (len(v)+len(nb)+ori.imp_hcount==3) and (ori.number==6 or ori.number==7) : #thiagggg
    # sp2 case
    elif (len(neighboursPlaced) + len(neighboursToPlace) + nbOfGhostsToPlace == 3) and (ori.number in [6,7,16]) :
	    # sys.stderr.write("nbOfGhostsToPlace %d\n" % nbOfGhostsToPlace)
	    if len(neighboursPlaced)==2 : # a on cxx or nxx  
		    if verbose:
			    sys.stderr.write("setCoords: a on cxx or nxx\n")
		    crds.append([0,0,0])
		    # sys.stderr.write("setCoords: a on cxx or nxx (2)\n")
		    rh = dist(ori,neighboursToPlace[0]) - 0.01
		    # sys.stderr.write("setCoords: a on cxx or nxx (2)\n")
		    th = sens*(2.*math.pi - angle(neighboursPlaced[0],ori,neighboursPlaced[1]))/2.
		    # sys.stderr.write("setCoords: a on cxx or nxx (3)\n")
		    tau = math.pi
		    cosdir(ori,neighboursPlaced[1],neighboursPlaced[0],A,B,C)
		    # sys.stderr.write("setCoords: a on cxx or nxx (3)\n")
		    xn = rh * math.cos(th)
		    yn = rh * math.sin(th) * math.cos(tau)
		    zn = -rh * math.sin(th) * math.sin(tau)
		    crds[0][0] = A[0] * xn + A[1] * yn + A[2] * zn + ori.x
		    crds[0][1] = B[0] * xn + B[1] * yn + B[2] * zn + ori.y
		    crds[0][2] = C[0] * xn + C[1] * yn + C[2] * zn + ori.z
		    # sys.stderr.write("setCoords: a on cxx or nxx (done)\n")
	    elif len(neighboursPlaced)==1 : # a2 on cx -> a2 on nx
		    if verbose:
			    sys.stderr.write("setCoords: a2 on cx -> a2 on nx\n")
		    for i in range(2 - nbOfGhostsToPlace) :
			    crds.append([0,0,0])
		    for i2 in neighboursPlaced[0].oatoms :
			    if i2 == ori :
				continue
			    elif i2.x != None :
				break
		    # Repere orthonorme dont l'axe X est le vecteur ori-neighboursPlaced[0]
		    cosdir(ori,neighboursPlaced[0],i2,A,B,C)
		    th  = math.radians(120)
		    if sens == 1:
			    tau = -math.pi # Plug sens here !! since cos(tau) gives sign !!
		    else:
			    tau = 0.
		    for i in range(2 - nbOfGhostsToPlace) :
			    if i+1>len(neighboursToPlace) :
				rh = 1.09
			    else :
				rh = dist(ori,neighboursToPlace[i])
			    xn = rh * math.cos(th)
			    yn = rh * math.sin(th) * math.cos(tau)
			    zn = -rh * math.sin(th) * math.sin(tau)
			    crds[i][0] = A[0] * xn + A[1] * yn + A[2] * zn + ori.x
			    crds[i][1] = B[0] * xn + B[1] * yn + B[2] * zn + ori.y
			    crds[i][2] = C[0] * xn + C[1] * yn + C[2] * zn + ori.z
			    tau += math.pi
				    
	    # elif len(v)+len(nb)+ori.hcount==4 : #and (ori.number==6 or ori.number==7) : # carbone/azote tetragonal
	    # elif (len(v)+len(nb)+ori.imp_hcount==4): #thiagggggg

    #
    # sp3 case
    #
    elif len(neighboursPlaced) + len(neighboursToPlace)  + nbOfGhostsToPlace == 4:
	    if len(neighboursPlaced)==3 : # a on cxxx
		    if verbose:
			    sys.stderr.write("setCoords: a on cxxx\n")
		    # sys.stderr.write("setCoords: a on cxxx\n")
		    crds.append([0,0,0])
		    # Vecteur atm-i1
		    # sys.stderr.write("setCoords: a on cxxx (0) %s %s\n" % (str(neighboursPlaced[0]), str(neighboursPlaced[0].x)))
		    A[0] = neighboursPlaced[0].x - ori.x
		    A[1] = neighboursPlaced[0].y - ori.y
		    A[2] = neighboursPlaced[0].z - ori.z
		    # sys.stderr.write("setCoords: a on cxxx (1)\n")
		    ar = math.sqrt(A[0]*A[0] + A[1]*A[1] + A[2]*A[2])
		    # sys.stderr.write("Yes\n")
		    # Vecteur atm-i2
		    B[0] = neighboursPlaced[1].x - ori.x
		    B[1] = neighboursPlaced[1].y - ori.y
		    B[2] = neighboursPlaced[1].z - ori.z
		    br = math.sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2])
		    # Vecteur atm-i3
		    C[0] = neighboursPlaced[2].x - ori.x
		    C[1] = neighboursPlaced[2].y - ori.y
		    C[2] = neighboursPlaced[2].z - ori.z
		    cr = math.sqrt(C[0]*C[0] + C[1]*C[1] + C[2]*C[2])
		    # Vecteur C = atm-i2 - atm-i1 (normes)
		    dx = B[0]/br - A[0]/ar
		    dy = B[1]/br - A[1]/ar
		    dz = B[2]/br - A[2]/ar
		    # Vecteur D = atm-i3 - atm-i1 (normes)
		    ex = C[0]/cr - A[0]/ar
		    ey = C[1]/cr - A[1]/ar
		    ez = C[2]/cr - A[2]/ar
		    # Vecteur CxD (Pointe dans la bonne direction)
		    B[0] = dy * ez - dz * ey
		    B[1] = dz * ex - dx * ez
		    B[2] = dx * ey - dy * ex
		    dot = A[0]*B[0] + A[1]*B[1] + A[2]*B[2]
		    br = math.sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2])
		    if dot < 0. :
			    epsi = 1.
		    else :
			    epsi = -1.
		    fac = epsi*dist(ori,neighboursToPlace[0])/br
		    crds[0][0] = ori.x + B[0] * fac
		    crds[0][1] = ori.y + B[1] * fac
		    crds[0][2] = ori.z + B[2] * fac
		    # sys.stderr.write("setCoords: a on cxxx (done)\n")
	    elif len(neighboursPlaced)==2 : # a2 on cxx
		    if verbose:
			    sys.stderr.write("setCoords: a2 on cxx %d \n" % (nbOfGhostsToPlace) )
		    for i in range(2 - nbOfGhostsToPlace) :
			    crds.append([0,0,0])
		    # Vecteur atm-i1
		    A[0] = neighboursPlaced[0].x - ori.x
		    A[1] = neighboursPlaced[0].y - ori.y
		    A[2] = neighboursPlaced[0].z - ori.z
		    ar = math.sqrt(A[0] * A[0] + A[1] * A[1] + A[2] * A[2])
		    # Vecteur atm-i2
		    B[0] = neighboursPlaced[1].x - ori.x
		    B[1] = neighboursPlaced[1].y - ori.y
		    B[2] = neighboursPlaced[1].z - ori.z
		    br = math.sqrt(B[0] * B[0] + B[1] * B[1] + B[2] * B[2])
		    # sys.stderr.write("Yes\n")
		    # Vecteur C = atm-i2 - atm-i1 (normes)
		    C[0] = -A[0] / ar - B[0] / br
		    C[1] = -A[1] / ar - B[1] / br
		    C[2] = -A[2] / ar - B[2] / br
		    cr = math.sqrt(C[0] * C[0] + C[1] * C[1] + C[2] * C[2])
		    d0 = dist(ori,neighboursToPlace[0])
		    fach = d0 * math.cos(math.radians(109.5)/2) / cr
		    # Vecteur D = atm-i1 x atm-i2
		    dx = A[1] * B[2] - A[2] * B[1]
		    dy = A[2] * B[0] - A[0] * B[2]
		    dz = A[0] * B[1] - A[1] * B[0]
		    dr = math.sqrt(dx * dx + dy * dy + dz * dz)
		    # - sens stands for compatibility with sp2 case and for counter clockwise for @
		    facv = - sens * d0 * math.sin(math.radians(109.5)/2) / dr

		    crds[0][0] = ori.x + C[0] * fach + dx * facv
		    crds[0][1] = ori.y + C[1] * fach + dy * facv
		    crds[0][2] = ori.z + C[2] * fach + dz * facv
		    if (verbose):
			    sys.stderr.write("crds0 %f %f %f\n" % (crds[0][0], crds[0][1], crds[0][2])  )
		    if nbOfGhostsToPlace == 0:
			    d1 = dist(ori,neighboursToPlace[1])
			    if d0 == d1:
				    crds[1][0] = ori.x + C[0] * fach - dx * facv
				    crds[1][1] = ori.y + C[1] * fach - dy * facv
				    crds[1][2] = ori.z + C[2] * fach - dz * facv
			    else:
				    crds[1][0] = ori.x + (C[0] * fach - dx * facv) * (d1 / d0)
				    crds[1][1] = ori.y + (C[1] * fach - dy * facv) * (d1 / d0)
				    crds[1][2] = ori.z + (C[2] * fach - dz * facv) * (d1 / d0)
			    if (verbose):
				    sys.stderr.write("crds1 %f %f %f\n" % (crds[1][0], crds[1][1], crds[1][2])  )

	    elif len(neighboursPlaced)==1 : # a3 on cx
		    if verbose:
			    sys.stderr.write("setCoords: a3 on cx %d\n" % (nbOfGhostsToPlace) )
		    # sys.stderr.write("setCoords: a3 on cx\n")
		    for i in range(3 - nbOfGhostsToPlace) :
			    crds.append([0,0,0])
		    for i2 in neighboursPlaced[0].oatoms :
			    if i2 == ori :
				    continue
			    elif i2.x != None:
				    break
		    # Build an ortho,normalized reference on ori
		    cosdir(ori,neighboursPlaced[0],i2,A,B,C)
		    th = math.radians(109.47)
		    tau = math.pi
		    for i in range(3 - nbOfGhostsToPlace) :
			    if i+1>len(neighboursToPlace) :
				    rh = 1.09
			    else :
				    rh = dist(ori,neighboursToPlace[i])
			    xn = rh * math.cos(th)
			    yn = rh * math.sin(th) * math.cos(tau)
			    zn = -rh * math.sin(th) * math.sin(tau)
			    crds[i][0] = A[0] * xn + A[1] * yn + A[2] * zn + ori.x
			    crds[i][1] = B[0] * xn + B[1] * yn + B[2] * zn + ori.y
			    crds[i][2] = C[0] * xn + C[1] * yn + C[2] * zn + ori.z
			    # tau -= math.radians(120)
			    tau = tau - sens * math.radians(120)
			    if verbose:
				    sys.stderr.write("tau is %f\n" % tau)
            
    #elif len(v)+len(nb)+ori.hcount==2 and (ori.number==8 or ori.number==16 or ori.number==7 or ori.number==6) : # x-o-  x-s- x=n- x#c-
    #elif len(v)+len(nb)+ori.imp_hcount==2 and (ori.number==8 or ori.number==16 or ori.number==7 or ori.number==6) : #thiagggggg
    elif (len(neighboursPlaced) + len(neighboursToPlace) == 2) and (ori.number in [6,7,8,16,17]):
        if verbose:
		sys.stderr.write("special linear case for type %d linked to %d Will place type %d\n" % (ori.number, neighboursPlaced[0].number, neighboursToPlace[0].number))
        #print "in set coord: X-O- ou X-S- ou X=N- ou X#C-", ori, neighboursPlaced, neighboursToPlace
        #print "ori .x: ", ori.x, ".y: ", ori.y, ".z: ", ori.z
        #print "neighboursPlaced[0].x: ",neighboursPlaced[0].x, ".y: ", neighboursPlaced[0].y, ".z: ", neighboursPlaced[0].z
        # sys.stderr.write("setCoords: special linear\n")
        crds.append([0,0,0])
        if ori.number == 8 :
            # th = math.radians(105.)
            th = math.radians(109.5)
        elif ori.number == 16 :
            th = math.radians(95.)
        #elif ori.number == 7 :
        #    th = math.radians(120)
        #else :
        #    print "\n\n OIUI\n\n"
        #    th = math.pi
        else:
            rh = dist(ori,neighboursToPlace[0])
            distt = dist(ori,neighboursToPlace[0])
            crds[0][0] = ori.x + ((ori.x-neighboursPlaced[0].x)/distt)*rh
            crds[0][1] = ori.y + ((ori.y-neighboursPlaced[0].y)/distt)*rh
            crds[0][2] = ori.z + ((ori.z-neighboursPlaced[0].z)/distt)*rh
            return crds
        tau = math.pi
        rh = dist(ori,neighboursToPlace[0])
        for i2 in neighboursPlaced[0].oatoms :
            if i2 == ori :
                continue
            elif i2.x != None :
                #print "YO"
                break

        cosdir(ori,neighboursPlaced[0],i2,A,B,C)
        xn = rh * math.cos(th)
        yn = rh * math.sin(th) * math.cos(tau)
        zn = -rh * math.sin(th) * math.sin(tau)

        crds[0][0] = A[0] * xn + A[1] * yn + A[2] * zn + ori.x
        crds[0][1] = B[0] * xn + B[1] * yn + B[2] * zn + ori.y
        crds[0][2] = C[0] * xn + C[1] * yn + C[2] * zn + ori.z

    else:
	sys.stderr.write("no geometry parameters for atom %s\n" % ori.symbol)
        raise Exception("no geometry parameters for atom " + ori.symbol)
    return crds

def cosdir(i,j,k,A,B,C) :
	"""
	@param i,j,k : 3 atoms
	@param A,B,C : 3 atoms (output), a reference based on i
	[0] = VECTEUR i-j
	[1] = (i-j ^ i-k) ^ i-j
	[2] = i-j ^ i-k
	"""
	# Vecteur ij norme
	A[0] = j.x - i.x
	B[0] = j.y - i.y
	C[0] = j.z - i.z
	pm = math.sqrt(A[0]*A[0] + B[0]*B[0] + C[0]*C[0])
	A[0] /= pm
	B[0] /= pm
	C[0] /= pm

	# Vecteur jk
	a = k.x - j.x
	b = k.y - j.y
	c = k.z - j.z

	# Vecteur ij^jk norme
	A[2] = B[0] * c - C[0] * b
	B[2] = C[0] * a - A[0] * c
	C[2] = A[0] * b - B[0] * a
	pm = math.sqrt(A[2]*A[2] + B[2]*B[2] + C[2]*C[2])
	try:
		A[2] /= pm
		B[2] /= pm
		C[2] /= pm
	except ZeroDivisionError:
		if (a != 0) and (b != 0):
			a = 0
		elif (b != 0) and (c != 0):
			b = 0
		elif (a != 0) and (c != 0):
			a = 0
		elif a == 0:
			a = 1
		elif b == 0:
			b = 1
		else:
			c = 1
		A[2] = B[0] * c - C[0] * b
		B[2] = C[0] * a - A[0] * c
		C[2] = A[0] * b - B[0] * a
		pm = math.sqrt(A[2]*A[2] + B[2]*B[2] + C[2]*C[2])
		A[2] /= pm
		B[2] /= pm
		C[2] /= pm

	# Vecteur (ij^jk)^ij norme
	A[1] = B[2] * C[0] - C[2] * B[0]
	B[1] = C[2] * A[0] - A[2] * C[0]
	C[1] = A[2] * B[0] - B[2] * A[0]
	pm = math.sqrt(A[1]*A[1] + B[1]*B[1] + C[1]*C[1])
	A[1] /= pm
	B[1] /= pm
	C[1] /= pm
	
		
def setLinkAppCoord(indexLinkApp,data,ori,imolecule, saveindex = None, verbose = 0):
	"""
	setLinkAppCoord: recursively build the coordinates for a linker/appendix
	
	@param indexLinkApp : ith connection in data.
	@param data : the graph connecting rings and linkers.
	@param ori : linker atom instance 
	@param imolecule : the iMolecule to buildd instance
	"""
	isCOUCOU = 0
	if verbose:
		sys.stderr.write("setLinkAppCoord: ori atom type %d\n" % ori.number)
		if ori.symbol == "N":
			sys.stderr.write("COUCOUUOUOUOUOUU\n")
			isCOUCOU = 1
	i = indexLinkApp
	v = [] # liste des voisins de 'ori' ayant des coordonnees
	nb = [] # liste des voisins de 'ori' qui n'ont pas encore de coordonnees
	for oatom in ori.oatoms:
		if oatom.x != None:
			v.append(oatom)
		else:
			nb.append(oatom)
	
	# special case ??
	if (len(data[0]) == 4) and (ori == data[i][2].molecule.atoms[data[i][3]]):
		ori = data[i][0].molecule.atoms[data[i][1]]
		for oatom in data[i][0].molecule.atoms[data[i][1]].oatoms:
			if oatom.x:
				v.append(oatom)
	if (not nb) or (not v):
		return			
	
	# We get the index of the ori atom in the imolecule instance
	indexOri = getIndexAtom(ori,imolecule.molecule)
	nbOfGhostsToPlace = len(imolecule.molecule.atoms[indexOri].oatoms) - len(nb) - len(v)
	if verbose:
		sys.stderr.write("setLinkAppCoord: nbOfGhostsToPlace %d\n" % nbOfGhostsToPlace)

	# C=N- We need to take bond cardinality into account.
	if (ori.number == 7):
		valence = 0
		for oatom in v:
			indexOatom = getIndexAtom(oatom,imolecule.molecule)
			for b in imolecule.molecule.atoms[indexOri].bonds:
				if b.bondtype == 2:
					nbOfGhostsToPlace = 1
			

	if isCOUCOU:
		sys.stderr.write("ori ghosts: %d\n" % nbOfGhostsToPlace)
	#print ori.symbol, imolecule.molecule.atoms[indexOri].index, len(imolecule.molecule.atoms[indexOri].oatoms), len(nb), len(v)
	# Non planar azote
	if (ori.number == 7) and (len(imolecule.molecule.atoms[indexOri].oatoms) == 3) and not (imolecule.molecule.atoms[indexOri].aromatic):
		if isCOUCOU:
			sys.stderr.write("Non planar Azote\n")
		for oat in imolecule.molecule.atoms[indexOri].oatoms:
			if oat.symbol == 'C' and len(oat.oatoms) == 3:
				nbOfGhostsToPlace -= 1
			elif oat.symbol == 'C' and oat.aromatic:
				nbOfGhostsToPlace -= 1
		nbOfGhostsToPlace += 1
	# Sulfur
	if (ori.number == 16) and (len(ori.oatoms) == 3) and not (ori.aromatic):
		nbOfGhostsToPlace += 1
		
	if nbOfGhostsToPlace < 0:
		nbOfGhostsToPlace = 0

	# gestion des positions cis/trans autour des doubles liaisons
	for bond in ori.bonds:
		oatom = bond.xatom(ori)
		if bond.dbo:
			dbo = bond.dbo[0]
			if (len(data[0]) == 4) and (oatom.handle == data[i][0].molecule.atoms[data[i][1]].handle) and (nb):
				oatom = data[i][0].molecule.atoms[data[i][1]]
				bondToPlace = 0
				if dbo[0].atoms[0].handle in [at.handle for at in data[i][0].molecule.atoms]:
					bondToPlace = 1
				if dbo[(bondToPlace+1)%2].atoms[0].handle == data[i][0].molecule.atoms[data[i][1]].handle:
					handlePt1 = dbo[(bondToPlace+1)%2].atoms[1].handle
				else:
					handlePt1 = dbo[(bondToPlace+1)%2].atoms[0].handle
				for at in data[i][0].molecule.atoms:
					if at.handle == handlePt1:
						pt1 = at
				for at in nb:
					if at.handle in [at4.handle for at4 in dbo[bondToPlace].atoms]:
						pt4 = at
				if dbo[2] == 2:
					dbo = ["", "", 2]
				else:
					dbo = ["", "", 1]
			elif (oatom in v) and (nb):
				if dbo[0] in ori.bonds:
					pt1 = dbo[1].xatom(oatom)
					pt4 = dbo[0].xatom(ori)
				else:
					pt4 = dbo[1].xatom(ori)
					pt1 = dbo[0].xatom(oatom)
			else:
				continue
			if (pt1.x != None and pt4.x == None) or (pt1.x == None and pt4.x != None):
				crds = setCoord(ori,v,nb,nbOfGhostsToPlace, saveindex = saveindex, verbose  = verbose)
				try:
				  	index4 = nb.index(pt4)
                                        v12 = [oatom.x - pt1.x, oatom.y - pt1.y, oatom.z - pt1.z]
				except:
					index4 = nb.index(pt1)
					v12 = [oatom.x - pt4.x, oatom.y - pt4.y, oatom.z - pt4.z]
				v34 = [crds[index4][0] - ori.x, crds[index4][1] - ori.y, crds[index4][2] - ori.z]	
    			        scalar = v12[0]*v34[0] + v12[1]*v34[1] + v12[2]*v34[2]
				# print "scalar:", scalar, "dbo:", dbo, nb
				if not (((dbo[2] == 1) and (scalar > 0)) or ((dbo[2] == 2) and (scalar < 0))):
					#if not (((dbo[2] == 1) and (scalar < 0)) or ((dbo[2] == 2) and (scalar > 0))):
					#nb[0], nb[1] = nb[1], nb[0]
					if len(nb) == 2:
						nb[0], nb[1] = nb[1], nb[0]
						crds = setCoord(ori,v,nb, saveindex = saveindex, verbose  = verbose)
					elif len(nb) == 1:
						if isCOUCOU:
							sys.stderr.write("Forcing -1 \n")
	
						# crds = setCoord(ori,v,nb,-1, verbose  = verbose)
						sens = 1
						if dbo[2] == 2:
							sens = -1
						crds = setCoord(ori,v,nb,nbOfGhostsToPlace, sens, saveindex = saveindex, verbose  = verbose)
				break
		elif (len(data[0]) == 4) and (ori == data[i][0].molecule.atoms[data[i][1]]) and (oatom.symbol != 'H') and \
		     (oatom.x  == None) and (imolecule.axialEquatorial):	
			descript = imolecule.axialEquatorial[getIndexAtom(oatom,imolecule.molecule)]
			if descript in ['A','E']:
				if oatom != nb[0]:
					if verbose:
						sys.stderr.write("nb Swap on AE, case 1\n")
					nb[0], nb[1] = nb[1], nb[0]
				crds = setCoord(ori,v,nb,nbOfGhostsToPlace, saveindex = saveindex, verbose  = verbose)
			
				if not isAxEqValidPosition(ori,data[i][0].molecule.atoms,crds,descript):	
					if verbose:
						sys.stderr.write("nb Swap on AE, case 2\n")
					nb[0], nb[1] = nb[1], nb[0]
					crds = setCoord(ori,v,nb,nbOfGhostsToPlace, saveindex = saveindex, verbose  = verbose)
				break
	else:
		# Chiral sense
		sens = 1
		if ori.chiral_class == 2:
			sens = -1
		# Not a good idea to swap here ...
		# if (ori.chiral_class == 2) and (len(nb) > 1):
		# 	if len(v) == 2:
		# 		if getIndexAtom(v[0],imolecule.molecule) > getIndexAtom(v[1],imolecule.molecule):
		# 			#print "getIndexAtom(v[0],imolecule.molecule) < getIndexAtom(v[1],imolecule.molecule)"
		# 			v[0], v[1] = v[1], v[0]
		# 	else:
		# 		if verbose:
		# 			sys.stderr.write("nb swap on chiral_class 2\n")
		# 		nb[0], nb[1] = nb[1], nb[0]
		
		# elif (ori.symbol == 'N') and (len(ori.oatoms) == 3):
		# 	for bond in ori.bonds:
		# 		if (bond.xatom(ori).symbol == "C") and (bond.bondtype == 1) and \
		# 		   (bond.xatom(ori) in v) and (len(bond.xatom(ori).oatoms) == 3) and not (bond.xatom(ori).rings):
		# 			if len(nb) > 1:
		# 				if verbose:
		# 					sys.stderr.write("nb swap on N case\n")
		# 				nb[0], nb[1] = nb[1], nb[0]
		# 			else:
		# 				sens = -1
		# elif (ori.symbol == 'C') and (len(ori.oatoms) == 3):
		# 	for bond in ori.bonds:
		# 		if (bond.xatom(ori).symbol == "N") and (bond.bondtype == 1) and \
		# 		   (bond.xatom(ori) in v) and (len(bond.xatom(ori).oatoms) == 3) and not (bond.xatom(ori).rings):
		# 			if len(nb) > 1:
		# 				if verbose:
		# 					sys.stderr.write("nb Swap on C3 \n")
		# 				nb[0], nb[1] = nb[1], nb[0]
		# 			else:
		# 				sens = -1
		crds = setCoord(ori,v,nb,nbOfGhostsToPlace, sens, saveindex = saveindex, verbose  = verbose) # calcul des coordonnees -> crds
				
	for k in range(len(nb)) :
		# affectation des coordonnees aux atomes
		nb[k].x = crds[k][0]
		nb[k].y = crds[k][1]
		nb[k].z = crds[k][2]
		if verbose: 
			sys.stderr.write("setLinkAppCoord: for %s %d : %f %f %f\n" % (nb[k].symbol, nb[k].index, nb[k].x, nb[k].y, nb[k].z))
					
	for atom in nb:
		setLinkAppCoord(i,data,atom,imolecule, saveindex = saveindex, verbose = verbose)


def dist(ori,neighbour):
	#print "dist----\n", ori, ori.symbol, ori.oatoms, [oat.symbol for oat in ori.oatoms]
	#print neighbour, neighbour.symbol, neighbour.oatoms, [oat.symbol for oat in neighbour.oatoms], "---dist\n"
	try:
		if ori.symbol == "H" and neighbour.symbol in ["C", "N", "O", "S"]:
			return HBONDS[neighbour.symbol + "-" + ori.symbol]
		if neighbour.symbol == "H" and ori.symbol in ["C", "N", "O", "S"]:
			return HBONDS[ori.symbol + "-" + neighbour.symbol]
	except KeyError:
		pass
		
	atOri = getEnghHuberAtomType(ori)
	atNeigh = getEnghHuberAtomType(neighbour)
	if atOri and atNeigh:
		try:
			return BONDS1[atOri + "-" + atNeigh]
		except KeyError:
			try:
				return BONDS1[atNeigh + "-" + atOri]
			except KeyError:
				pass
	if (ori.symbol == "C" and neighbour.symbol in ["Br", "I", "Cl", "F", "P"]) or \
	   (neighbour.symbol == "C" and ori.symbol in ["Br", "I", "Cl", "F", "P"]):
		if ori.findbond(neighbour).bondtype == 1:
			l = "-"
		elif ori.findbond(neighbour).bondtype == 2 or ori.findbond(neighbour).bondorder == 2:
			l = "="
		elif ori.findbond(neighbour).bondtype == 3:
			l = "#"
		try:
			return BONDS2[ori.symbol + l + neighbour.symbol]
		except KeyError:
			try:
				return BONDS2[neighbour.symbol + l + ori.symbol]	
			except KeyError:
				return dist2(ori,neighbour)
	return dist2(ori,neighbour)
"""
def distBBB(ori,neighbour):
	res = ["",""]
	if (ori.number == 1) or (neighbour.number == 1):
		if (ori.symbol in ["C","N","O","S"]) or (neighbour.symbol in ["C","N","O","S"]):
			res[0] = ori.symbol
			res[1] = neighbour.symbol
		else:
			return dist2(ori,neighbour)
	elif ori.number == 6:
		if len(ori.oatoms) == 3:
			res[0] = "CA"
		else :
			res[0] = "C"
		if neighbour.number == 6:
			if len(neighbour.oatoms) == 3:
				res[1] = "CA"
			else :
				res[1] = "C"
		else:
			res[1] = neighbour.symbol
	elif ori.number in [7,8,16]:
		res[0] = ori.symbol
		if (neighbour.number == 6) and (len(neighbour.oatoms) == 3):
			res[1] = "CA"
		else :
			res[1] = neighbour.symbol
	else :
		res[0] = ori.symbol
		if (neighbour.symbol == 6) and (len(neighbour.oatoms) == 3):
			res[1] = "CA"
		else :
			res[1] = neighbour.symbol
	try:
		return BONDS1[res[0] + "-" + res[1]]
	except KeyError:
		try:
			return BONDS1[res[1] + "-" + res[0]]
		except KeyError:
			return dist2(ori,neighbour)
"""
    
def dist2(ori,neighbour):
	try:
		tmp = ori.mmffAtomType
	except:
		ori.mmffAtomType = getMMFFAtomType(ori)		
	try:
		tmp = neighbour.mmffAtomType
	except:	
		neighbour.mmffAtomType = getMMFFAtomType(neighbour)
	atomsTypes = [ori.mmffAtomType, neighbour.mmffAtomType]
	atomsTypes.sort()
	try:
		return BONDS3[str(atomsTypes[0]) + "-" + str(atomsTypes[1])]
	except KeyError:
		atomsTypes = [getSimplestMMFFAtomType(ori), getSimplestMMFFAtomType(neighbour)]
		atomsTypes.sort()
		try:
			return BONDS3[str(atomsTypes[0]) + "-" + str(atomsTypes[1])]
		except:
			#raise MissingParameter, "Distance " + str(atomsTypes[0]) + "-" + str(atomsTypes[1])
			allAtomsTypes = [getAllMMFFAtomsTypesForSymbol(ori), getAllMMFFAtomsTypesForSymbol(neighbour)]
			for atomT1 in allAtomsTypes[0]:
				for atomT2 in allAtomsTypes[1]:
					atomsTypes = [atomT1, atomT2]
					atomsTypes.sort()
					try:
						return BONDS3[str(atomsTypes[0]) + "-" + str(atomsTypes[1])]
					except:
						pass
	

def angle(neighbour1,ori,neighbour2):
	atNeigh1 = getEnghHuberAtomType(neighbour1)
	atOri = getEnghHuberAtomType(ori)
	atNeigh2 = getEnghHuberAtomType(neighbour2)

	if atNeigh1 and atOri and atNeigh2:
		try:
			return (ANGLES1[atNeigh1+"-"+atOri+"-"+atNeigh2]*math.pi)/180.
		except KeyError:
			try:
				return (ANGLES1[atNeigh2+"-"+atOri+"-"+atNeigh1]*math.pi)/180.
			except KeyError:
				return angle2(neighbour1,ori,neighbour2)
	
	return angle2(neighbour1,ori,neighbour2)


def angle2(neighbour1,ori,neighbour2):
	try:
		tmp = ori.mmffAtomType
	except:
		ori.mmffAtomType = getMMFFAtomType(ori)		
	try:
		tmp = neighbour1.mmffAtomType
	except:	
		neighbour1.mmffAtomType = getMMFFAtomType(neighbour1)
	try:
		tmp = neighbour2.mmffAtomType
	except:	
		neighbour2.mmffAtomType = getMMFFAtomType(neighbour2)		

	atomsTypesOri = [neighbour1.mmffAtomType, ori.mmffAtomType, neighbour2.mmffAtomType]
	if atomsTypesOri[0] > atomsTypesOri[2]:
		atomsTypesOri[0], atomsTypesOri[2] = atomsTypesOri[2], atomsTypesOri[0]
	atomsTypes = atomsTypesOri[:]
	for i in range(2):
		for j in range(9):
			try:
				res = ANGLES2[str(j)+"-"+str(atomsTypes[0])+"-"+str(atomsTypes[1])+"-"+str(atomsTypes[2])]
				if res:
					return (res*math.pi)/180.
			except:
				pass
		if i == 0:
			atomsTypes = [0,atomsTypesOri[1],0]


	atomsTypes = [getSimplestMMFFAtomType(neighbour1), getSimplestMMFFAtomType(ori), getSimplestMMFFAtomType(neighbour2)]
	if atomsTypes[0] > atomsTypes[2]:
		atomsTypes[0], atomsTypes[2] = atomsTypes[2], atomsTypes[0]
	atomsTypes2 = atomsTypes[:]
	for i in range(2):
		for j in range(9):
			try:
				#print str(j)+"-"+str(atomsTypes[0])+"-"+str(atomsTypes[1])+"-"+str(atomsTypes[2])
				res = ANGLES2[str(j)+"-"+str(atomsTypes2[0])+"-"+str(atomsTypes2[1])+"-"+str(atomsTypes2[2])]
				if res:
					return (res*math.pi)/180.
			except:
				pass
		if i == 0:
			atomsTypes2 = [0,atomsTypes2[1],0]


	allAtomsTypes = [getAllMMFFAtomsTypesForSymbol(neighbour1), getAllMMFFAtomsTypesForSymbol(ori), getAllMMFFAtomsTypesForSymbol(neighbour2)]
	for atomT1 in allAtomsTypes[0]:
		for atomT2 in allAtomsTypes[1]:
			for atomT3 in allAtomsTypes[2]:
				atomsTypes = [atomT1, atomT2, atomT3]
				if atomsTypes[0] > atomsTypes[2]:
					atomsTypes[0], atomsTypes[2] = atomsTypes[2], atomsTypes[0]
					for i in range(2):
						for j in range(9):
							try:
								res = ANGLES2[str(j)+"-"+str(atomsTypes[0])+"-"+str(atomsTypes[1])+"-"+str(atomsTypes[2])]
								if res:
									return (res*math.pi)/180.
						       	except:
								pass
					       	if i == 0:
							atomsTypes = [0,atomsTypes[1],0]
					
	
	raise MissingParameter, "Angle "+str(atomsTypesOri[0])+"-"+str(atomsTypesOri[1])+"-"+str(atomsTypesOri[2])


def getDihedralAngleValues(atom1,atom2,atom3,atom4, verbose = 0):
	try:
		tmp = atom1.mmffAtomType
	except:
		atom1.mmffAtomType = getMMFFAtomType(atom1)
	try:
		tmp = atom2.mmffAtomType
	except:
		atom2.mmffAtomType = getMMFFAtomType(atom2)
	try:
		tmp = atom3.mmffAtomType
	except:
		atom3.mmffAtomType = getMMFFAtomType(atom3)
	try:
		tmp = atom4.mmffAtomType
	except:
		atom4.mmffAtomType = getMMFFAtomType(atom4)	

	# MMFF atom types Not enough. MUST LOOK AT BOND VALENCE.
	for bond in atom2.bonds:
		if bond.xatom(atom2) == atom3:
			bondType = bond.bondtype
			if verbose:
				sys.stderr.write("Bond type : %d\n" % bondType)
			break
	# if (bondType == 1) and (atom2.mmffAtomType == 2):
	if (bondType == 1) and (atom2.mmffAtomType in [2, 37, 78]):
		atom2.mmffAtomType = 1
	# if (bondType == 1) and (atom3.mmffAtomType == 2):
	if (bondType == 1) and (atom3.mmffAtomType in [2, 37, 78]):
		atom3.mmffAtomType = 1
	
	if verbose:
		sys.stderr.write("%s-%d %s-%d %s-%d %s-%d : starting with %d %d %d %d\n" % (atom1.symbol, atom1.index, atom2.symbol, atom2.index, atom3.symbol, atom3.index ,atom4.symbol, atom4.index,  atom1.mmffAtomType, atom2.mmffAtomType, atom3.mmffAtomType, atom4.mmffAtomType))
	atomsTypesOri = [atom1.mmffAtomType, atom2.mmffAtomType, atom3.mmffAtomType, atom4.mmffAtomType]
	atomsTypes = atomsTypesOri[:]
	for i in range(8):
		for j in [0,1,2,4,5]:
			try:
				if verbose:
					sys.stderr.write("%d %d %d %d\n" % (atomsTypes[0], atomsTypes[1], atomsTypes[2], atomsTypes[3]))
				guess = str(j) + "-" + str(atomsTypes[0]) + "-" + str(atomsTypes[1]) + "-" + \
				                             str(atomsTypes[2]) + "-" + str(atomsTypes[3])
				res = DIEDRES[guess]
				if(res):
					if len(res) == 3 and (atom2.rings or atom3.rings):
						rres = list(res)
						# for pos in range(len(rres)):
						# 	aVal = rres[pos][0]+180.
						# 	if aVal > 180.:
						# 		aVal = aVal - 360.
						# 	rres.append((aVal,rres[pos][1]))
						rres.append((90.,0.0))
						rres.append((-90.,0.0))
						# rres.append((0.,0.0))
						res = tuple(rres)
							
					if verbose:
						sys.stderr.write("Found type %s for i %d j %d : %s\n" % (str(res), i, j, guess))
					
					return res
			except:
				pass
		if i%2 == 0:
			atomsTypes.reverse()
		elif i == 1:
			atomsTypes = atomsTypesOri[:len(atomsTypesOri)-1]
			atomsTypes.extend([0])
		elif i == 3:
			atomsTypes = [0]
			atomsTypes.extend(atomsTypesOri[1:])
		elif i == 5:
			atomsTypes = [0,atomsTypesOri[1],atomsTypesOri[2],0]	


	atomsTypesOri = [getSimplestMMFFAtomType(atom1), getSimplestMMFFAtomType(atom2), getSimplestMMFFAtomType(atom3), getSimplestMMFFAtomType(atom4)]	
	atomsTypes = atomsTypesOri[:]
	if verbose:
		sys.stderr.write("%s-%d %s-%d %s-%d %s-%d : trying with %d %d %d %d\n" % (atom1.symbol, atom1.index, atom2.symbol, atom2.index, atom3.symbol, atom3.index ,atom4.symbol, atom4.index,  atomsTypes[0], atomsTypes[1], atomsTypes[2], atomsTypes[3] ))
	for i in range(8):
		for j in [0,1,2,4,5]:
			try:
				res = DIEDRES[str(j) + "-" + str(atomsTypes[0]) + "-" + str(atomsTypes[1]) + "-" + \
			                             str(atomsTypes[2]) + "-" + str(atomsTypes[3])]
				if(res):
					return res
			except:
				pass
		if i%2 == 0:
			atomsTypes.reverse()
		elif i == 1:
			atomsTypes = atomsTypesOri[:len(atomsTypesOri)-1]
			atomsTypes.extend([0])
		elif i == 3:
			atomsTypes = [0]
			atomsTypes.extend(atomsTypesOri[1:])
		elif i == 5:
			atomsTypes = [0,atomsTypesOri[1],atomsTypesOri[2],0]

	
        allAtomsTypes = [getAllMMFFAtomsTypesForSymbol(atom1), getAllMMFFAtomsTypesForSymbol(atom2), getAllMMFFAtomsTypesForSymbol(atom3), getAllMMFFAtomsTypesForSymbol(atom4)]
	if verbose:
		sys.stderr.write("%s-%d %s-%d %s-%d %s-%d : trying with even more generic\n" % (atom1.symbol, atom1.index, atom2.symbol, atom2.index, atom3.symbol, atom3.index, atom4.symbol, atom4.index ))

        for atomT1 in allAtomsTypes[0]:
                for atomT2 in allAtomsTypes[1]:
                        for atomT3 in allAtomsTypes[2]:
                                for atomT4 in allAtomsTypes[3]:
                                        atomsTypesOri = [atomT1, atomT2, atomT3, atomT4]
                                        atomsTypes = atomsTypesOri[:]
	                                for i in range(8):
                                                for j in [0,1,2,4,5]:
                                                        try:
                                                                res = DIEDRES[str(j) + "-" + str(atomsTypes[0]) + "-" + str(atomsTypes[1]) + "-" + \
			                                                      str(atomsTypes[2]) + "-" + str(atomsTypes[3])]
                                                                if(res):
                                                                        return res
                                                        except:
                                                                pass
                                                if i%2 == 0:
                                                        atomsTypes.reverse()
                                                elif i == 1:
                                                        atomsTypes = atomsTypesOri[:len(atomsTypesOri)-1]
                                                        atomsTypes.extend([0])
                                                elif i == 3:
                                                        atomsTypes = [0]
                                                        atomsTypes.extend(atomsTypesOri[1:])
                                                elif i == 5:
                                                        atomsTypes = [0,atomsTypesOri[1],atomsTypesOri[2],0]
                        
	raise MissingParameter, "Dihedral "+ str(atomsTypesOri[0]) + "-" + str(atomsTypesOri[1]) + "-" + \
	                                       str(atomsTypesOri[2]) + "-" + str(atomsTypesOri[3])

def getMMFFPartialAtomicCharge(atom):
	fAc = atom.charge #formal atomic charge, a ameliorer afin de tenir compte 
	                  #des charges comme pour l'azote dans guanidium de charge 1/3
        try:
		atomType = atom.mmffAtomType
	except:
		atom.mmffAtomType = getMMFFAtomType(atom)
		atomType = atom.mmffAtomType

	neighboursAtomsTypes = []
	neighboursSimplestAtomsTypes = []
	for oatom in atom.oatoms:
		try:
			neighboursAtomsTypes.append(oatom.mmffAtomType)
		except:
			oatom.mmffAtomType = getMMFFAtomType(oatom)
			neighboursAtomsTypes.append(oatom.mmffAtomType)
			neighboursSimplestAtomsTypes.append(getSimplestMMFFAtomType(oatom))
	
	bondChargeIncrements = 0.
	for i, neighbourAtomType in enumerate(neighboursAtomsTypes):
		excepted = False
		if neighbourAtomType > atomType:
			try:
				bondChargeIncrements += -BCI[str(atomType)+"-"+str(neighbourAtomType)]
			except:
				excepted = True;
		else:
			try:
				bondChargeIncrements += BCI[str(neighbourAtomType)+"-"+str(atomType)]
			except:
				excepted = True;
		if excepted:
			try:
				pbci1 = PBCI[str(atomType)][0]
			except:
				try:
					pbci1 = PBCI[getSimplestMMFFAtomType(atom)][0]
				except:
					pbci1 = 0
			try:
				pbci2 = PBCI[str(neighbourAtomType)][0]
			except:
				try:
					pbci2 = PBCI[neighboursSimplestAtomsTypes[i]][0]
				except:
					pbci2 = 0
			bondChargeIncrements += pbci1 - pbci2
	return fAc + bondChargeIncrements

def getNextAtoms(atomPrec, atom):
	res = []
	for at in atom.oatoms:
		if at != atomPrec:
			res.append(at)
	return res
	
def getMol2AtomType(atom):
	"""get the type of an atom used in Mol2"""
	if atom.symbol in ['F','Cl','Br','I','Cu','Zn','Se','Mo','Sn','Fe','Mn','Ca','K','Si','Al','Mg','Na','Li']:
		return atom.symbol
	if atom.symbol == 'Cr':
		if len(atom.oatoms) == 4:
			return 'Cr.th'
		if len(atom.oatoms) == 8:
			return 'Cr.oh'
	if (atom.symbol == 'P') and (len(atom.oatoms) == 4):
		return 'P.3'
	if (atom.symbol == 'Co') and (len(atom.oatoms) == 8):
		return 'Co.oh'
	if atom.symbol == 'C':
		if atom.aromatic:
			return 'C.ar'
		if atom.mmffAtomType == 57:
			return 'C.cat'
		if len(atom.oatoms) == 4:
			return 'C.3'
		if len(atom.oatoms) == 3:
			return 'C.2'
		if len(atom.oatoms) == 2:
			return 'C.1'
	if atom.symbol == 'N':
		if atom.aromatic:
			return 'N.ar'
		if (atom.charge == 1) or (len(atom.oatoms) == 4):
			return 'N.4'
		nbDouble = 0
		for bond in atom.bonds:
			if (bond.bondtype == 1) and (bond.xatom(atom).symbol == 'C'):
				for bbond in bond.xatom(atom).bonds:
					if (bbond.bondtype == 2) and (bbond.xatom(bond.xatom(atom)).symbol in ['O','S']):
						return 'N.am'
			if bond.bondtype == 2:
				nbDouble += 1
		if len(atom.oatoms) == 3:
			nbN = simpleO = 0
			for bond in atom.bonds:
				if (bond.bondtype == 1) and (bond.xatom(atom).symbol in ['S','P']):
					if len(bond.xatom(atom).oatoms) == 4:
						doubleO = 0
						for bbond in bond.xatom(atom).bonds:
							if (bbond.bondtype == 2) and (bbond.xatom(bond.xatom(atom)).symbol == 'O'):
								doubleO += 1
						if doubleO == 2:
							#NITROGEN IN SULFONAMIDES or PHOSPHONAMIDES
							for bbbond in atom.bonds:
								if bbbond.xatom(atom).aromatic:
									return 'N.pl3'
							return 'N.3'	
		if len(atom.oatoms) >= 3:
			if nbDouble == 1:
				return 'N.pl3'
			if atom.charge == 1:
				return 'N.3'
			if getNbOfHInOatoms(atom) >= 2:
				return 'N.3'
			return 'N.pl3'
		if len(atom.oatoms) == 2:
			return 'N.2'
		return 'N.1'
	if atom.symbol == 'O':	
		if atom.aromatic:
			return 'O.2'
		if (len(atom.oatoms) == 1) and (atom.bonds[0].xatom(atom).symbol in ['C','S','P']) and \
		   (atom.bonds[0].bondtype <= 2):
			doubleO = simpleO = 0
			for bond in atom.bonds[0].xatom(atom).bonds:
				if bond.xatom(atom.bonds[0].xatom(atom)).symbol == 'O':
					if bond.bondtype == 2:
						doubleO += 1
					if (bond.bondtype == 1) and ((getNbOfHInOatoms(bond.xatom(atom.bonds[0].xatom(atom))) == 1) or \
												 (bond.xatom(atom.bonds[0].xatom(atom)).charge == -1)):
						simpleO += 1
			if (doubleO >= 1) and (simpleO >= 1):
				#if ((simpleO == 3) and (len(atom.oatoms) == 4)) or ((simpleO >= 1) and (len(atom.oatoms) == 3)):
				return 'O.co2'
		if len(atom.oatoms) == 2:
			for bond in atom.bonds:
				if bond.xatom(atom).symbol == 'N':
					for bbond in bond.xatom(atom).bonds:
						if (len(bond.xatom(atom).oatoms) == 3) and (bbond.bondtype == 2):
							return 'O.2'
			return 'O.3'
		#if atom.charge == -1:
		#	return '0.3'
		return 'O.2'
	if atom.symbol == 'S':
		if atom.aromatic:
			return 'S.2'
		doubleO = simpleO = 0
		for bond in atom.bonds:
			if (bond.bondtype == 2) and (bond.xatom(atom).symbol == 'O'):
				doubleO += 1
			if (bond.bondtype == 1) and (bond.xatom(atom).symbol == 'O'):
				simpleO += 1
		if (doubleO == 2) and (simpleO == 0):
			return 'S.o2'
		if doubleO == 1:
			return 'S.O'
		if len(atom.oatoms) >= 2:
			return 'S.3'
		return 'S.2'
	if atom.symbol == 'H':
		return 'H'
	return 'Hev'
	
MMFFAtomTypes = {
	"C" : [1, 2, 3, 4, 20, 22, 30, 37, 41, 57,  60, 63, 64, 78, 80],
	"O" : [6, 7, 32, 35, 49, 51, 59, 70],
	"N" : [8, 9, 10, 34, 38, 39, 40, 42, 43, 45, 46, 47, 48, 53, 54, 55, 56, 58, 61, 62, 65, 66, 67, 68, 69, 76, 79, 81, 82],
	"S" : [15, 16, 17, 18, 44, 72, 73, 74],
	"F" : [11, 87, 88, 89],
	"CL": [12, 77, 90],
	"BR": [13, 91],
	"I" : [14],
	"SI": [19],
	"P" : [25, 26, 75],
	"LI": [92],
	"NA": [93],
	"K" : [94],
	"ZN": [95],
	"CA": [96],
	"CU": [97, 98],
	"MG": [99],
}
	
def getCarbonMMFFAtomType(atom):
	"""
	possible types are:
	1: alkyl
	2: sp2
	3: hetero sp2
	4: sp
	20: aliphatic in 4-ring
	22: aliphatic in 3-ring
	30: olefinic in 4 ring (e.g. C1-C=C-C1)
	37: aromatic
	41: carboxylate anion
	57: guanidium (charged +)
	60: isonitrile
	63, 64: not managed (5 ring alpha, resp. beta to O,N or S)
	78: hetero aromatic 5 ring
	80: imidazole
	"C" : [1, 2, 3, 4, 20, 22, 30, 37, 41, 57,  60, 63, 64, 78, 80],	
	"""
	if not atom.rings:
		if len(atom.oatoms) == 4:
			# ALKYL CARBON, SP3
			return 1
		if len(atom.oatoms) == 3:
			o1ats = []
			o2ats = []
			for bond in atom.bonds:
				if bond.bondtype in [2,4]:
					if bond.xatom(atom).symbol not in o2ats:
						o2ats.append(bond.xatom(atom).symbol)
				elif bond.bondtype in [1]:
					if bond.xatom(atom).symbol not in o1ats:
						o1ats.append(bond.xatom(atom).symbol)
				if (bond.xatom(atom).charge == 1) and (len(o2ats) == 1) and (o2ats[0] == "N"):
					# GUANIDINIUM CARBON or C IN +N=C-N RESONANCE STRUCTURES
					return 57
				if (bond.xatom(atom).charge == -1) and (("O" in o2ats) or ("S" in o2ats)):
					# CARBOXYLATE ANION CARBON or CARBON IN THIOCARBOXYLATE ANION
					return 41
						# CARBOXYL, THIAMIDE, CARBONYL, CARBONYL-LIKE ETC CARBON
				if len(o2ats) and (("O" in o2ats) or ("N" in o2ats) or ("S" in o2ats) or ("P" in o2ats)):
					return 3
			return 2
		# if (atom.charge == 0) and (len(atom.oatoms) == 2):
		if (len(atom.oatoms) == 2):
			# ACETYLENIC or ALLENIC CARBON
			return 4
		if (len(atom.oatoms) == 1) and (atom.bonds[0].xatom(atom).symbol == 'N'):
			# ISONITRILE CARBON
			return 60
		
# 	"C" : [30],
	else: # in ring
		if atom.aromatic:
			for ring in atom.rings:
				if len(ring.atoms) == 5:
					for oatom in ring.atoms:
						if oatom.symbol != "C":
							# GENERAL CARBON IN 5-MEMBERED HETEROAROMATIC RING
							return 78
			# GENERIC AROMATIC CARBON IN RING AS IN BENZENE, PYRROLE
			# We do not manage 63 and 64 types
			return 37
		if len(atom.oatoms) == 4:
			for ring in atom.rings:
				if len(ring.atoms) == 3:
					return 22
				if len(ring.atoms) == 4:
					for oatom in ring.atoms:
						if (oatom != atom) and (oatom.symbol not in  ["C", "H"]):
							return 20
			nbN = 0
			for bond in atom.bonds:
				if bond.xatom(atom).symbol == 'N':
					nbN += 1
			if nbN == 2:
				#C IN N-C-N IN IMIDAZOLIUM ION
				return 80
			return 1
		elif len(atom.oatoms) == 3:
			for ring in atom.rings:
				if len(ring.atoms) == 4:
					return 30
			for bond in atom.bonds:
				if (bond.xatom(atom).symbol not in  ["C", "H"]):
					return 3
			return 2
		elif len(atom.oatoms) == 2:
			# We should not get there but ...
			return 3

	return None

def getOxygenMMFFAtomType(atom):
	"""
	Here, we are interested in some cases only.
	-O- : 6
	O=x : 7

	Possible types:
	6: divalent oxygen as in x-O-x
	7: monovalent oxygen as in O=x
	37: O=P
	"O" : [6, 7, 32, 35, 49, 51, 59, 70],

	"""
	if len(atom.oatoms) == 2:
		return 6
	if len(atom.oatoms) == 1:
		for bond in atom.bonds:
			if (bond.xatom(atom).symbol == "P"):
				return 32
		return 7
	return None

def getSulfurMMFFAtomType(atom):
	"""
	"S" : [15, 16, 17, 18, 44, 72, 73, 74],
	Here, we are interested in some cases only.
	-S- : 15
	S-x : 72

	Possible types:
	15: divalent sulfur as in x-S-x
	16: monovalent as in S=C-
	17: -S=O or -S=N
	72: terminal sulfur 
	74: sulfinyl =S=O
	"S" : [15, 16, 17, 18, 44, 72, 73, 74],

	"""
	if len(atom.oatoms) == 4:
		return 18

	if len(atom.oatoms) == 3:
		return 17
		# for bond in atom.bonds:
		# 	if bond.bondtype == 2:
		# 		if (bond.xatom(atom).symbol in ["O", "N"]):
		# 			return 17
		
	if len(atom.oatoms) == 2:
		rs = 15
		for bond in atom.bonds:
			if bond.bondtype == 2:
				if (bond.xatom(atom).symbol == "O"):
					rs = 74
				rs = 15
		return 15

	if len(atom.oatoms) == 1:
		for bond in atom.bonds:
			if bond.bondtype == 2:
				if (bond.xatom(atom).symbol == "C"):
					return 16
			elif (bond.xatom(atom).symbol == "P"):
				return 72
		return 16
	return None

def betaAtoms(atom):
	"""
	List atoms in beta position
	"""
	rs = []
	for oatom in atom.oatoms:
		for ooatom in oatom.oatoms:
			if ooatom != atom and ooatom not in rs:
				rs.append(ooatom)
	return rs

def isAmideLike(atom, verbose = 0):
	"""
	is amide or thioamide (nitrogen)
	"""
	if verbose:
		sys.stderr.write("isAmideLike:\n")
	for bond in atom.bonds:

		if verbose:
			sys.stderr.write("%d %s: " % ( atom.index, atom.symbol) )
		if (bond.xatom(atom).symbol == "C") and (bond.bondtype == 1):
			if verbose:
				sys.stderr.write("%d %s " % ( bond.xatom(atom).index, bond.xatom(atom).symbol) )
			oatom = bond.xatom(atom)
			for obond in oatom.bonds:
				if (obond.bondtype == 2 ) and (obond.xatom(oatom).symbol in ["O","S"]):
					if verbose:
						sys.stderr.write("isAmideLike: done ...\n")
					return True
		if (bond.xatom(atom).symbol == "N") and (bond.bondtype == 1):
			oatom = bond.xatom(atom)
			for obond in oatom.bonds:
				if (obond.bondtype == 2 ) and (obond.xatom(oatom).symbol in ["C","N"]):
					if verbose:
						sys.stderr.write("isAmideLike: done ...\n")
					return True
	if verbose:
		sys.stderr.write("isAmideLike: done ...\n")
	return False
		
def isEnamineLike(atom):
	"""
	is enamine like (nitrogen)
	-N-C=C, -N-C=N, -N-C=P, -N-C=-C
	"""
	for bonds in atom.bonds:
		if (bond.xatom(atom).symbol == "C") and (bond.bondtype == 1):
			oatom = bond.xatom(atom)
			for obond in oatom.bonds:
				if (obond.bondtype == 2 ) and (obond.xatom(oatom).symbol in ["C","N"]):
					return True
				if (obond.bondtype == 3 ) and (obond.xatom(oatom).symbol in ["C"]):
					return True
	return False
		
def getNitrogenMMFFAtomType(atom):
	"""
	"N" : [8, 9, 10, 34, 38, 39, 40, 42, 43, 45, 46, 47, 48, 53, 54, 55, 56, 58, 61, 62, 65, 66, 67, 68, 69, 76, 79, 81, 82]

	Possible types:
	8 : amine (3 voisins, simple liaisons)
	9 : imine (C=N-), azo (N=N-)
	10: amide like (-N-C=[O,S], -N-N=[C,N])
	34: quaternary N
	38: aromatic
	39: aromatic 5-ring
	40: enamine like (-N-C=[C,N,P])
	42: triply bonded
	47: azido
	53:  C=N=N or N=N=N
	58: arromatic in pyridinium
	81: positive N in 5-ring
	82: neutral in 5-ring
	"""
	o1ats = []
	o2ats = []
	o3ats = []
	for bond in atom.bonds:
		if bond.bondtype in [2,4]:
			o2ats.append(bond.xatom(atom).symbol)
		elif bond.bondtype in [1]:
			o1ats.append(bond.xatom(atom).symbol)
		elif bond.bondtype in [3]:
			o3ats.append(bond.xatom(atom).symbol)
	betas = betaAtoms(atom)
	
	if len(atom.oatoms) == 1:
		if len(o3ats) and (o3ats[0].symbol == "N"):
			return 47
		# TRIPLY BONDED
		return 42
	if (len(atom.oatoms) == 4) and (atom.charge == 1):
		# quaternary nitrogen sp3
		return 34
	if len(atom.oatoms) == 3:
		if (atom.rings):
			if (len(atom.rings) == 1) and (len(atom.rings[0].atoms) == 5):
				if (atom.charge == 1): 
					return 81
				else:
					return 82
			if (len(atom.rings) == 1) and (atom.aromatic) and (len(atom.rings[0].atoms) == 6) \
				    and ((len(atom.oatoms) == 3) or (atom.charge == 1)):
				# PYRIDINIUM-TYPE NITROGEN - FORMAL CHARGE=1
				return 58
			
		if isAmideLike(atom):
			return 10
		if isEnamineLike(atom):
			return 40
		# standard sp3 case
		return 8
	if len(atom.oatoms) == 2:
		if not atom.rings:
			# IMINE OR AZO
			if ("C" in o2ats) or ("N" in o2ats):
				return 9
		if atom.aromatic:
			if (atom.rings) and (len(atom.rings) == 1):
				if len(atom.rings[0].atoms) == 5:
				        # NITROGEN, AS IN PYRROLE
					# Note: we do not consider type 79
					return 39
			# Default aromatic
			return 38
		if (len(o2ats) == 2) and (o2ats[0].symbol == "N") or (o2ats[1].symbol == "N"):
			return 53
		if len(o3ats):
			return 42
		return 9	


	return None


def getMMFFAtomType_new(atom):
	rs = None
	if atom.symbol == "C":
		rs = getCarbonMMFFAtomType(atom)
	if atom.symbol == "O":
		rs = getOxygenMMFFAtomType(atom)
	if atom.symbol == "S":
		rs = getSulfurMMFFAtomType(atom)
	if atom.symbol == "N":
		rs = getNitrogenMMFFAtomType(atom)

	if atom.symbol == 'F':
		if atom.charge == -1:
			#FLUORIDE ANION
			rs = 89
		#FLUORINE
		rs = 11
	
	if atom.symbol == 'Cl':
		if atom.charge == -1:
			if len(atom.oatoms) == 3:
				for bond in atom.bonds:
					if not ((bond.xatom(atom).symbol == 'O') and (bond.bondtype == 2)):
						break
				else:
					#CHLORINE IN PERCHLORATE ANION, CLO4(-)
					rs = 77
			#CHLORIDE ANION
			rs = 90
		#CHLORINE
		rs = 12
			
	if atom.symbol == 'Br':
		if atom.charge == -1:
			#BROMIDE ANION
			rs = 91
		#BROMINE
		rs = 13
	
	if atom.symbol == 'I':
		#IODINE
		rs = 14
	
	if atom.symbol == 'Si':
		#SILICONE
		rs = 19
	
	if atom.symbol == 'Fe':
		if atom.charge == 2:
			#IRON +2 CATION
			rs = 87	
		if atom.charge == 3:
			#IRON +3 CATION
			rs = 88
		
	if (atom.symbol == 'Li') and (atom.charge == 1):
		#LITHIUM CATION
		rs = 92

	if (atom.symbol == 'Na') and (atom.charge == 1):
		#SODIUM CATION
		rs = 93
	
	if (atom.symbol == 'K') and (atom.charge == 1):
		#POTASSIUM CATION
		rs = 94
	
	if (atom.symbol == 'Zn') and (atom.charge == 2):
		#PDIPOSITIVE ZINC
		rs = 95
	
	if (atom.symbol == 'Ca') and (atom.charge == 2):
		#DIPOSITIVE CALCIUM
		rs = 96
	
	if atom.symbol == 'Cu':
		if atom.charge == 1:
			#MONOPOSITIVE COPPER
			rs = 97	
		if atom.charge == 2:
			#DIPOSITIVE COPPER
			rs = 98
	
	if (atom.symbol == 'Mg') and (atom.charge == 2):
		#DIPOSITIVE MAGNESIUM CATION
		rs = 99
	
	if atom.symbol == 'P':
		if len(atom.oatoms) == 4:
			#GENERAL TETRACOORDINATE PHOSPHORUS
			rs = 25
		if len(atom.oatoms) == 3:
			#TRICOORDINATE P, AS IN PHOSPHINES
			rs = 26
		if len(atom.oatoms) == 2:
			for bond in atom.bonds:
				if ((bond.xatom(atom).symbol == 'C') and (bond.bondtype == 2)):
					#PHOSPHOROUS DOUBLY BONDED TO CARBON
					rs = 75
	if not rs:
		rs = getSimplestMMFFAtomType(atom)
	return rs
		
def getMMFFAtomType(atom):
	"""get the type of an atom used in MMFF"""
	atomTypes = []
	if atom.symbol == 'C':
		rs = getCarbonMMFFAtomType(atom)
		if not rs:
			rs = getSimplestMMFFAtomType(atom)
		return rs
	if atom.symbol == 'O':
		rs = getOxygenMMFFAtomType(atom)
		if not rs:
			rs = getSimplestMMFFAtomType(atom)
		return rs
	if atom.symbol == 'F':
		if atom.charge == -1:
			#FLUORIDE ANION
			return 89
		#FLUORINE
		return 11
	
	if atom.symbol == 'Cl':
		if atom.charge == -1:
			if len(atom.oatoms) == 3:
				for bond in atom.bonds:
					if not ((bond.xatom(atom).symbol == 'O') and (bond.bondtype == 2)):
						break
				else:
					#CHLORINE IN PERCHLORATE ANION, CLO4(-)
					return 77
			#CHLORIDE ANION
			return 90
		#CHLORINE
		return 12
			
	if atom.symbol == 'Br':
		if atom.charge == -1:
			#BROMIDE ANION
			return 91
		#BROMINE
		return 13
	
	if atom.symbol == 'I':
		#IODINE
		return 14
	
	if atom.symbol == 'Si':
		#SILICONE
		return 19
	
	if atom.symbol == 'Fe':
		if atom.charge == 2:
			#IRON +2 CATION
			return 87	
		if atom.charge == 3:
			#IRON +3 CATION
			return 88
		
	if (atom.symbol == 'Li') and (atom.charge == 1):
		#LITHIUM CATION
		return 92

	if (atom.symbol == 'Na') and (atom.charge == 1):
		#SODIUM CATION
		return 93
	
	if (atom.symbol == 'K') and (atom.charge == 1):
		#POTASSIUM CATION
		return 94
	
	if (atom.symbol == 'Zn') and (atom.charge == 2):
		#PDIPOSITIVE ZINC
		return 95
	
	if (atom.symbol == 'Ca') and (atom.charge == 2):
		#DIPOSITIVE CALCIUM
		return 96
	
	if atom.symbol == 'Cu':
		if atom.charge == 1:
			#MONOPOSITIVE COPPER
			return 97	
		if atom.charge == 2:
			#DIPOSITIVE COPPER
			return 98
	
	if (atom.symbol == 'Mg') and (atom.charge == 2):
		#DIPOSITIVE MAGNESIUM CATION
		return 99
	
	if atom.symbol == 'P':
		if len(atom.oatoms) == 4:
			#GENERAL TETRACOORDINATE PHOSPHORUS
			return 25
		if len(atom.oatoms) == 3:
			#TRICOORDINATE P, AS IN PHOSPHINES
			return 26
		if len(atom.oatoms) == 2:
			for bond in atom.bonds:
				if ((bond.xatom(atom).symbol == 'C') and (bond.bondtype == 2)):
					#PHOSPHOROUS DOUBLY BONDED TO CARBON
					return 75
				
	if atom.symbol == 'S':
		if (atom.rings) and (len(atom.rings) == 1) and (atom.aromatic) and (len(atom.rings[0].atoms) == 5):
			#SULFUR AS IN THIOPHENE
			return 44
		if len(atom.oatoms) == 1:
			if atom.charge == -1:
				#TERMINAL SULFUR - FORMAL CHARGE=-1
				return 72
			if (atom.bonds[0].bondtype == 2) and (atom.bonds[0].xatom(atom).symbol == 'C'):
				#TERMINAL SULFUR DOUBLY BONDED TO CARBON
				return 16	
		if len(atom.oatoms) == 2:
			#print "2"
			simple = simpleP = doubleC = doubleO = simpleOMoins = 0
			for bond in atom.bonds:
				if bond.bondtype == 1:
					if bond.xatom(atom).symbol == 'P':
						simpleP += 1
					elif (bond.xatom(atom).symbol == 'O') and (bond.xatom(atom).charge == -1):
						simpleOMoins += 1
					else:
						simple += 1
				elif bond.bondtype == 2:
					if bond.xatom(atom).symbol == 'O':
						doubleO += 1
					elif bond.xatom(atom).symbol == 'C':
						doubleC += 1
			if simple == 2:
				#SULFUR IN THIOETHERS AND MERCAPTANS
				return 15
			if (simple == 1) and (simpleP == 1):
				#TERMINAL SULFUR BONDED TO PHOSPHORUS
				return 72
			#if (simple == 1) and (doubleC == 1):
			#	#TERMINAL SULFUR DOUBLY BONDED TO CARBON
			#	return 16
			if (simpleOMoins == 1) and (doubleO == 1):
				#SULFUR IN NEGATIVELY CHARGED SULFINATE GROUP
				return 73
			if (doubleC == 1) and (doubleO == 1):
				#SULFINYL SULFUR, EG. IN C=S=O
				return 74
			if doubleC == 1:
				#TERMINAL SULFUR DOUBLY BONDED TO CARBON
				return 16
		if len(atom.oatoms) == 3:
			simple = simpleS = doubleC = doubleO = doubleN = simpleOMoins = 0
			for bond in atom.bonds:
				if bond.bondtype == 1:
					if bond.xatom(atom).symbol == 'S':
						simpleS += 1
					elif (bond.xatom(atom).symbol == 'O') and (bond.xatom(atom).charge == -1):
						simpleOMoins += 1
					else:
						simple += 1
				elif bond.bondtype == 2:
					if bond.xatom(atom).symbol == 'O':
						doubleO += 1	
					elif bond.xatom(atom).symbol == 'N':
						doubleN += 1
					elif bond.xatom(atom).symbol == 'C':
						doubleC += 1		
			if (doubleC == 1) and (doubleO + simpleOMoins == 2):
				#SULFONE SULPHER DOUBLY BONDED TO CARBON
				return 18
			if (simple == 2) and ((doubleN == 1) or (doubleO == 1)):
				#SULFUR, TRICOORD, DOUBLY BONDED TO N, or in SULFOXIDES
				return 17
			if (simpleS == 1) and (simpleOMoins == 1) and (doubleO == 1):
				#TRICOORD SULFUR IN THIOSULFINATE GROUP
				return 73
		if len(atom.oatoms) == 4:
			simple = simpleO = simpleN = doubleO = doubleN = 0
			for bond in atom.bonds:
				if bond.bondtype == 1:
					if bond.xatom(atom).symbol == 'O':
						simpleO += 1
					elif bond.xatom(atom).symbol == 'N':
						simpleN += 1
					else:
						simple += 1
				elif bond.bondtype == 2:
					if bond.xatom(atom).symbol == 'O':
						doubleO += 1	
					elif bond.xatom(atom).symbol == 'N':
						doubleN += 1
			if (simple == 2) and (doubleO + doubleN == 2):
				return 18
			if (simpleO + doubleO >= 2) and (simpleN >= 1):
				return 18
			if simpleO + doubleO >= 2:
				return 18
	if atom.symbol == 'H':
		if (atom.bonds[0].xatom(atom).symbol == 'C') or (atom.bonds[0].xatom(atom).symbol == 'Si'):
			#H  ATTACHED TO C or SI
			return 5
 		if atom.bonds[0].xatom(atom).symbol == 'O':
			if atom.bonds[0].xatom(atom).charge == 1:
				if len(atom.bonds[0].xatom(atom).oatoms) == 2:
					for bond in atom.bonds[0].xatom(atom).bonds:
						if bond.bondtype == 2:
							#HYDROGEN ON OXENIUM OXYGEN
							return 52
				#HYDROGEN ON O+ OXYGEN
				return 50
			if atom.bonds[0].xatom(atom).charge == -1:
				#HYDROGEN IN HYDROXIDE ANION
				return 21
			for oatom in atom.oatoms[0].oatoms:
				if oatom != atom:
					break
			if oatom.symbol == 'H':
				#HYDROGEN IN H2O
				return 31
			if oatom.symbol == 'S':
				#H ON OXYGEN ATTACHED TO SULFUR
				return 33
			if oatom.symbol == 'P':
				#HYDROGEN ON OXYGEN ATTACHED TO PHOSPHOROUS
				return 24
			if oatom.symbol == 'C':
				#if oatom.imp_hcount + oatom.explicit_hcount <= 2: #atom.bonds[0].xatom(atom) >= 2:
				#	#HYDROGEN IN ALCOHOLS
				#	return 21
				if oatom.rings:
					#H-O IN ENOLS AND PHENOLS
					return 29
				for bond in oatom.bonds:
					if bond.bondtype == 2:
						break
				if bond.xatom(oatom).symbol == 'O':
					#H-O IN CARBOXYLIC ACIDS
					return 24
				if bond.xatom(oatom).symbol == 'N':
					#H-O IN HO-C=N
					return 29
			#GENERAL H ON OXYGEN
			return 21
		if atom.bonds[0].xatom(atom).symbol == 'N':
			nType = getMMFFAtomType(atom.bonds[0].xatom(atom))
			if nType in [10, 40, 43]:
				return 28
			elif nType == 9:
				return 27
			elif len(atom.bonds[0].xatom(atom).oatoms) == 2:
				return 28
			elif atom.bonds[0].xatom(atom).charge == 1:
				return 36
			#GENERAL H ON NITROGEN
			return 23
		if atom.bonds[0].xatom(atom).symbol in ['P','S']:
			return 71
		
	if atom.symbol == 'N':
		if (atom.rings) and (len(atom.rings) == 1) and (atom.aromatic) and (len(atom.rings[0].atoms) == 6) \
		   and ((len(atom.oatoms) == 3) or (atom.charge == 1)):
			#PYRIDINIUM-TYPE NITROGEN - FORMAL CHARGE=1
			return 58
		if len(atom.oatoms) == 1:
			if atom.bonds[0].bondtype == 3:
				#NITROGEN, TRIPLE BONDED
				return 42
			if (atom.bonds[0].bondtype == 2) and (atom.oatoms.symbol in ['N','C']):
				if (atom.charge == 0.5) and (atom.oatoms.symbol == 'C'):
					if (atom.bonds[0].bondtype == 1) and (atom.oatoms.symbol == 'N'):
							#N IN +N=C-N RESONANCE STRUCTURES - FORMAL CHARGE=1/2
							return 55
				elif (atom.charge < 0) and (atom.oatoms.symbol == 'N') and (atom.oatoms.charge > 0):
					for bond in atom.oatoms.bonds:
						if (bond.bondtype == 2) and (bond.xatom(atom.oatoms) != atom) and (bond.xatom(atom.oatoms).symbol == 'N'):
							#TERMINAL NITROGEN IN AZIDO OR DIAZO GROUP
							return 47
				
				elif atom.charge == 1:
					#POSITIVELY CHARGED NITROGEN DOUBLE-BONDED TO N, or IMINIUM NITROGEN
					return 54
		if len(atom.oatoms) == 2:
			for bond in atom.bonds:
				if (bond.bondtype == 3) and (bond.xatom(atom).symbol in ['C','N']):
					#ISONITRILE NITROGEN [FC = 0] OR DIAZO NITROGEN [FC = 1]
					return 61
			if atom.charge == 0:
				nOAtoms = cOAtoms = 0
				for bond in atom.bonds:
					if (bond.bondtype in [2,4]) and (bond.xatom(atom).symbol == 'N'):
						nOAtoms += 1
					if (bond.bondtype in [2,4]) and (bond.xatom(atom).symbol == 'C'):
						cOAtoms += 1
					if (bond.bondtype in [2,4]) and (bond.xatom(atom).symbol == 'O'):
						#NITROSO NITROGEN
						return 46
					if (bond.bondtype in [2,4]) and (bond.xatom(atom).symbol == 'S'):
						for bbond in bond.xatom(atom).bonds:
							if (bbond.bondtype in [2,4]) and (bbond.xatom(bond.xatom(atom)).symbol == 'O'):
								#DIVALENT NITROGEN REPLACING MONOVALENT O IN SO2 GROUP
								return 48
				if (nOAtoms == 1) and (cOAtoms == 1):
					#NITROGEN IN C=N=N OR -N=N=N
					return 53
				if (nOAtoms == 1) or (cOAtoms == 1):
					#NITROGEN IN IMINES or AZO COMPOUNDS
					return 9
			elif atom.charge == -1:
				for oatom in atom.oatoms:
					if oatom.symbol == 'S':
						if getMMFFAtomType(oatom) == 18:
							#DEPROTONATED SULFONAMIDE N-; FORMAL CHARGE=-1
							return 62
		if len(atom.oatoms) == 3:
			nbN = simpleO = ddoubleO = 0
			for bond in atom.bonds:
				if (bond.bondtype in [1,4]) and (bond.xatom(atom).symbol == 'N'):
					for bbond in bond.xatom(atom).bonds:
						if (bbond.bondtype in [2,4]) and (bbond.xatom(bond.xatom(atom)).symbol in ['C','N']):
							#NITROGEN IN N-N=C or N-N=N
							atomTypes.append(10)
				if (bond.bondtype in [1,4]) and (bond.xatom(atom).symbol == 'C'):
					for bbond in bond.xatom(atom).bonds:
						if bbond == bond:
							continue
						if (bbond.bondtype in [2,4]) and (bbond.xatom(bond.xatom(atom)).symbol in ['O','S']):
							#NITROGEN IN AMIDES, N-C=S or  THIOAMIDE
							return 10
							if bbond.xatom(bond.xatom(atom)).symbol == 'O':
								return 10
				if (bond.bondtype in [1,4]) and (bond.xatom(atom).symbol == 'O'):
					simpleO += 1
				if (bond.bondtype in [2,4]) and (bond.xatom(atom).symbol == 'O'):
					ddoubleO += 1
			if ((simpleO == 0) and (ddoubleO == 1)) or (simpleO + ddoubleO >= 2):
				#NITRO or NITRATE GROUP NITROGEN
				atomTypes.append(45)
			for bond in atom.bonds:
				if (bond.bondtype in [1,4]) and (bond.xatom(atom).symbol in ['S','P']):
					if len(bond.xatom(atom).oatoms) == 4:
						doubleO = 0
						for bbond in bond.xatom(atom).bonds:
							if (bond.bondtype in [2,4]) and (bond.xatom(atom).symbol == 'O'):
								doubleO += 1
						if doubleO == 2:
							#NITROGEN IN SULFONAMIDES or PHOSPHONAMIDES
							atomTypes.append(43)
				if (bond.bondtype in [1,4]) and (bond.xatom(atom).symbol == 'C'):
					for bbond in bond.xatom(atom).bonds:
						if bbond == bond:
							continue
						if (bbond.bondtype in [2,4]) and (bbond.xatom(bond.xatom(atom)).symbol in ['C','N','P']):
							#NITROGEN ON N-C=C or N-C=N or N-C=P
							atomTypes.append(40)
						if bbond.bondtype == 3:
							if bbond.xatom(bond.xatom(atom)).symbol == 'N':
								#NITROGEN ATTACHED TO CYANO GROUP
								atomTypes.append(43)
							if bbond.xatom(bond.xatom(atom)).symbol == 'C':
								#NITROGEN ATTACHED TO C-C TRIPLE BOND
								atomTypes.append(40)
						if (bbond.xatom(bond.xatom(atom)) != atom) and (bbond.bondtype == 1) \
						   and (bbond.xatom(bond.xatom(atom)).symbol == 'N'):
							if (len(bbond.xatom(bond.xatom(atom)).oatoms) == 3) \
							   and (getNbOfHInOatoms(bbond.xatom(bond.xatom(atom))) >=1 ):
								nbN += 1
			#	if (bond.bondtype == 1) and (bond.xatom(atom).symbol == 'O'):
			#		simpleO += 1
			#	if (bond.bondtype == 2) and (bond.xatom(atom).symbol == 'O'):
			#		ddoubleO += 1
			#if ((simpleO == 0) and (ddoubleO == 1)) or (simpleO + ddoubleO >= 2):
			#	#NITRO or NITRATE GROUP NITROGEN
			#	return 45
			if (nbN == 2) and (getNbOfHInOatoms(atom) >= 1):
				#GUANIDINIUM-TYPE NITROGEN
				atomTypes.append(56)
				return 56
			#NITROGEN IN ALIPHATIC AMINES
			return 8
		if (len(atom.oatoms) == 4) and (atom.charge == 1):
			for bond in atom.bonds:
				if bond.bondtype != 1:
					break
			else:
				#QUATERNARY NITROGEN, SP3, POSITIVELY CHARGED
				return 34
                if len(atom.oatoms) == 4:
                        return 8
		if (atom.rings) and (len(atom.rings) == 1):
			if atom.aromatic:
				if len(atom.rings[0].atoms) == 5:
					for oatom in atom.oatoms:
						if (oatom not in atom.rings[0].atoms) and (oatom.symbol == 'H'):
							#NITROGEN, AS IN PYRROLE
							atomTypes.append(39)
					#AROM HETEROCYCLIC 5-RING  NITROGEN
					atomTypes.append(65)
				if len(atom.rings[0].atoms) == 6:
					#if atom.charge == 1:
					if len(atom.oatoms) == 2:	
						for oatom in atom.oatoms:
							#if (oatom not in atom.rings[0].atoms) and (oatom.symbol == 'H'):
							#	#PYRIDINIUM-TYPE NITROGEN - FORMAL CHARGE=1
							#	return 58
							if (oatom not in atom.rings[0].atoms) and (oatom.symbol == 'O'):
								#PYRIDINE N-OXIDE NITROGEN
								atomTypes.append(69)
						#NITROGEN, AS IN PYRIDINE
						atomTypes.append(38)
					elif (len(atom.oatoms) == 3) or (atom.charge == 1):
						#PYRIDINIUM-TYPE NITROGEN - FORMAL CHARGE=1
						atomTypes.append(58)
			else:
				if len(atom.rings[0].atoms) == 5:
					if atom.charge <= 0:
						nbN = 0
						for atom in atom.rings[0].atoms:
							if atom.symbol == 'N':
								nbN += 1
						if nbN >= 3:
							#NEGATIVELY CHARGED N IN, E.G, TRI- OR TETRAZOLE ANION
							atomTypes.append(76)
					for oatom in atom.oatoms:
						if (oatom not in atom.rings[0].atoms) and (oatom.symbol == 'O'):
							#N-OXIDE NITROGEN IN GENERAL 5-RING POSITION
							atomTypes.append(82)
					#GENERAL NITROGEN IN 5-MEMBERED HETEROCYCLIC RING
					return 79
		#print "atomtypes", atomTypes
		if 45 in atomTypes:
			return 45
		if 40 in atomTypes:
			return 40
		if 9 in atomTypes:
			return 9
		if 10 in atomTypes:
			return 10
	if atomTypes: 
		if len(atomTypes) == 1:
			return atomTypes[0]
		return atomTypes
	return getSimplestMMFFAtomType(atom)
	raise ValueError, "MMFF Atom type not distinguished "+atom.symbol+" "+str(atom.index)


def getMMFFAtomType_ori(atom):
	"""get the type of an atom used in MMFF"""
	atomTypes = []
	if atom.symbol == 'C':
		if atom.rings:
			for ring in atom.rings:
				if len(ring.atoms) == 3:
					#CARBON IN A 3-MEMBERED RING
					if 22 not in atomTypes:
						atomTypes.append(22)
				if len(ring.atoms) == 4:
					for bond in atom.bonds:
						if (bond.xatom(atom) not in ring.atoms) and (bond.xatom(atom).symbol == 'C') and (bond.bondtype == 2):
							#OLEFINIC CARBON IN 4-MEMBERED RINGS
							if 30 not in atomTypes:
								atomTypes.append(30)
					#CARBON IN A 4-MEMBERED RING
					if 20 not in atomTypes:
						atomTypes.append(20)
				if (atom.aromatic) and (len(ring.atoms) == 6):
					#CARBON AS IN BENZENE, PYRROLE
					if 37 not in atomTypes:
						atomTypes.append(37)
				if (atom.aromatic) and (len(ring.atoms) == 5):
					#GENERAL CARBON IN 5-MEMBERED HETEROAROMATIC RING
					if 78 not in atomTypes:
						atomTypes.append(78)
		if (atom.charge == 0) and (len(atom.oatoms) == 2):
			#ACETYLENIC or ALLENIC CARBON
			atomTypes.append(4)
		if len(atom.oatoms) == 3:
			for bond in atom.bonds:
				if bond.bondtype in [2,4]:
					if bond.xatom(atom).symbol in ['S','P','N','O']:
						for bondd in atom.bonds:
							if (bondd.bondtype in [1,4]):
								if (bond.xatom(atom).charge == 1) and (bondd.xatom(atom).symbol == 'N'):
									#GUANIDINIUM CARBON or C IN +N=C-N RESONANCE STRUCTURES
									if 57 not in atomTypes:
										atomTypes.append(57)
								if (bond.xatom(atom).charge == -1) and (bondd.xatom(atom).symbol in ['O','S']):
									#CARBOXYLATE ANION CARBON or CARBON IN THIOCARBOXYLATE ANION
									if 41 not in atomTypes:
										atomTypes.append(41)
						#CARBONYL or CARBONYL-LIKE CARBON
						if 3 not in atomTypes:
							atomTypes.append(3)
						if atom.rings:
							for ring in atom.rings:
								if bond.xatom(atom) in ring.atoms:
									break
							else:
								#print "yoyoyo"
								return 3
					
			else:
				#VINYLIC or GENERIC SP2 CARBON
				atomTypes.append(2)
		if (len(atom.oatoms) == 1) and (atom.bonds[0].xatom(atom).symbol == 'N'):
			#ISONITRILE CARBON
			atomTypes.append(60)
		if len(atom.oatoms) == 4:
			nbN = 0
			for bond in atom.bonds:
				if bond.xatom(atom).symbol == 'N':
					nbN += 1
			if nbN == 2:
				#C IN N-C-N IN IMIDAZOLIUM ION
				atomTypes.append(80)
		if len(atom.oatoms) == 4:
			#ALKYL CARBON, SP3
			atomTypes.append(1)
		
		if len(atomTypes) == 1:
			return atomTypes[0]
		#if (37 in atomTypes) and (3 in atomTypes):
		#	return 3
		#print atomTypes
		if 37 in atomTypes:
			return 37
		if 3 in atomTypes:
			return 3
		if 30 in atomTypes:
			return 30
		if 20 in atomTypes:
			return 20
		if 22 in atomTypes:
			return 22
		if 2 in atomTypes:
			return 2
                if 80 in atomTypes:
                        return 80
	if atom.symbol == 'F':
		if atom.charge == -1:
			#FLUORIDE ANION
			return 89
		#FLUORINE
		return 11
	
	if atom.symbol == 'Cl':
		if atom.charge == -1:
			if len(atom.oatoms) == 3:
				for bond in atom.bonds:
					if not ((bond.xatom(atom).symbol == 'O') and (bond.bondtype == 2)):
						break
				else:
					#CHLORINE IN PERCHLORATE ANION, CLO4(-)
					return 77
			#CHLORIDE ANION
			return 90
		#CHLORINE
		return 12
			
	if atom.symbol == 'Br':
		if atom.charge == -1:
			#BROMIDE ANION
			return 91
		#BROMINE
		return 13
	
	if atom.symbol == 'I':
		#IODINE
		return 14
	
	if atom.symbol == 'Si':
		#SILICONE
		return 19
	
	if atom.symbol == 'Fe':
		if atom.charge == 2:
			#IRON +2 CATION
			return 87	
		if atom.charge == 3:
			#IRON +3 CATION
			return 88
		
	if (atom.symbol == 'Li') and (atom.charge == 1):
		#LITHIUM CATION
		return 92

	if (atom.symbol == 'Na') and (atom.charge == 1):
		#SODIUM CATION
		return 93
	
	if (atom.symbol == 'K') and (atom.charge == 1):
		#POTASSIUM CATION
		return 94
	
	if (atom.symbol == 'Zn') and (atom.charge == 2):
		#PDIPOSITIVE ZINC
		return 95
	
	if (atom.symbol == 'Ca') and (atom.charge == 2):
		#DIPOSITIVE CALCIUM
		return 96
	
	if atom.symbol == 'Cu':
		if atom.charge == 1:
			#MONOPOSITIVE COPPER
			return 97	
		if atom.charge == 2:
			#DIPOSITIVE COPPER
			return 98
	
	if (atom.symbol == 'Mg') and (atom.charge == 2):
		#DIPOSITIVE MAGNESIUM CATION
		return 99
	
	if atom.symbol == 'P':
		if len(atom.oatoms) == 4:
			#GENERAL TETRACOORDINATE PHOSPHORUS
			return 25
		if len(atom.oatoms) == 3:
			#TRICOORDINATE P, AS IN PHOSPHINES
			return 26
		if len(atom.oatoms) == 2:
			for bond in atom.bonds:
				if ((bond.xatom(atom).symbol == 'C') and (bond.bondtype == 2)):
					#PHOSPHOROUS DOUBLY BONDED TO CARBON
					return 75
				
	if atom.symbol == 'S':
		if (atom.rings) and (len(atom.rings) == 1) and (atom.aromatic) and (len(atom.rings[0].atoms) == 5):
			#SULFUR AS IN THIOPHENE
			return 44
		if len(atom.oatoms) == 1:
			if atom.charge == -1:
				#TERMINAL SULFUR - FORMAL CHARGE=-1
				return 72
			if (atom.bonds[0].bondtype == 2) and (atom.bonds[0].xatom(atom).symbol == 'C'):
				#TERMINAL SULFUR DOUBLY BONDED TO CARBON
				return 16	
		if len(atom.oatoms) == 2:
			#print "2"
			simple = simpleP = doubleC = doubleO = simpleOMoins = 0
			for bond in atom.bonds:
				if bond.bondtype == 1:
					if bond.xatom(atom).symbol == 'P':
						simpleP += 1
					elif (bond.xatom(atom).symbol == 'O') and (bond.xatom(atom).charge == -1):
						simpleOMoins += 1
					else:
						simple += 1
				elif bond.bondtype == 2:
					if bond.xatom(atom).symbol == 'O':
						doubleO += 1
					elif bond.xatom(atom).symbol == 'C':
						doubleC += 1
			if simple == 2:
				#SULFUR IN THIOETHERS AND MERCAPTANS
				return 15
			if (simple == 1) and (simpleP == 1):
				#TERMINAL SULFUR BONDED TO PHOSPHORUS
				return 72
			#if (simple == 1) and (doubleC == 1):
			#	#TERMINAL SULFUR DOUBLY BONDED TO CARBON
			#	return 16
			if (simpleOMoins == 1) and (doubleO == 1):
				#SULFUR IN NEGATIVELY CHARGED SULFINATE GROUP
				return 73
			if (doubleC == 1) and (doubleO == 1):
				#SULFINYL SULFUR, EG. IN C=S=O
				return 74
			if doubleC == 1:
				#TERMINAL SULFUR DOUBLY BONDED TO CARBON
				return 16
		if len(atom.oatoms) == 3:
			simple = simpleS = doubleC = doubleO = doubleN = simpleOMoins = 0
			for bond in atom.bonds:
				if bond.bondtype == 1:
					if bond.xatom(atom).symbol == 'S':
						simpleS += 1
					elif (bond.xatom(atom).symbol == 'O') and (bond.xatom(atom).charge == -1):
						simpleOMoins += 1
					else:
						simple += 1
				elif bond.bondtype == 2:
					if bond.xatom(atom).symbol == 'O':
						doubleO += 1	
					elif bond.xatom(atom).symbol == 'N':
						doubleN += 1
					elif bond.xatom(atom).symbol == 'C':
						doubleC += 1		
			if (doubleC == 1) and (doubleO + simpleOMoins == 2):
				#SULFONE SULPHER DOUBLY BONDED TO CARBON
				return 18
			if (simple == 2) and ((doubleN == 1) or (doubleO == 1)):
				#SULFUR, TRICOORD, DOUBLY BONDED TO N, or in SULFOXIDES
				return 17
			if (simpleS == 1) and (simpleOMoins == 1) and (doubleO == 1):
				#TRICOORD SULFUR IN THIOSULFINATE GROUP
				return 73
		if len(atom.oatoms) == 4:
			simple = simpleO = simpleN = doubleO = doubleN = 0
			for bond in atom.bonds:
				if bond.bondtype == 1:
					if bond.xatom(atom).symbol == 'O':
						simpleO += 1
					elif bond.xatom(atom).symbol == 'N':
						simpleN += 1
					else:
						simple += 1
				elif bond.bondtype == 2:
					if bond.xatom(atom).symbol == 'O':
						doubleO += 1	
					elif bond.xatom(atom).symbol == 'N':
						doubleN += 1
			if (simple == 2) and (doubleO + doubleN == 2):
				return 18
			if (simpleO + doubleO >= 2) and (simpleN >= 1):
				return 18
			if simpleO + doubleO >= 2:
				return 18
	if atom.symbol == 'H':
		if (atom.bonds[0].xatom(atom).symbol == 'C') or (atom.bonds[0].xatom(atom).symbol == 'Si'):
			#H  ATTACHED TO C or SI
			return 5
 		if atom.bonds[0].xatom(atom).symbol == 'O':
			if atom.bonds[0].xatom(atom).charge == 1:
				if len(atom.bonds[0].xatom(atom).oatoms) == 2:
					for bond in atom.bonds[0].xatom(atom).bonds:
						if bond.bondtype == 2:
							#HYDROGEN ON OXENIUM OXYGEN
							return 52
				#HYDROGEN ON O+ OXYGEN
				return 50
			if atom.bonds[0].xatom(atom).charge == -1:
				#HYDROGEN IN HYDROXIDE ANION
				return 21
			for oatom in atom.oatoms[0].oatoms:
				if oatom != atom:
					break
			if oatom.symbol == 'H':
				#HYDROGEN IN H2O
				return 31
			if oatom.symbol == 'S':
				#H ON OXYGEN ATTACHED TO SULFUR
				return 33
			if oatom.symbol == 'P':
				#HYDROGEN ON OXYGEN ATTACHED TO PHOSPHOROUS
				return 24
			if oatom.symbol == 'C':
				#if oatom.imp_hcount + oatom.explicit_hcount <= 2: #atom.bonds[0].xatom(atom) >= 2:
				#	#HYDROGEN IN ALCOHOLS
				#	return 21
				if oatom.rings:
					#H-O IN ENOLS AND PHENOLS
					return 29
				for bond in oatom.bonds:
					if bond.bondtype == 2:
						break
				if bond.xatom(oatom).symbol == 'O':
					#H-O IN CARBOXYLIC ACIDS
					return 24
				if bond.xatom(oatom).symbol == 'N':
					#H-O IN HO-C=N
					return 29
			#GENERAL H ON OXYGEN
			return 21
		if atom.bonds[0].xatom(atom).symbol == 'N':
			nType = getMMFFAtomType(atom.bonds[0].xatom(atom))
			if nType in [10, 40, 43]:
				return 28
			elif nType == 9:
				return 27
			elif len(atom.bonds[0].xatom(atom).oatoms) == 2:
				return 28
			elif atom.bonds[0].xatom(atom).charge == 1:
				return 36
			#GENERAL H ON NITROGEN
			return 23
		if atom.bonds[0].xatom(atom).symbol in ['P','S']:
			return 71
		
	if atom.symbol == 'O':
		if len(atom.oatoms) == 2:
			if getNbOfHInOatoms(atom) == 2:
				#OXYGEN ON WATER
				return 70
			if atom.aromatic:
				#AROMATIC OXYGEN AS IN FURAN
				return 59
			if atom.charge == 1:
				for bond in atom.bonds:
					if bond.bondtype == 2:
						#POSITIVELY CHARGED OXENIUM (DICOORDINATE) OXYGEN
						return 51
			#DIVALENT OXYGEN 
			return 6
		if (len(atom.oatoms) == 3) and (atom.charge == 1):
			#POSITIVELY CHARGED OXONIUM (TRICOORDINATE) OXYGEN
			return 49
		if (len(atom.oatoms) == 1) and (atom.oatoms[0].symbol in ['C','N']):
			if (atom.charge == -1) and (len(atom.oatoms[0].oatoms) == 3):
				if atom.oatoms[0].symbol == 'C':
					for bond in atom.oatoms[0].bonds:
						if (bond.bondtype == 2) and (bond.xatom(atom.oatoms[0]).symbol == 'O'):
							#OXYGEN IN CARBOXYLATE ANION
							return 32
				if (atom.oatoms[0].symbol == 'N') and (atom.oatoms[0].charge != 0):
					return 32
				#OXIDE OXYGEN ON SP2 CARBON, NEGATIVELY CHARGED
				return 35
			if (atom.charge == -1):
				#ALKOXIDE OXYGEN, NEGATIVELY CHARGED
				return 35
			if (atom.bonds[0].bondtype == 2) and (atom.oatoms[0].symbol == 'C'):
				#GENERAL C=O, OR CARBONYL-LIKE OXYGEN
				return 7
		if (len(atom.oatoms) == 1) and (atom.oatoms[0].symbol == 'S'):
			if atom.bonds[0].bondtype == 2:
				if len(atom.oatoms[0].oatoms) == 3:
					nbSimpleBonds = 0
					for bond in atom.oatoms[0].bonds:
						if (bond.bondtype == 1):
							nbSimpleBonds += 1
					if nbSimpleBonds == 2:
						#O=S IN SULFOXIDES
						return 7
				if len(atom.oatoms[0].oatoms) == 2:
					for bond in atom.oatoms[0].bonds:
						if (bond.bondtype == 2) and (bond.xatom(atom.oatoms[0]).symbol == 'C'):
							#O=S ON SULFUR DOUBLY BONDED TO, E.G., CARBON
							return 7
		return 32
		
	if atom.symbol == 'N':
		if (atom.rings) and (len(atom.rings) == 1) and (atom.aromatic) and (len(atom.rings[0].atoms) == 6) \
		   and ((len(atom.oatoms) == 3) or (atom.charge == 1)):
			#PYRIDINIUM-TYPE NITROGEN - FORMAL CHARGE=1
			return 58
		if len(atom.oatoms) == 1:
			if atom.bonds[0].bondtype == 3:
				#NITROGEN, TRIPLE BONDED
				return 42
			if (atom.bonds[0].bondtype == 2) and (atom.oatoms.symbol in ['N','C']):
				if (atom.charge == 0.5) and (atom.oatoms.symbol == 'C'):
					if (atom.bonds[0].bondtype == 1) and (atom.oatoms.symbol == 'N'):
							#N IN +N=C-N RESONANCE STRUCTURES - FORMAL CHARGE=1/2
							return 55
				elif (atom.charge < 0) and (atom.oatoms.symbol == 'N') and (atom.oatoms.charge > 0):
					for bond in atom.oatoms.bonds:
						if (bond.bondtype == 2) and (bond.xatom(atom.oatoms) != atom) and (bond.xatom(atom.oatoms).symbol == 'N'):
							#TERMINAL NITROGEN IN AZIDO OR DIAZO GROUP
							return 47
				
				elif atom.charge == 1:
					#POSITIVELY CHARGED NITROGEN DOUBLE-BONDED TO N, or IMINIUM NITROGEN
					return 54
		if len(atom.oatoms) == 2:
			for bond in atom.bonds:
				if (bond.bondtype == 3) and (bond.xatom(atom).symbol in ['C','N']):
					#ISONITRILE NITROGEN [FC = 0] OR DIAZO NITROGEN [FC = 1]
					return 61
			if atom.charge == 0:
				nOAtoms = cOAtoms = 0
				for bond in atom.bonds:
					if (bond.bondtype in [2,4]) and (bond.xatom(atom).symbol == 'N'):
						nOAtoms += 1
					if (bond.bondtype in [2,4]) and (bond.xatom(atom).symbol == 'C'):
						cOAtoms += 1
					if (bond.bondtype in [2,4]) and (bond.xatom(atom).symbol == 'O'):
						#NITROSO NITROGEN
						return 46
					if (bond.bondtype in [2,4]) and (bond.xatom(atom).symbol == 'S'):
						for bbond in bond.xatom(atom).bonds:
							if (bbond.bondtype in [2,4]) and (bbond.xatom(bond.xatom(atom)).symbol == 'O'):
								#DIVALENT NITROGEN REPLACING MONOVALENT O IN SO2 GROUP
								return 48
				if (nOAtoms == 1) and (cOAtoms == 1):
					#NITROGEN IN C=N=N OR -N=N=N
					return 53
				if (nOAtoms == 1) or (cOAtoms == 1):
					#NITROGEN IN IMINES or AZO COMPOUNDS
					return 9
			elif atom.charge == -1:
				for oatom in atom.oatoms:
					if oatom.symbol == 'S':
						if getMMFFAtomType(oatom) == 18:
							#DEPROTONATED SULFONAMIDE N-; FORMAL CHARGE=-1
							return 62
		if len(atom.oatoms) == 3:
			nbN = simpleO = ddoubleO = 0
			for bond in atom.bonds:
				if (bond.bondtype in [1,4]) and (bond.xatom(atom).symbol == 'N'):
					for bbond in bond.xatom(atom).bonds:
						if (bbond.bondtype in [2,4]) and (bbond.xatom(bond.xatom(atom)).symbol in ['C','N']):
							#NITROGEN IN N-N=C or N-N=N
							atomTypes.append(10)
				if (bond.bondtype in [1,4]) and (bond.xatom(atom).symbol == 'C'):
					for bbond in bond.xatom(atom).bonds:
						if bbond == bond:
							continue
						if (bbond.bondtype in [2,4]) and (bbond.xatom(bond.xatom(atom)).symbol in ['O','S']):
							#NITROGEN IN AMIDES, N-C=S or  THIOAMIDE
							return 10
							if bbond.xatom(bond.xatom(atom)).symbol == 'O':
								return 10
				if (bond.bondtype in [1,4]) and (bond.xatom(atom).symbol == 'O'):
					simpleO += 1
				if (bond.bondtype in [2,4]) and (bond.xatom(atom).symbol == 'O'):
					ddoubleO += 1
			if ((simpleO == 0) and (ddoubleO == 1)) or (simpleO + ddoubleO >= 2):
				#NITRO or NITRATE GROUP NITROGEN
				atomTypes.append(45)
			for bond in atom.bonds:
				if (bond.bondtype in [1,4]) and (bond.xatom(atom).symbol in ['S','P']):
					if len(bond.xatom(atom).oatoms) == 4:
						doubleO = 0
						for bbond in bond.xatom(atom).bonds:
							if (bond.bondtype in [2,4]) and (bond.xatom(atom).symbol == 'O'):
								doubleO += 1
						if doubleO == 2:
							#NITROGEN IN SULFONAMIDES or PHOSPHONAMIDES
							atomTypes.append(43)
				if (bond.bondtype in [1,4]) and (bond.xatom(atom).symbol == 'C'):
					for bbond in bond.xatom(atom).bonds:
						if bbond == bond:
							continue
						if (bbond.bondtype in [2,4]) and (bbond.xatom(bond.xatom(atom)).symbol in ['C','N','P']):
							#NITROGEN ON N-C=C or N-C=N or N-C=P
							atomTypes.append(40)
						if bbond.bondtype == 3:
							if bbond.xatom(bond.xatom(atom)).symbol == 'N':
								#NITROGEN ATTACHED TO CYANO GROUP
								atomTypes.append(43)
							if bbond.xatom(bond.xatom(atom)).symbol == 'C':
								#NITROGEN ATTACHED TO C-C TRIPLE BOND
								atomTypes.append(40)
						if (bbond.xatom(bond.xatom(atom)) != atom) and (bbond.bondtype == 1) \
						   and (bbond.xatom(bond.xatom(atom)).symbol == 'N'):
							if (len(bbond.xatom(bond.xatom(atom)).oatoms) == 3) \
							   and (getNbOfHInOatoms(bbond.xatom(bond.xatom(atom))) >=1 ):
								nbN += 1
			#	if (bond.bondtype == 1) and (bond.xatom(atom).symbol == 'O'):
			#		simpleO += 1
			#	if (bond.bondtype == 2) and (bond.xatom(atom).symbol == 'O'):
			#		ddoubleO += 1
			#if ((simpleO == 0) and (ddoubleO == 1)) or (simpleO + ddoubleO >= 2):
			#	#NITRO or NITRATE GROUP NITROGEN
			#	return 45
			if (nbN == 2) and (getNbOfHInOatoms(atom) >= 1):
				#GUANIDINIUM-TYPE NITROGEN
				atomTypes.append(56)
				return 56
			#NITROGEN IN ALIPHATIC AMINES
			return 8
		if (len(atom.oatoms) == 4) and (atom.charge == 1):
			for bond in atom.bonds:
				if bond.bondtype != 1:
					break
			else:
				#QUATERNARY NITROGEN, SP3, POSITIVELY CHARGED
				return 34
                if len(atom.oatoms) == 4:
                        return 8
		if (atom.rings) and (len(atom.rings) == 1):
			if atom.aromatic:
				if len(atom.rings[0].atoms) == 5:
					for oatom in atom.oatoms:
						if (oatom not in atom.rings[0].atoms) and (oatom.symbol == 'H'):
							#NITROGEN, AS IN PYRROLE
							atomTypes.append(39)
					#AROM HETEROCYCLIC 5-RING  NITROGEN
					atomTypes.append(65)
				if len(atom.rings[0].atoms) == 6:
					#if atom.charge == 1:
					if len(atom.oatoms) == 2:	
						for oatom in atom.oatoms:
							#if (oatom not in atom.rings[0].atoms) and (oatom.symbol == 'H'):
							#	#PYRIDINIUM-TYPE NITROGEN - FORMAL CHARGE=1
							#	return 58
							if (oatom not in atom.rings[0].atoms) and (oatom.symbol == 'O'):
								#PYRIDINE N-OXIDE NITROGEN
								atomTypes.append(69)
						#NITROGEN, AS IN PYRIDINE
						atomTypes.append(38)
					elif (len(atom.oatoms) == 3) or (atom.charge == 1):
						#PYRIDINIUM-TYPE NITROGEN - FORMAL CHARGE=1
						atomTypes.append(58)
			else:
				if len(atom.rings[0].atoms) == 5:
					if atom.charge <= 0:
						nbN = 0
						for atom in atom.rings[0].atoms:
							if atom.symbol == 'N':
								nbN += 1
						if nbN >= 3:
							#NEGATIVELY CHARGED N IN, E.G, TRI- OR TETRAZOLE ANION
							atomTypes.append(76)
					for oatom in atom.oatoms:
						if (oatom not in atom.rings[0].atoms) and (oatom.symbol == 'O'):
							#N-OXIDE NITROGEN IN GENERAL 5-RING POSITION
							atomTypes.append(82)
					#GENERAL NITROGEN IN 5-MEMBERED HETEROCYCLIC RING
					return 79
		#print "atomtypes", atomTypes
		if 45 in atomTypes:
			return 45
		if 40 in atomTypes:
			return 40
		if 9 in atomTypes:
			return 9
		if 10 in atomTypes:
			return 10
	if atomTypes: 
		if len(atomTypes) == 1:
			return atomTypes[0]
		return atomTypes
	return getSimplestMMFFAtomType(atom)
	raise ValueError, "MMFF Atom type not distinguished "+atom.symbol+" "+str(atom.index)



def getSimplestMMFFAtomType(atom):
	"""get the simplest type of an atom used in MMFF"""
	if atom.symbol == 'C':
		if len(atom.oatoms) == 2:
			return 4
		if len(atom.oatoms) == 3:
			return 2
		if len(atom.oatoms) == 4:
			return 1
		return 1
		
	if atom.symbol == 'F':
		if atom.charge == -1:
			return 89
		return 11
	
	if atom.symbol == 'Cl':
		return 12
	
	if atom.symbol == 'Br':
		if atom.charge == -1:
			return 91
		return 13
	
	if atom.symbol == 'I':
		return 14
	
	if atom.symbol == 'Si':
		return 19
	
	if atom.symbol == 'Fe':
		if atom.charge < 3:
			return 87	
		return 88
		
	if atom.symbol == 'Li':
		return 92

	if atom.symbol == 'Na':
      		return 93
	
	if atom.symbol == 'K':
		return 94
	
	if atom.symbol == 'Zn':
		return 95
	
	if atom.symbol == 'Ca':
		return 96
	
	if atom.symbol == 'Cu':
		if atom.charge <= 1:
			return 97	
		if atom.charge >= 2:
			return 98
	
	if atom.symbol == 'Mg':
		return 99
	
	if atom.symbol == 'P':
		if len(atom.oatoms) == 4:
			return 25
		if len(atom.oatoms) == 3:
			return 26
		if len(atom.oatoms) == 2:
			return 75
				
	if atom.symbol == 'S':
		if len(atom.oatoms) == 1:
			return 16
		if len(atom.oatoms) == 2:
			return 16
		if len(atom.oatoms) == 3:
			return 18
		return 18
			
	if atom.symbol == 'H':
		if atom.bonds[0].xatom(atom).symbol in ['C','Si']:
			return 5
 		if atom.bonds[0].xatom(atom).symbol == 'O':
			return 21
		if atom.bonds[0].xatom(atom).symbol == 'N':
			return 23
		if atom.bonds[0].xatom(atom).symbol in ['P','S']:
			return 71
		return 5
		
	if atom.symbol == 'O':
		if len(atom.oatoms) == 2:
			return 6
		if len(atom.oatoms) == 1:
			return 7
		return 6
		
	if atom.symbol == 'N':
		if len(atom.oatoms) == 1:
			return 42
		if len(atom.oatoms) == 2:
			return 9
		if len(atom.oatoms) == 3:
			return 8
		if len(atom.oatoms) == 4:
			return 34
		return 34

	raise ValueError, "MMFF Atom type not distinguished"


def getAllMMFFAtomsTypesForSymbol(atom):
	if atom.symbol == 'C':
                return [1,2,3,4,20,22,30,37,41,57,60,63,64,78,80]
		
	if atom.symbol == 'F':
                return [11,89]
	
	if atom.symbol == 'Cl':
		return [12]
	
	if atom.symbol == 'Br':
                return [13,91]
	
	if atom.symbol == 'I':
		return [14]
	
	if atom.symbol == 'Si':
		return [19]
	
	if atom.symbol == 'Fe':
                return [87,88]
		
	if atom.symbol == 'Li':
		return [92]

	if atom.symbol == 'Na':
      		return [93]
	
	if atom.symbol == 'K':
		return [94]
	
	if atom.symbol == 'Zn':
		return [95]
	
	if atom.symbol == 'Ca':
		return [96]
	
	if atom.symbol == 'Cu':
                return [97,98]
	
	if atom.symbol == 'Mg':
		return [99]
	
	if atom.symbol == 'P':
                return [25,26,75]
				
	if atom.symbol == 'S':
                return [15,16,17,18,44,72,73,74]
			
	if atom.symbol == 'H':
		return [5,21,23,24,27,28,29,31,33,36,50,52,71]
		
	if atom.symbol == 'O':
               return [6,7,32,35,49,51,59,70]
		
	if atom.symbol == 'N':
                return [8,9,10,34,38,39,40,42,43,45,46]

	raise ValueError, "MMFF Atom type not distinguished"



def getEnghHuberAtomType(atom):

	if atom.symbol == 'O':
		if len(atom.oatoms) == 2:
			if getNbOfHInOatoms(atom) == 1 and (atom.oatoms[0].symbol == 'C' or atom.oatoms[1].symbol == 'C'):
				return "OH1"
		elif len(atom.oatoms) == 1 and atom.oatoms[0].symbol == 'C':
			for oatom in atom.oatoms[0].oatoms:
				if oatom != atom and oatom.symbol == 'O' and oatom.findbond(atom.oatoms[0]).bondtype == 1:
					return "OC"
			if len(atom.oatoms[0].oatoms) == 3:
				return "O"

	
	if atom.symbol == 'S':
		if len(atom.oatoms[0].oatoms) == 2 and getNbOfHInOatoms(atom) == 1:
			return "SH1E"
		if len(atom.oatoms) > 1 and len(atom.oatoms[0].oatoms) == 2 and atom.oatoms[0].symbol == 'C' and atom.oatoms[1].symbol == 'C':
			return "SM"
		return "S"

	if atom.symbol == 'N':
		if getNbOfHInOatoms(atom) == 3:
			return "NH3"
		if getNbOfHInOatoms(atom) == 2:
			return "NH2"
		if getNbOfHInOatoms(atom) == 1:
			return "NH1"
		return "N"

	if atom.symbol == 'C':
		if getNbOfHInOatoms(atom) == 1 and atom.aromatic:
			return "CR1E"
		if getNbOfHInOatoms(atom) == 3 and len(atom.oatoms) == 4:
			return "CH3E"
		if getNbOfHInOatoms(atom) == 2 and len(atom.oatoms) == 4:
			return "CH2E"
		if getNbOfHInOatoms(atom) == 1 and len(atom.oatoms) == 4:
			return "CH1E"
		if len(atom.oatoms) == 3:
			simpleO = doubleO = 0
			soc = doc = False
			for bond in atom.bonds:
				if bond.bondtype == 1 and bond.xatom(atom).symbol == 'O':
					simpleO += 1
					if bond.xatom(atom).charge == 0:
						soc = True
				if bond.bondtype == 2 and bond.xatom(atom).symbol == 'O':
					doubleO += 1
					if bond.xatom(atom).charge == 0:
						doc = True
			if simpleO == 1:
				if doubleO == 1 and soc and doc:
					return "CN"
				return "C"
				
	return None

