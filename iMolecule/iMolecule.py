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


from frowns import Smiles
from frowns import Atom
from frowns import Bond
import copy
#import  tools
from Exceptions import *
from Brique import *
from frowns.Utils import MolecularWeight
from frowns.Utils import Components
from tools import *
from BestFit import *
from Diedre import *
from AxeRot import *
from math import sqrt, fabs
import math
import os
import time, random
from subprocess import Popen, PIPE
from parametres import *

NB_MAX_CONFS_TO_EXPLORE = 1500

from Config import *
from AMMOS import *

def babel_convert(fname1, fmt1, fname2, fmt2, verbose = 0):
	# cmd = "%s -i%s  -omol2 %s/%s.mol2" % (BABEL_PATH, fmt1, fname1, fmt2, fname2)
	cmd = "%s -i%s %s -o%s %s 2> /dev/null" % (BABEL_PATH, fmt1, fname1, fmt2, fname2)
	if verbose:
		sys.stderr.write("%s\n" % cmd)
	status = os.system(cmd)

def babel_reprotonate_mol2_old(filePrefix,oriExt=".ori.mol2", endExt=".mol2", verbose = False):
	"""
	A wrapper to delete hydrogens then regenerate them
	Will work for openbabel < 2.2.0 (e.g. 2.0.2)
	PT 13/01/2011 based on suggestion by DOM 23/09/2010

	NOTE: NOT COMPATIBLE WITH SEVERAL CALLS OF THIS IN SAME DIR ... SHOULD USE TEMPFILES
	"""
	# Delete hydrogens
	cmd = "%s -d -imol2 %s -omol2 %s 2> /dev/null" % (BABEL_PATH,filePrefix+oriExt, filePrefix+".babelskip.mol2")
	if not verbose:
		cmd = "%s 2> /dev/null" % cmd
	os.system(cmd)
	# Build hydrogens
	cmd = "%s -p -imol2 %s -omol2 %s 2> /dev/null" % (BABEL_PATH,filePrefix+".babelskip.mol2", filePrefix+endExt)
	if not verbose:
		cmd = "%s 2> /dev/null" % cmd
	os.system(cmd)
	

def babel_reprotonate_mol2(filePrefix,oriExt=".ori.mol2", endExt=".mol2", verbose = False):
	"""
	A wrapper to delete hydrogens then regenerate them
	Will work for openbabel > 2.2.0
	PT 13/01/2011 absed on suggestion by DOM 23/09/2010

	NOTE: NOT COMPATIBLE WITH SEVERAL CALLS OF THIS IN SAME DIR ... SHOULD USE TEMPFILES
	"""
	cmd = "%s -imol2 %s -omol2 %s 2> /dev/null" % (BABEL_PATH,filePrefix+oriExt, filePrefix+".babelskip.mol2")
	if not verbose:
		cmd = "%s 2> /dev/null" % cmd
	os.system(cmd)
	cmd = "%s -h -imol2 %s -omol2 %s 2> /dev/null" % (BABEL_PATH,filePrefix+".babelskip.mol2", filePrefix+endExt)
	if not verbose:
		cmd = "%s 2> /dev/null" % cmd
	os.system(cmd)


class iMolecule :
	
	def __init__(self, smiles=None, sdf=None, sdfMolecule =None, id=None, axialEquatorial=None, createFromSdf=False, extractCoordsFromSdf=False, sampleStereo = False, verbose = False):
		"""
		Create an instance of iMolecule.
		@param smiles      : generate from smiles
		@param sdf         : generate from sdf
		@param sdfMolecule : generate from sdf frowns molecule instance.
		@param id          : compound identifier (if None, will be assessed from smiles or sdf)
		@param axialEquatorial: flag axialEquatorial information for further generation
		@param createFromSdf : convert smiles into sdf and back to smiles (solves some problematic cases)
		@param extractCoordsFromSdf: if True, we just cconsider sdf is 3D sdf and we do not revert to smiles

		This will generate a graph() of the compound.

		@return iMolecule instance.
		"""
		self.graphCreated = False
		self.diedres = []
		if smiles and not createFromSdf:
			self.smiles = smiles
			try:
				if verbose:
					sys.stderr.write("iMolecule: Attempting molecule instance from %s\n" % smiles)
				aMolecule = Smiles.smilin(smiles)
				self.molecule = aMolecule
				# self.molecule = Smiles.smilin(smiles) 
				if verbose:
					sys.stderr.write("iMolecule: molecule instance OK\n")
				# sys.exit(0)
			except Exception, err:
				sys.stderr.write("iMolecule: smilin failed\n")
				# raise BadSmile, "Some problem occurred with the smiles " + self.smiles + " (Might be not conform)"
				
			if verbose:
				sys.stderr.write("iMolecule: axialequatorial\n")
			self.axialEquatorial = axialEquatorial
		elif smiles:
			# Convert into sdf and then back to smiles
			sdfFile = open("tmpSmiForSdf.smi", "w")
			sdfFile.write(smiles)
			if id:
				sdfFile.write(" " + id + "\n")
			else:
				sdfFile.write("\n")
			sdfFile.close()
			os.system(BABEL_PATH + " -ismi tmpSmiForSdf.smi -osdf tmpSdfFromSmi.sdf 2> /dev/null")
			os.system("echo '$$$$' >> tmpSdfFromSmi.sdf")
			sdfFile = open("tmpSdfFromSmi.smi", "w")
			reader = MDL.sdin(sdfFile)
			self.smiles = reader.next()[0].cansmiles()
			self.molecule = Smiles.smilin(self.smiles)
			sdfFile.close()
			os.system("rm -f tmpSdfFromSmi.sdf")
		elif sdf:
			fileSdf = open(sdf)
			reader = MDL.sdin(fileSdf)
			try:
				if not extractCoordsFromSdf:
					self.smiles = reader.next()[0].cansmiles()
					self.molecule = Smiles.smilin(self.smiles)
				else:
					self.molecule = reader.next()[0]
					self.smiles = self.molecule.cansmiles()
			
			except:
				raise BadSdf, "some problem occurred while readin the sdf (check its format)"
			
			fileSdf.close()
			self.axialEquatorial = axialEquatorial		
		elif sdfMolecule:
			# We have a frowns Molecule instance: we propagate.
			try:
				self.molecule = sdfMolecule
				self.smiles = self.molecule.cansmiles()
				self.diedres = self.getDiedres()
				self.defaultDihedrals(setAngle = False) # Do not affect dihedrals ...
				# for at in self.molecule.atoms:
				# 	print at.x, at.y, at.z
				# sys.exit(0)
				self.graphCreated = True
			except:
				raise BadSdf, "some problem occurred while readin the sdf (check its format)"
		
		if id:
			self.id = id
		else:
                        self.id = self.molecule.cansmiles()
			
		# for at in self.molecule.atoms:
		# 	print at.x, at.y, at.z
		# sys.exit(0)

		if not sdfMolecule:
			self.numberOfNonHAtoms = len(self.molecule.atoms)
			self.numberOfBondsWithoutH = len(self.molecule.bonds)

			#for at in self.molecule.atoms:
			#	print at, at.symbol, at.hcount, at.oatoms
			for i in range(self.numberOfNonHAtoms):
				hcount = self.molecule.atoms[i].hcount
				#print hcount
				at = self.molecule.atoms[i]
				if at.symbol == "C" and len(at.oatoms) < 4:
					oats = 0
					for bond in at.bonds:
						oats += bond.bondtype
					if at.charge <= -1:
						hcount = 4 - oats - at.charge
					else:
						hcount = 4 - oats + at.charge


				#print self.molecule.atoms[i], self.molecule.atoms[i].hcount, self.molecule.atoms[i].oatoms
				j = 0
				while j < hcount:
					h = Atom.Atom()
					h.set_symbol("H")
					self.molecule.add_atom(h)
					newBond = Bond.Bond()
					self.molecule.add_bond(newBond,self.molecule.atoms[i],h)
					j += 1
				if (hcount == 1) and (self.molecule.atoms[i].chiral_class != 0):
					hPos = 1
					if i == 0:
						hPos = 0
					self.molecule.atoms[i].oatoms.insert(hPos,self.molecule.atoms[i].oatoms.pop())
					self.molecule.atoms[i].bonds.insert(hPos,self.molecule.atoms[i].bonds.pop())
				self.molecule.atoms[i].hcount = 0

			for i in range(self.numberOfBondsWithoutH):
				nBondsZE = 0
				if self.molecule.bonds[i].bondtype == 2 or self.molecule.bonds[i].bondtype == 4 :
					bonds = []
					for atom in self.molecule.bonds[i].atoms:
						for bond in atom.bonds:
							if bond.symbol in ['\\','/']:
								bonds.append((bond, bond.symbol))
								#bondss.append((bond, bond.symbol))
					if len(bonds) == 2:
						stereo = 2 #trans
						if bonds[0][1] == bonds[1][1]:
							stereo = 1 #cis
						self.molecule.bonds[i].dbo = []
						self.molecule.bonds[i].dbo.append((bonds[0][0], bonds[1][0], stereo))
						nBondsZE += 1
		"""nSlash = 0
		for i in range(len(self.smiles)):
			if self.smiles[i] in ['\\','/']:
				nSlash += 1
		if nSlash and not nBondsZE:
			self.smiles = self.smiles.replace("/","")
			self.smiles = self.smiles.replace("\\","")
		"""	

		# Sample stereo sites (so not to get one side substituents)
		if sampleStereo:
			toggle = 0
			asa = self.ambiguousStereoAtoms(verbose = verbose)
			if verbose:
				sys.stderr.write("sampleStero for :")
			for a in asa:
				# a.chiral_class = random.randint(1,2)
				a.chiral_class = toggle + 1
				toggle = (toggle + 1) % 2
				if verbose:
					sys.stderr.write("%s %d -> %d" % (a.symbol, a.index, a.chiral_class))
			if verbose:
				sys.stderr.write("\n")
		for at in self.molecule.atoms:
			at.mmffAtomType = getMMFFAtomType(at)
		self.brickgraph = self.graph()
		self.toxic = None
		if verbose:
			sys.stderr.write( "iMolecule __INIT__ DONE\n")
			for atom in self.molecule.atoms:
				sys.stderr.write("%s %s: chiral %s symorder %s\n" % (str(atom.symbol), str(atom.index), str(atom.chiral_class), str(atom.symorder)))
			for bond in self.molecule.bonds:
				if bond.dbo:
					sys.stderr.write("Bond %s-%s-%s %s-%s-%s %s : symbol %s dbo %s\n" % (str(bond.atoms[0].symbol), str(bond.atoms[0].index), str(bond.atoms[0].handle), str(bond.atoms[1].symbol), str(bond.atoms[1].index), str(bond.atoms[1].handle), str(bond.bondtype), str(bond.symbol), str(bond.dbo)) )
				else:
					sys.stderr.write("Bond %s-%s-%s %s-%s-%s %s : symbol %s\n" % (str(bond.atoms[0].symbol), str(bond.atoms[0].index), str(bond.atoms[0].handle), str(bond.atoms[1].symbol), str(bond.atoms[1].index), str(bond.atoms[1].handle), str(bond.bondtype), str(bond.symbol)) )
			sys.stderr.write("Smiles: %s\n" % self.smiles)
			# sys.exit(0)

		# Report about ring size
		rings = self.detectFusedRings()
		largestRingSize = 0
		naromatic = 0
		for ring in rings:
			if len(ring) > largestRingSize:
				largestRingSize = len(ring)
				naromatic = 0
				for atm in ring:
					if atm.aromatic:
						naromatic += 1

		# Save correspondance between index and handle
		self.saveindex = {}
		for atom in self.molecule.atoms:
			self.saveindex[atom.handle] = atom.index

		# if largestRingSize > 18:
		self.bigRingWarning = False
		if (largestRingSize > MAXWARNRINGSIZE) and (naromatic < MAXWARNRINGAROMATIC):
			self.bigRingWarning = True
		if verbose:
			sys.stdout.write("%s Largest Ring is %d aromatic %d\n" % (self.id, largestRingSize, naromatic))
		# sys.exit(0)

	def getSaveIndex(handle):
		try:
			return self.saveindex[handle]
		except:
			return None
				
	def molecular_weight(self) :
		"""
		calcul poids moleculaire
		"""
		return MolecularWeight.computeMW(self.molecule)[1]

	"""
	def donnors(self) :

	def acceptor(self) :

	def logP(self) :
	"""
	
	def tpsa(self) :
		"""
		calculation of polar surface area based on fragment contributions (TPSA)
		"""
		mol = self.molecule
		atoms = mol.atoms
		rings = mol.cycles

		psa = 0.
		for atom in atoms :
			an = atom.number
			if not an == 8 and not an == 7 :
				continue
			nh = atom.hcount
			q = atom.charge
			isAromatic = atom.aromatic
			# checking whether this atom is in a 3-membered ring
			isIn3ring = 0
			for ring in rings :
				if len(ring) > 3 :
					continue
				ratoms = ring.atoms
				for ratom in ratoms :
					if atom == ratom :
						isIn3ring = 1
						break
				if isIn3ring :
					break
			nv = 0 # number of neighbours
			# finding the number of -,=,#,: bonds originating from the atom
			nsingle = ndouble = ntriple = naromatic = 0
			bonds = atom.bonds
			for bond in bonds :
				nv += 1
				bt = bond.bondtype
				if bt == 1 :
					nsingle += 1
				elif bt == 2 :
					ndouble += 1
				elif bt == 3 :
					ntriple += 1
				else :
					naromatic += 1
			# finding the psa conrtibution for this fragment (atom with hydrogens)
			p = -1
			if an == 7 :
				if nv == 1 :
					if nh == 0 and ntriple == 1 and q == 0 :
						p = 23.79 # N#
					elif nh == 1 and ndouble == 1 and q == 0 :
						p = 23.85 # [NH]=
					elif nh == 2 and nsingle == 1 and q == 0 :
						p = 26.02 # [NH2]-
					elif nh == 2 and ndouble == 1 and q == 1 :
						p = 25.59 # [NH2+]=
					elif nh == 3 and nsingle == 1 and q == 1 :
						p = 27.64 # [NH3+]-
				elif nv == 2 :
					if nh == 0 and nsingle == 1 and ndouble == 1 and q == 0 :
						p = 12.36 # =N-
					elif nh == 0 and ntriple == 1 and ndouble == 1 and q ==0 :
						p = 13.6 # =N#
					elif nh == 1 and nsingle == 2 and q == 0 and not isIn3ring :
						p = 12.03 # -[NH]-
					elif nh == 1 and nsingle == 2 and q == 0 and isIn3ring :
						p = 21.94 # -[NH]-r3
					elif nh == 0 and ntriple == 1 and nsingle == 1 and q == 1 :
						p = 4.36 # -[N+]#
					elif nh == 1 and ndouble == 1 and nsingle == 1 and q == 1 :
						p = 13.97 # -[NH+]=
					elif nh == 2 and nsingle == 2 and q == 1 :
						p = 16.61 # -[NH2+]-
					elif nh == 0 and naromatic == 2 and q == 0 :
						p = 12.89 # :[n]:
					elif nh == 1 and naromatic == 2 and q == 0 :
						p = 15.79 # :[nH]:
					elif nh == 1 and naromatic == 2 and q == 1 :
						p = 14.14 # :[nH+]:
				elif nv == 3 :
					if nh == 0 and nsingle == 3 and q == 0 and not isIn3ring :
						p = 3.24 # -N(-)-
					elif nh == 0 and nsingle == 3 and q == 0 and inIn3ring :
						p = 3.01 # -N(-)-r3
					elif nh == 0 and nsingle == 1 and ndouble == 2 and q == 0 :
						p = 11.68 # -N(=)=
					elif nh == 0 and nsingle == 2 and ndouble == 1 and q == 1 :
						p = 3.01 # =[N+](-)-
					elif nh == 1 and nsingle == 3 and q == 1 :
						p = 4.44 # -[NH+](-)-
					elif nh == 0 and naromatic == 3 and q == 0 :
						p = 4.41 # :[n](:):
					elif nh == 0 and nsingle == 1 and naromatic == 2 and q == 0 :
						p = 4.93 # -:[n](:):
					elif nh == 0 and ndouble == 1 and naromatic == 2 and q == 0 :
						p = 8.39 # =:[n](:):
					elif nh == 0 and naromatic == 3 and q == 1 :
						p = 4.10 # :[n+](:):
					elif nh == 0 and nsingle == 1 and naromatic == 2 and q == 1 :
						p = 3.88 # -:[n+](:):
				elif nv == 4 :
					if nh == 0 and nsingle == 4 and q == 1 :
						p = 0.0 # -[N+](-)(-)-
				if p < 0. : # N with non-standard valency
					p = 30.5 - (nv * 8.2) + (nh + 1.5)
					if p < 0. :
						p = 0.
			elif an == 8 :
				if nv == 1 :
					if nh == 0 and ndouble == 1 and q == 0 :
						p = 17.07 # O=
					elif nh == 1 and nsingle == 1 and q == 0 :
						p = 20.23 # [OH]-
					elif nh == 0 and nsingle == 1 and q == 1 :
						p = 23.06 # [O-]-
				elif nv == 2 :
					if nh == 0 and nsingle == 2 and q == 0 and not isIn3ring :
						p = 9.23 # -O-
					elif nh ==0 and nsingle == 2 and q == 0 and isIn3ring :
						p = 12.53 # -O-r3
					elif nh == 0 and naromatic == 2 and q == 0 :
						p = 13.14 # :o:
				if p < 0. : # O with non-standard valency
					p = 28.5 - (nv * 8.6) + (nh * 1.5)
					if p < 0. :
						p = 0.
			psa += p
		return psa
				
	"""    
	def is_toxic(self) :
	"""
	
	def detectFusedRings(self) :
		"""
		trouve les cycles fusionnes dans une molecule.
		on retourne une liste de listes d'atomes par cycle fusionne.
		"""
		res=[] # liste des cycles fusionnes trouve
		index=[] # indice des cycles deja traites
		for i in range(len(self.molecule.cycles)) : # on parcours la liste des cycles de la molecule
			if i in index : # si on a deja traite ce cycle on passe au suivant
				continue
			else : # sinon 
				x=self.molecule.cycles[i].atoms # on recupere tout les atomes du cycle -> x
				index.append(i)
				j=i+1
				while j<len(self.molecule.cycles) : # on regarde les autres cycles de la molecule
					if j in index :
						j+=1
						continue
					else :
						fr = conc(x,self.molecule.cycles[j].atoms, 1) # on verifie si les 2 cycles consideres on au moins 1 atome en commun.
																# Si oui, on renvoie la liste des atomes des 2 cycles sans les doublons -> fr
						if not fr :
							j+=1
							continue
						else :
							x=fr
							index.append(j)
							j=i+1
				res.append(x) # on ajoute la liste des atomes -> res
		# print res
		return res

	def _atomsToRemoveFromBrique(self,brique):
		return [atom for atom in brique.molecule.atoms \
	                     if atom.oatoms == []]

	def bricks(self) :
		"""
		decomposer la molecule en briques -> liste de cycles + liste de linkers et appendices
		"""
		cycles=[] # liste des cycles de la molecule
		linkApp=[] # liste des linkers + appendices de la molecule
		bond = [] # num des atoms des cycles ayant linker ou appendice
		
		### LES CYCLES
		brique = copy.deepcopy(self)
		fr = brique.detectFusedRings() # on recupere une liste des cycles fusionnes de la molecule
		# print "Detected %d rings" % len(fr)
		for i in range(len(fr)) : # on parcours la liste des cycles fusionnes
			for atom in fr[i] : # pr chacun des atomes de ces cycles
				for v in atom.oatoms : # on regarde leurs voisins
					if v not in fr[i] : # si un des voisins n'est pas ds le cycle que l'on considere
						if atom.handle not in bond : # ca veut dire que l'atome du cycle que l'on considere est lie a un linker/appendice, si ca n'est pas un H ...
							for oatom in atom.oatoms:
								if (oatom not in fr[i]) and not (oatom.symbol == 'H'):
									bond.append(atom.handle) # on recupere son 'handle' (de l'atome)
									break
							#bond.append(atom.handle) # on recupere son 'handle' (de l'atome)
		#print fr[i]
			brique.molecule.pruneToAtoms(fr[i]) # on isole le supercycle du reste de la molecule
			# la Brique que l'on vient d'isoler possede tj ds ses attributs la liste de tout les cycles de la molecule de depart
			# ca va causer des pb pr la suite, il faut donc laisser ds cette liste seulement les cycles du supercycle que l'on vient de traiter
			# et retirer tout les autres
			z = 0
			index = []
			while z < len(brique.molecule.cycles) :
				for atom in brique.molecule.cycles[z].atoms:
					if atom.index in index :
						brique.molecule.cycles.remove(brique.molecule.cycles[z])
						index=[]
						break
					else :
						index.append(atom.index)
				else :
					index = []
					z += 1
			#print "\n\n len(cycles) ", len(brique.molecule.cycles), "\n\n"
			cycles.append(Brique(brique.molecule,bond)) # creation d'une instance de Brique -> ajoute a cycles
			# on met a jour les variables et c reparti pr un tour
			brique = copy.deepcopy(self)
			fr = brique.detectFusedRings()
			bond = []
			
		### LES LINKERS - APPENDICES
		# pre-traitement de la molecule -> enlever les cycles		
		for cycle in fr : # pr chaque supercycle 
			# on enleve les liaisons entre les atomes du cycle
			bondsToDel = []
			for atom in cycle:
				for bond in atom.bonds:
					if bond.xatom(atom) in cycle: 
						bondsToDel.append(bond)
			i = 0
			for bond in bondsToDel:
				if bond not in bondsToDel[:i]:
					brique.molecule.remove_bond(bond)
				i += 1
		# on enleve tous les atomes du cycle
		atomsToRemove = self._atomsToRemoveFromBrique(brique)
		for atom in atomsToRemove:
			brique.molecule.remove_atom(atom)
			
		# recherche des linkers + appendices
		toIsolate = [] # liste des atomes a isoler
		la2 = []
		linker = copy.deepcopy(brique)
		la = Components.components(linker.molecule) # composants connectes de la molecule sous la forme [(atoms, bonds), (atoms, bonds)] -> la
		                                            # en fait on a une liste de tuple ou chaque tuple correspond
							    # a une liste des atomes et une liste des liaisons d'un linker/appendice de la molecule
		# print "la", la
		ii = 0
		while ii < len(la):
			atCycle = 0
			for at in la[ii][0]:
				if (at.symbol == 'H'):
				        #if (at.symbol == 'H') and not (at.oatoms[0].rings):
					continue
				if at.rings:
					atCycle += 1
					if atCycle == 1:
						continue
				break
			else:
				pass
			ii += 1
			
		for j in range(len(la)):
			la2.append(la[j][0])
		listIndex = sortIndex(la2) # creation d'une liste d'indice des atomes des linkers et appendices
		for i in range(len(la)):
			for index in listIndex[i] :
				toIsolate.append(linker.molecule.atoms[index]) # pr chaque tuple, on recupere la liste des atomes -> toIsolate
			linker.molecule.pruneToAtoms(toIsolate) # on isole ces atomes de la molecule
			linker.molecule.cycles = [] # le fragment isole n'est pas un cycle donc il faut enlever les cycles de l'attribut 'cycles'
			nbH = 0
			for atom in linker.molecule.atoms:
				if atom.symbol == 'H':
					nbH += 1
			if len(linker.molecule.atoms) - nbH > 1:
				linkApp.append(Brique(linker.molecule,[])) # creation d'une instance de Brique -> ajoute a linkApp
			toIsolate = []
			linker = copy.deepcopy(brique)
		
		return [cycles,linkApp]


	def graph(self, verbose = 0):
                """
                construction du graphe des briques
		return a graph: a list of connections between one ring and a linker/appendix
		(cycle (Brique instance), ring atom index, linker (Brique instance), linker atom index)
		The atom is shared by both ring and linker (i.e. is the ring atom to which linker is connected)
                """
                cycles, linkApp = self.bricks()  # Une liste de cycles fusionnes et de linkers 
		graph = []
		if not linkApp :
			return [cycles]
		if not cycles :
			return [linkApp]
	
		appToCycles = []
		for i in range(len(linkApp)):
			appToCycles.append([])

		# For all linkers / appendices
		# (on parcourt la liste des linker obtenus grace a la methode 'bricks()' )
		for i in range(len(linkApp)): 

			# Look for atoms belonging to a ring
			# Pour ce linker, on liste les atomes qui appartiennent a un cycle.
			# on cree un dictionnaire a base des "handle"
			num = {}
			for atom in linkApp[i].molecule.atoms: # pr chaque atome de cette brique
				if atom.rings: # c'est par cet atome que le linker/appendice est lie a un cycle
					num[atom.handle] = atom.index # on ajoute ds un dictionnaire l'index de cet atom que l'on retrouve grace a son num unique 'handle'
			if verbose:
				sys.stderr.write("Linker %d %s connected to %d rings\n" % (i, linkApp[i].molecule.cansmiles(), len(num)) )

			# For atoms belonging to a ring, look for ring
			traites = 0
			for j in range(len(cycles)):
				for atom in cycles[j].molecule.atoms:
					if atom.handle in num.keys():
						# store  information about cycle, index of atom belonging both to ring and appendix, twice ??)
						if verbose:
							sys.stderr.write("linker %d (%s) cycle %d (%s) linker atom is %d %s %d %d\n" % (i, linkApp[i].molecule.cansmiles(), j, cycles[j].molecule.cansmiles(), atom.index, atom.symbol, atom.handle, num[atom.handle]) )
						appToCycles[i].append((j, atom.index, num[atom.handle])) # (Num cycle, atom index in cycle, atom index in linker) 
						traites += 1
						break
				if traites == len(num):
					break

		# Now we have, for each linker, a list of cycles connected
                cyclesVisites = []
		if verbose:
			sys.stderr.write("appToCycles %d : %s\n" % (0, str(appToCycles[0])))
                for app_cycle in appToCycles[0]:
			# graph is a list of (cycle (Brique), ring atom index, linker (Brique), linker atom index)
			# sys.stderr.write("Cycle %s atom %s %d linked to linker %d atom %s %d\n" % (cycles[app_cycle[0]].cansmiles(), cycles[app_cycle[0]].molecule.atoms[app_cycle[1]].symbol, cycles[app_cycle[0]].molecule.atoms[app_cycle[1]].index, 0, linkApp[0].molecule.atoms[app_cycle[2]].symbol, linkApp[0].molecule.atoms[app_cycle[2]].handle)) 
                        graph.append([cycles[app_cycle[0]], app_cycle[1], linkApp[0], app_cycle[2]])
                        cyclesVisites.append(app_cycle[0])
		
                if len(appToCycles) == 1:
                        return graph
                while True:
                        trouve = False
			done = False
			i = 1
                        while i < len(appToCycles):
				if verbose:
					sys.stderr.write("appToCycles %d : %s\n" % (i, str(appToCycles[i])))
                                for j in range(len(appToCycles[i])):
					app_cycle = appToCycles[i][j]
                                        if appToCycles[i][j][0] in cyclesVisites:
                                                graph.append([cycles[appToCycles[i][j][0]], appToCycles[i][j][1], linkApp[i], appToCycles[i][j][2]])
						if verbose:
							sys.stderr.write("Cycle %s atom %s %d linked to linker %d atom %s %d\n" % (cycles[appToCycles[i][j][0]].cansmiles(), cycles[appToCycles[i][j][0]].molecule.atoms[appToCycles[i][j][1]].symbol, cycles[appToCycles[i][j][0]].molecule.atoms[appToCycles[i][j][1]].index, i, linkApp[i].molecule.atoms[appToCycles[i][j][2]].symbol, linkApp[i].molecule.atoms[appToCycles[i][j][2]].index) )
                                                if len(appToCycles[i]) == 1:
							done = True
						del(appToCycles[i][j])
                                                break
                                else:
                                        i += 1
                                        continue
				if not done:
					# sys.stderr.write("Not done ...\n")
                                	while len(appToCycles[i]) > 0:
                                        	if appToCycles[i][0][0] not in cyclesVisites:
                                        	        cyclesVisites.append(appToCycles[i][0][0])
                                        	graph.append([cycles[appToCycles[i][0][0]], appToCycles[i][0][1], linkApp[i], appToCycles[i][0][2]])
                                        	if len(appToCycles[i]) == 1:
							del(appToCycles[i][0])
							break
						else:
							del(appToCycles[i][0])
                                trouve = True
                                break
                        if not trouve:
                                break
                                  
		if verbose:
			for cnx in graph:
				sys.stderr.write("%s - %s-%d linked to %s - %s-%d\n" % (cnx[0].molecule.cansmiles(), cnx[0].molecule.atoms[cnx[1]].symbol, cnx[0].molecule.atoms[cnx[1]].index, cnx[2].molecule.cansmiles(),  cnx[2].molecule.atoms[cnx[3]].symbol, cnx[2].molecule.atoms[cnx[3]].index))
			# sys.exit(0)
		return graph                                                   



	def graph2molecule(self, workingPath, verbose = 0):
		"""
		Reconstruire la molecule a partir du graphe des briques + affectation de coordonnees a chaque atomes

		Genere une conformation 3D avec les valeurs de diedres canoniques.
		data: a brickgraph
		"""
		mol = copy.deepcopy(self)
		data = mol.brickgraph  # result of self.graph()

		if verbose:
			sys.stderr.write("graph2molecule: ...\n")

		#if there is no cycle in the molecule
		if (len(data[0]) == 1) and (not data[0][0].molecule.cycles):
			# calcul des coordonnees des atomes
			if len(data[0][0].molecule.atoms) >= 3:
				if len(data[0][0].molecule.atoms[0].oatoms) >= 2:
					i, j, k = 0, 1, data[0][0].molecule.atoms[0].oatoms[1].index
					if k == 1:
						j = data[0][0].molecule.atoms[0].oatoms[0].index
				elif len(data[0][0].molecule.atoms[1].oatoms) >= 2:
					i, j, k = 1, 0, 2
			else:
				data[0][0].molecule.atoms[0].x = data[0][0].molecule.atoms[0].y = data[0][0].molecule.atoms[0].z = 1.0
				if len(data[0][0].molecule.atoms) == 2:
					data[0][0].molecule.atoms[1].y = data[0][0].molecule.atoms[1].z = 1.0
					data[0][0].molecule.atoms[1].x = dist(data[0][0].molecule.atoms[0], data[0][0].molecule.atoms[1])
				return
				
			data[0][0].molecule.atoms[i].x = data[0][0].molecule.atoms[i].y = data[0][0].molecule.atoms[i].z = 1.0
			
			data[0][0].molecule.atoms[j].x = 1.0 + dist(data[0][0].molecule.atoms[i], data[0][0].molecule.atoms[j])
			data[0][0].molecule.atoms[j].y = data[0][0].molecule.atoms[j].z = 1.0
			
			dist_ik = dist(data[0][0].molecule.atoms[i], data[0][0].molecule.atoms[k])
			angle_jik = angle(data[0][0].molecule.atoms[j], data[0][0].molecule.atoms[i], \
		                      data[0][0].molecule.atoms[k])
			data[0][0].molecule.atoms[k].x = 1.0 + math.cos(angle_jik)*dist_ik
			data[0][0].molecule.atoms[k].y = 1.0 + math.sin(angle_jik)*dist_ik
			data[0][0].molecule.atoms[k].z = 1.0
						
			if len(data[0][0].molecule.atoms) > 3:
				for index in [i,j,k]:
					setLinkAppCoord(0,data,data[0][0].molecule.atoms[index],self, saveindex = self.saveindex, verbose = verbose)
			
			for atom in data[0][0].molecule.atoms:
				at = self.molecule.atoms[getIndexAtom(getAtom(atom.handle,self.molecule),self.molecule)]
				at.x = atom.x
				at.y = atom.y
				at.z = atom.z
						
			if len(data[0][0].molecule.atoms) > 3:
				self.diedres = self.getDiedres() # recherche des liaisons pouvant tourner
				self.defaultDihedrals()
			self.graphCreated = True
			return
		
		#if the molecule is just a cycle
		elif len(data[0]) == 1:
			if verbose:
				sys.stderr.write("graph2molecule: Just a cyle\n")
			try:
				data[0][0].map(workingPath)
			except CycleUnknown, errMsg:
				#print "\nERROR: " + str(errMsg)
				raise CycleUnknown, errMsg
			for atom in data[0][0].molecule.atoms:
				at = self.molecule.atoms[getIndexAtom(getAtom(atom.handle,self.molecule),self.molecule)]
				at.x = atom.x
				at.y = atom.y
				at.z = atom.z
			for cycle in self.molecule.cycles:
			        for atom in cycle.atoms:
				        v = []
				        nb = []
		          		for oatom in atom.oatoms:
				        	if oatom.symbol == 'H':
					                nb.append(oatom)
					        else:
						        v.append(oatom)
				        else:
					        crds = setCoord(atom,v,nb, saveindex = self.saveindex)
					        for i in range(len(nb)):
						        nb[i].x = crds[i][0]
						        nb[i].y = crds[i][1]
						        nb[i].z = crds[i][2]
			self.graphCreated = True
			return
			
		if verbose:
			sys.stderr.write("graph2molecule: will manage rings ...\n")

		# si la molecule contient au-moins un cycle
		if verbose:
			sys.stderr.write("graph2molecule: at least one cyle\n")
		# data is the graph() result: list of connections cycle/linker
		for i in range(len(data)):
			if verbose:
				sys.stderr.write("graph2molecule: ===> Considering %s versus %s\n" % (data[i][0].molecule.cansmiles(), data[i][2].molecule.cansmiles()))
			#print "for i in data"
			### CYCLE
			if data[i][0].molecule.atoms[0].x == None : # mapping du cycle si c pas encore fait
				if verbose:
					sys.stderr.write("graph2molecule: will manage a ring ...\n")
				try:
					if verbose:
						sys.stderr.write("graph2molecule: will map cycle %s\n" % workingPath)
					# get coordinates for the ring (from library)
					data[i][0].map(workingPath, verbose = verbose)
					# sys.stderr.write("Did map cycle %s\n" % workingPath)					
				except CycleUnknown, errMsg:
					#print "\nERROR: " + str(errMsg)
					if verbose:
						sys.stderr.write("graph2molecule: missing ring ...\n")

					raise CycleUnknown, errMsg
				except:
					pass
					
					
			### LINKER - APPENDICE
			for bond in data[i][2].molecule.bonds : # pr chaque liaison du linker
				if verbose:
					sys.stderr.write("graph2molecule: will manage an appendix for %s ...\n" % data[i][2].molecule.cansmiles())

				for la in bond.atoms : # pr chaque atome de la liaison du linker
					if la == data[i][2].molecule.atoms[data[i][3]] : # si c'est l'atome qui lie la brique a un cycle
						# mise a jour des voisins et des liaisons de l'atome du cycle portant la brique
						# We add to the ring atom the atom of the linker
						data[i][0].molecule.atoms[data[i][1]].oatoms.append(bond.xatom(la)) 
						data[i][0].molecule.atoms[data[i][1]].bonds.append(bond)
						# mise a jour des voisins du voisins de l'atome (!!!) liant la brique au cycle
						# The linker is not connected to ring any longer (avoid several identical connections to ring)
						bond.xatom(la).oatoms.remove(la)
						bond.xatom(la).oatoms.append(data[i][0].molecule.atoms[data[i][1]])
						# mise a jour des atomes de la liaison
						bond.atoms.remove(la)
						bond.atoms.append(data[i][0].molecule.atoms[data[i][1]])
						break
			
			# From here, the ring knows about the linker atoms.
			# 
			# We look if linker coordinates are set.
			linkAppCreated = True
			for at in data[i][2].molecule.atoms:
				if at.x == None:
					linkAppCreated = False
					break
			if linkAppCreated: #if data[i][2].molecule.atoms[data[i][3]].x != None :  bouger coordonnees du cycle pr le mettre au bon endroit si le linker est deja construit
				# coordonnees sur lesquelles on va superposer les coordonnees du cycle a deplacer, pr creer la matrice de transformation
				if verbose:
					sys.stderr.write("graph2molecule: Appendix coordinates are there ...\n")
				qxyz = [[data[i][2].molecule.atoms[data[i][3]].x, data[i][2].molecule.atoms[data[i][3]].y, \
				         data[i][2].molecule.atoms[data[i][3]].z]]
				r = []
				l = []
				for atom in data[i][0].molecule.atoms[data[i][1]].oatoms:
					if atom.rings:
						if atom.findbond(data[i][0].molecule.atoms[data[i][1]]).rings:
							r.append(atom)
						else:
							l.append(atom)
					else:
						l.append(atom)
				if (data[i][0].molecule.atoms[data[i][1]].chiral_class == 2) and (len(r) > 1):
					r[0], r[1] = r[1], r[0]
				nbOfGhostsToPlace = len(self.molecule.atoms[getIndexAtom(data[i][2].molecule.atoms[data[i][3]], \
											                self.molecule)].oatoms) - len(r) - len(l)
				ori = self.molecule.atoms[getIndexAtom(data[i][2].molecule.atoms[data[i][3]], self.molecule)]
				if (ori.number == 7) and (len(ori.oatoms) == 3) and not (ori.aromatic):
					for oat in ori.oatoms:
						if oat.symbol == 'C' and len(oat.oatoms) == 3:
							nbOfGhostsToPlace -= 1
						elif oat.symbol == 'C' and oat.aromatic:
							nbOfGhostsToPlace -= 1
					nbOfGhostsToPlace += 1

				if (ori.number == 16) and (len(ori.oatoms) == 3) and not (ori.aromatic):
					nbOfGhostsToPlace += 1
					
				if nbOfGhostsToPlace < 0:
					nbOfGhostsToPlace = 0
				# sys.stderr.write("LE BON ??\n")
				# sys.stderr.write("Handle is %d\n" % data[i][2].molecule.atoms[data[i][3]].handle)
				sens = 1
				if l != []:
					# sys.stderr.write("Other Handle is %d\n" % l[0].handle)
					for abond in data[i][2].molecule.atoms[data[i][3]].bonds:
						# sys.stderr.write("%s %s\n" % (abond.atoms[0].symbol, abond.atoms[1].symbol))
						# sys.stderr.write("other atom is %s\n" % str(abond.xatom(data[i][2].molecule.atoms[data[i][3]]).handle) )
						if abond.xatom(data[i][2].molecule.atoms[data[i][3]]) == l[0] :
							if abond.dbo:
								# sys.stderr.write("DBO: %s\n" % str(abond.dbo))
								# sys.stderr.write("DBO: %s\n" % str(abond.dbo[0][2]))
								sens = abond.dbo[0][2]
								if sens == 2:
									sens = -1
				crds = setCoord(data[i][2].molecule.atoms[data[i][3]],l,r,nbOfGhostsToPlace, sens, saveindex = self.saveindex, verbose = verbose)
				# sys.stderr.write("APRES LE BON ??\n")
				#print len(self.molecule.atoms[getIndexAtom(data[i][2].molecule.atoms[data[i][3]], self.molecule)].oatoms), len(r), len(l)
				for z in range(len(r)) :
					qxyz.append([crds[z][0],crds[z][1],crds[z][2]])
				
				# coordonnees du cycle a deplacer
				txyz = [[data[i][0].molecule.atoms[data[i][1]].x, data[i][0].molecule.atoms[data[i][1]].y, \
					     data[i][0].molecule.atoms[data[i][1]].z]]
				for atom in data[i][0].molecule.atoms[data[i][1]].oatoms :
					if atom.rings :
						if atom.findbond(data[i][0].molecule.atoms[data[i][1]]).rings :
							txyz.append([atom.x,atom.y,atom.z])
				# creation de la matrice de transformation
				bf = bestFit(3,qxyz,txyz)
				bf.bestFit()
				# deplace de cycle a la bonne place grace a la matrice de transformation
				TM2(data[i][0].molecule.atoms,bf.M1,bf.M2,bf.M3,bf.M4)
				#print "before:"
				#for at in data[i][0].molecule.atoms:
				#	print at.x, at.y, at.z
				#print [at.symbol for at in data[i][2].molecule.atoms]
				#print [ll.symbol for ll in l]
				#print [rr.symbol for rr in r]
				#print crds
				#check for axial/equatorial position
				if self.axialEquatorial:
					
					if verbose:
						sys.stderr.write("graph2molecule: axialEquatorial consideration\n")

					toCheck = False
					crdsL = []
					for atom in l:
						if atom.symbol != 'H':
							descript = self.axialEquatorial[getIndexAtom(atom,self.molecule)]
							if descript in ['A','E']:
								crdsL.insert(0,[atom.x,atom.y,atom.z])
								#crdsL.insert([atom.x,atom.y,atom.z])
								toCheck = True
							else:
								crdsL.append([atom.x,atom.y,atom.z])
						else:
							crdsL.append([atom.x,atom.y,atom.z])
					#print crdsL
					#if (len(crdsL) == 1) and (not self.molecule.atoms[getIndexAtom(data[i][2].molecule.atoms[data[i][3]], \
					#						                self.molecule)].aromatic):
					#	crdsL.append([crds[1][0], crds[1][1], crds[1][2]]) 
					
					#if (len(crdsL) == 2) and (toCheck) and (not self.molecule.atoms[getIndexAtom(data[i][0].molecule.atoms[data[i][1]], self.molecule)].aromatic) and not (isAxEqValidPosition(data[i][2].molecule.atoms[data[i][3]], data[i][0].molecule.atoms,crdsL,descript)):
					if (len(crdsL) == 2) and (toCheck) and not (isAxEqValidPosition(data[i][2].molecule.atoms[data[i][3]], data[i][0].molecule.atoms,crdsL,descript)):
						#print len(crdsL)
						qxyz[1], qxyz[2] = qxyz[2], qxyz[1]
						# coordonnees du cycle a deplacer
						txyz = [[data[i][0].molecule.atoms[data[i][1]].x, data[i][0].molecule.atoms[data[i][1]].y, \
						         data[i][0].molecule.atoms[data[i][1]].z]]
						for atom in data[i][0].molecule.atoms[data[i][1]].oatoms :
							if atom.rings :
								if atom.findbond(data[i][0].molecule.atoms[data[i][1]]).rings :
									txyz.append([atom.x,atom.y,atom.z])
						# creation de la matrice de transformation
						bf = bestFit(3,qxyz,txyz)
						bf.bestFit()
						# deplace de cycle a la bonne place grace a la matrice de transformation
						TM2(data[i][0].molecule.atoms,bf.M1,bf.M2,bf.M3,bf.M4)

				#check for Z/E position
				for bond in data[i][0].molecule.atoms[data[i][1]].bonds:
					if (bond.dbo) and (bond not in data[i][0].molecule.bonds):
						#print "ori:", ori.x, ori.y, ori.z
						oatom = data[i][2].molecule.atoms[data[i][3]]
						dbo = bond.dbo[0]
						
						oatom = data[i][0].molecule.atoms[data[i][1]]
						bondToPlace = 1
						if dbo[0].atoms[0].handle in [at.handle for at in data[i][0].molecule.atoms]:
							bondToPlace = 0
						if dbo[(bondToPlace+1)%2].atoms[0].handle == data[i][2].molecule.atoms[data[i][3]].handle:
							handlePt1 = dbo[(bondToPlace+1)%2].atoms[1].handle
						else:
							handlePt1 = dbo[(bondToPlace+1)%2].atoms[0].handle
						pt1 = None
						for at in data[i][2].molecule.atoms:
							if at.handle == handlePt1:
								pt1 = at
						if not pt1:
							continue
						if dbo[bondToPlace].atoms[0].handle == data[i][0].molecule.atoms[data[i][1]].handle:
							handlePt4 = dbo[bondToPlace].atoms[1].handle
						else:
							handlePt4 = dbo[bondToPlace].atoms[0].handle
						for at in data[i][0].molecule.atoms:
							if at.handle == handlePt4:
								pt4 = at
						v12 = [oatom.x - pt1.x, oatom.y - pt1.y, oatom.z - pt1.z]
						index4 = r.index(pt4)
						#print "cycle"
						#for at in data[i][0].molecule.atoms:
						#	print at, at.x
						#print "app"
						#for at in data[i][2].molecule.atoms:
						#	print at, at.x
						#try:
						#	v34 = [crds[index4][0] - ori.x, crds[index4][1] - ori.y, crds[index4][2] - ori.z]
						#except:
						ori2 = data[i][2].molecule.atoms[data[i][3]]
						#ori2 = data[i][0].molecule.atoms[data[i][1]]
						v34 = [crds[index4][0] - ori2.x, crds[index4][1] - ori2.y, crds[index4][2] - ori2.z]
						scalar = v12[0]*v34[0] + v12[1]*v34[1] + v12[2]*v34[2]
						#print "scalar:", scalar, dbo
						
						#if not (((dbo[2] == 2) and (scalar > 0)) or ((dbo[2] == 1) and (scalar < 0))):
						if not (((dbo[2] == 2) and (scalar < 0)) or ((dbo[2] == 1) and (scalar > 0))):
							#print "before:"
							#for at in data[i][0].molecule.atoms:
							#	print at.x, at.y, at.z
						        #if not (((dbo[2] == 2) and (scalar < 0)) or ((dbo[2] == 1) and (scalar > 0))):
							qxyz[1], qxyz[2] = qxyz[2], qxyz[1]
							# coordonnees du cycle a deplacer
							txyz = [[data[i][0].molecule.atoms[data[i][1]].x, data[i][0].molecule.atoms[data[i][1]].y, \
							         data[i][0].molecule.atoms[data[i][1]].z]]
							for atom in data[i][0].molecule.atoms[data[i][1]].oatoms:
								if atom.rings :
									if atom.findbond(data[i][0].molecule.atoms[data[i][1]]).rings:
										txyz.append([atom.x,atom.y,atom.z])
							# creation de la matrice de transformation
							bf = bestFit(3,qxyz,txyz)
							bf.bestFit()
							# deplace de cycle a la bonne place grace a la matrice de transformation
							TM2(data[i][0].molecule.atoms,bf.M1,bf.M2,bf.M3,bf.M4)
							#print "after ZE mvt:"
							#for at in data[i][0].molecule.atoms:
							#	print at.x, at.y, at.z
			
			#if the linker/appendice has not yet been created
			else:
				ori = data[i][2].molecule.atoms[data[i][3]] 
				if verbose:
					sys.stderr.write("graph2molecule: linkAppCreated false\n")
					sys.stderr.write("graph2molecule: ori %s\n" % ori)
					sys.stderr.write("Chiral sense. ori %s atom type %s\n" % (str(ori.chiral_class), str(ori.number)))
				# ith connection in data. ori is linker atom instance 
				setLinkAppCoord(i,data,ori,self, saveindex = self.saveindex, verbose = verbose)
				data[i][2].molecule.atoms[data[i][3]].x = data[i][0].molecule.atoms[data[i][1]].x
				data[i][2].molecule.atoms[data[i][3]].y = data[i][0].molecule.atoms[data[i][1]].y
				data[i][2].molecule.atoms[data[i][3]].z = data[i][0].molecule.atoms[data[i][1]].z

				
		#setting of the atomic coordinates of the molecule from the bricks 
		if verbose:
			sys.stderr.write("graph2molecule: assembling bricks ...\n")
		handleOfVisited = []
		for i in range(len(data)):
			for j in [0,2]:
				for atom in data[i][j].molecule.atoms:
					if atom.handle not in handleOfVisited:
						at = self.molecule.atoms[getIndexAtom(getAtom(atom.handle,self.molecule),self.molecule)]
						at.x = atom.x
						at.y = atom.y
						at.z = atom.z
						handleOfVisited.append(atom.handle)
						
		if verbose:
                    sys.stderr.write("graph2molecule: will manage ring hydrogens ...\n")

		#set the coordinates of hydrogens attached to cycles
		for cycle in self.molecule.cycles:
			for atom in cycle.atoms:
				v = []
				nb = []
				for oatom in atom.oatoms:
					for atCycle in atom.rings:
                                                if oatom in atCycle.atoms:
                                                        v.append(oatom)
                                                        break
                                        else:
                                                if (oatom.symbol == 'H') and (oatom.x == None):
                                                        nb.append(oatom)
                                                else:
                                                        break
				else:
					if verbose:
						sys.stderr.write("graph2molecule: atom %s v %s nb %s\nWill setCoord\n" % (str(atom),str(v),str(nb)))
						# print "atom:", atom, "v:",v,"nb:",nb
					crds = setCoord(atom,v,nb, saveindex = self.saveindex, verbose=verbose)
					if verbose:
						sys.stderr.write("graph2molecule: setCoord done ...\n")
					for i in range(len(nb)):
						nb[i].x = crds[i][0]
						nb[i].y = crds[i][1]
						nb[i].z = crds[i][2]
		#research of rotateable bond
		if verbose:
			sys.stderr.write("graph2molecule: bricks assembly done ...\n")
		self.diedres = self.getDiedres(verbose = verbose)
		if verbose:
			sys.stderr.write("graph2molecule: installing conformer default dihedral values...\n")
                self.defaultDihedrals()
		self.graphCreated = True
		if verbose:
			sys.stderr.write("graph2molecule: completed ...\n")
			
	def writeToSdf(self,file="mol",fournisseur = "") :
		"""
		ecrire une molecule en sdf
		"""
		mod = "w"
		if os.path.isfile(file) :
			mod = "a"
		### HEADER
		#result = [self.molecule.cansmiles()]
		result = [self.smiles] #thiaggggggggggggggggg
		result.append("  -ISIS-")
		result.append("")
		
		numatoms = len(self.molecule.atoms)
		numbonds = len(self.molecule.bonds)
		result.append("%3i%3i  0  0  0  0  0  0  0  0  1"%(numatoms, numbonds))
		### ATOME BLOC
		for atom in self.molecule.atoms :
			x = atom.x
			if not x : x = 0.
			y = atom.y
			if not y : y = 0.
			z = atom.z
			if not z : z = 0.
			result.append("%10.4f%10.4f%10.4f%3s%3i%3i%3i%3i%3i"%(x, y, z, atom.symbol, 0, 0, 0, 0, 0))
		### TABLE DE CONNECTIVITE
		for bond in self.molecule.bonds :
			result.append("%3i%3i%3i%3i%3i%3i"%(bond.atoms[0].index+1, bond.atoms[1].index+1, bond.bondtype,0, 0, 0))
		### FIN DU FICHIER
		result.append("M  END")
		result.append(">  <ID>")
		result.append(self.id)
		result.append("")
		if fournisseur :
			result.append(">  <FOURNISSEUR>")
			result.append(fournisseur)
			result.append("")
		result.append("$$$$")
		
		f = open(file,mod)
		f.write("\n".join(result))
		f.write("\n")
		f.close()

	def writeToPdb(self,file="mol",withH=True):
		"""
		ecrire une molecule en pdb
		"""
		mod = "w"
		if os.path.isfile(file):
			mod = "a"
		result = ["COMPND    "+self.smiles]
		result.append("AUTHOR    RPBS-TEAM")
		### ATOME BLOC
		for atom in self.molecule.atoms:
			if (not withH) and (atom.symbol == 'H'):
				continue
			x = atom.x
			if not x : x = 0.
			y = atom.y
			if not y : y = 0.
			z = atom.z
			if not z : z = 0.
			result.append("ATOM"+str(atom.index+1).rjust(7)+"  "+atom.symbol.upper().ljust(2)+"  UNK     0"+"%12.4f%8.4f%8.4f"%(x,y,z)+"  1.00  0.00"+atom.symbol.upper().rjust(12))
		### TABLE DE CONNECTIVITE
		for atom in self.molecule.atoms:
			if (not withH) and (atom.symbol == 'H'):
				continue
			tmp = "CONECT"+str(atom.index+1).rjust(5)
			for oatom in atom.oatoms:
				if (not withH) and (oatom.symbol == 'H'):
					continue
				tmp += str(oatom.index+1).rjust(5)
			result.append(tmp)
		result.append("END")
		
		f = open(file,mod)
		f.write("\n".join(result))
		f.write("\n")
		f.close()

	
	def getDiedres(self, verbose = 0) :
		"""
		renvoi listes des diedres pouvant tourner
		"""
		res = [] # liste des angles diedres pouvant tourner
		toExclude = []

		if verbose: 
			sys.stderr.write("getDiedres:\n")
		fRings = self.detectFusedRings() # on recupere une liste des cycles fusionnes de la molecule
		if verbose:
			for fr in fRings:
				for atm in fr:
					sys.stderr.write("%s %d " % (atm.symbol, atm.index))
				sys.stderr.write("\n")
		
		# Scan compound bonds.
		for bond in self.molecule.bonds: 
			if verbose: 
				sys.stderr.write("getDiedres: considering bond %s %d - %s %d (contacts %d rings)\n" % (bond.atoms[0].symbol, bond.atoms[0].index, bond.atoms[1].symbol, bond.atoms[1].index, len(bond.rings)))
				
			# We do not accept bonds involving hydrogens.
			if (bond.atoms[0].symbol == "H") or (bond.atoms[1].symbol == "H"):
				if verbose:
					sys.stderr.write("getDiedres: flush since hydrogens\n")
				continue


			# We do not accept bonds in a ring.
			inRing = 0
			for fr in fRings:
				if (bond.atoms[0] in fr) and (bond.atoms[1] in fr):
					if verbose:
						sys.stderr.write("getDiedres: flush since in same ring\n")
					inRing = 1
					break
					
					
			# si la liaison n'est pas ds un cycle et que c'est une liaison simple
			# if (not bond.rings) and (bond.bondtype == 1) and (bond not in toExclude): 
			if verbose:
				sys.stderr.write("getDiedres: inRing %s bondtype %d bond in toExclude %s\n" % (str(inRing), bond.bondtype, str(bond in toExclude)) )
			if (not inRing) and (bond.bondtype == 1) and (bond not in toExclude): 
				if verbose: 
					sys.stderr.write("getDiedres: considering further bond %s %d - %s %d\n" % (bond.atoms[0].symbol, bond.atoms[0].index, bond.atoms[1].symbol, bond.atoms[1].index))
				# si un des atomes n'a qu'un seul voisin, il n'y a pas assez d'atomes pr definir un angle diedre -> suivant
				if (len(bond.atoms[0].oatoms) == 1) or (len(bond.atoms[1].oatoms) == 1):
					if verbose: 
						sys.stderr.write("getDiedres: only one neighbour for one atom\n")
					continue
				
				for i in range(2):
					#if len(bond.atoms[i].oatoms) == 1:
					#	break
					nbH = 0
					for atom in bond.atoms[i].oatoms:
						if atom.symbol == 'H':
							nbH += 1
					# Other neighbours are hydrogens: we pass
					if len(bond.atoms[i].oatoms) - nbH <= 1:
						break
				else:
					# OK, we find the list of atoms defining the angle
					atoms = [] # liste des atomes definissant l'angle diedre
					first = True
					ok = True
					for atom in bond.atoms:
						max = 0
						plusLourd = 0
						# We look for the atom 
						for a in atom.oatoms:
							if (a.findbond(atom).bondtype == 3) or \
							   ((a.findbond(atom).bondtype == 2) and (len(atom.oatoms) == 2)):
								if (a.findbond(atom) not in toExclude) and (len(a.oatoms) > 1):
									aatoms = getNextAtoms(atom,a)
									if len(aatoms) == 1:
										aaatoms = getNextAtoms(a,aatoms[0])
										for aaa in aaatoms:
											if aaa.rings:
												plusLourd = aaa
												break
											elif plusLourd:
												if aaa.number > plusLourd:
													plusLourd == aaa
											else:
												plusLourd = aaa
										before = aatoms[0]
									else:
										for aa in aatoms:
											if aa.rings:
												plusLourd = aa
												break
											elif plusLourd:
												if aa.number > plusLourd:
													plusLourd == aa
											else:
												plusLourd = aa
										before = a
									toExclude.append(a.findbond(atom))	
									if first:
										atoms.append(plusLourd)
										atoms.append(before)
										first = False
									else:
										atoms.append(before)
										atoms.append(plusLourd)
								#	break
								#else:
								ok = False
								break
							if a.findbond(atom).bondtype == 2:
								for aa in atom.oatoms:
									if (aa != a) and (aa not in bond.atoms):
										plusLourd = aa
								break
							if a in bond.atoms :
								continue
							if a.rings and len(a.oatoms) > 2 :
								plusLourd = a
								break
							if len(a.oatoms) == 1 :
								if len(atom.oatoms) != 2 :
									if not plusLourd :
										plusLourd = a
									continue
							if a.number > max :
								plusLourd = a
								max = a.number
						if not ok:
							break
						if first:
							atoms.append(plusLourd)
							atoms.append(atom)
							first = False
						else:
							atoms.append(atom)
							atoms.append(plusLourd)
					else:
						atomsAfterBond = getDividedMoleculeAtoms(atoms[1].findbond(atoms[2]))[1]
						if verbose:
							sys.stderr.write("Appending dihedral\n")
						res.append(Diedre(atoms,atomsAfterBond))
		if verbose:
			sys.stderr.write("getDiedres: found %d dihedrals\n" % len(res))
		return res

	def defaultDihedrals(self, setAngle = True, verbose = 0) :
		"""
		Initialiser les valeurs canoniques des diedres qui ont ete identifies.
		@param setAngle : if True, we setup the first canonical angle values for the dihedrals.
		"""
		self.defaultDihedralValues = []
		for d in self.diedres :
			if verbose:
				sys.stderr.write("defaultDihedrals: diedre %d %d %d %d\n" % (d.atoms[0].index,d.atoms[1].index,d.atoms[2].index,d.atoms[3].index))
			vals = getDihedralAngleValues(d.atoms[0],d.atoms[1],d.atoms[2],d.atoms[3])
			if verbose:
				sys.stderr.write("defaultDihedrals: values: %s\n" % str(vals))
			self.defaultDihedralValues.append(vals)
			if setAngle:
				if verbose:
					sys.stderr.write("defaultDihedrals: Assigning angles\n")
				angleRot = (self.defaultDihedralValues[-1][0][0]*math.pi)/180 - d.value
				if verbose:
					sys.stderr.write("angle value: %f\n" % angleRot)
				# calcule de la matrice de transformation
				ar = AxeRot(d.atoms[1],d.atoms[2],angleRot)
				ar.axeRot()
				# on applique la matrice de rotation aux atomes
				# print "moving list:", d.movingList
				if verbose:
					sys.stderr.write("Will TM %s %s %s %s\n" % (ar.M1,ar.M2,ar.M3,ar.M4))
				TM2(d.movingList,ar.M1,ar.M2,ar.M3,ar.M4)
				# on recalcule l'angle diedre
				if verbose:
					sys.stderr.write("Will computeValue\n")
				d.setValue(d.computeValue())			

	def to1confWWW(self, fileName, format="mol2", refFileName = None, energeticBarrer=100, eIni = 300, nbBestResults=20, mcsteps=100, mode = "mono", minimize=False, ammp_energy = False, autoMinimize = False, clusterize = False, clus_trhld = 0.8, visuSdfFileName="/dev/null/", split = False, cleanup = 0, vibrate = 0, miniEach = False, verbose = 0):
		"""
		toconfWWW: the call to the external multiconf generator
		
		@param mode: one of "mono", "multi", to select for monoconf or multiconf.
		"""
		if verbose:
			sys.stderr.write("to1confWWW: \n")
                if format not in ["mol2", "sdf", "pdb"]:
			raise Exception("Unrecognized file format for toNconf(...)")

		# Name for the bad energy compounds
		if fileName.count(".%s" % format):
			badEnergyFileNamePrefix = fileName.replace(".%s" % format,"_badEnergy")
		else:
			badEnergyFileNamePrefix = "%s_badEnergy" % fileName

		try:
			fileNamePrefix = fileName[:fileName.rindex(".")]
		except:
			fileNamePrefix = fileName
			
		filePrefix = self.id

		# """
		# We need an internal representation of the compound
		# """
	        if not self.graphCreated:
			if verbose:
				sys.stderr.write("to1confWWW: Will generate 3D conformation ... \n")
			if fileName.find("/") != -1:
				self.graph2molecule(fileName[0:fileName.rfind("/")+1], verbose = verbose)
			else:
				self.graph2molecule("./", verbose = verbose)

		if verbose:
			sys.stderr.write("to1confWWW: dihedral identification: \n")

		# """
		# We need identified dihedrals
		# """
		if not self.diedres:
			if verbose:
				sys.stderr.write("to1confWWW: no dihedral. Direct output. \n")
	                coords = []
        	        for at in self.molecule.atoms:
                        	coords.append((at.x, at.y, at.z))
			coordsToMol2(self, [coords], ["no computed"], "%s.ori.mol2" % filePrefix)
			# PT 13/01/2011
			if BABEL_VERSION < 220:
				babel_reprotonate_mol2_old(filePrefix,".ori.mol2", ".mol2", verbose=verbose)
			else:
				babel_reprotonate_mol2(filePrefix,".ori.mol2", ".mol2", verbose=verbose)

			# Now we need to convert into proper format
			if format == "mol2":
				if not split:
					cmd = "cat %s.mol2 >> %s.mol2" % (filePrefix, fileNamePrefix)
					os.system(cmd)
			elif format == "pdb":
				babel_convert("%s.mol2" % filePrefix, "mol2", "%s.pdb" % filePrefix, "pdb", verbose = verbose)
				if not split:
					cmd = "cat %s.pdb >> %s.pdb" % (filePrefix, fileNamePrefix)
					os.system(cmd)
				else:
					# os.system("rm %s.mol2" % filePrefix)
					pass
			else:
				babel_convert("%s.mol2" % filePrefix, "mol2", "%s.sdf" % filePrefix, "sdf", verbose = verbose)
				if not split:
					cmd = "cat %s.sdf >> %s.sdf" % (filePrefix, fileNamePrefix)
					os.system(cmd)
				else:
					# os.system("rm %s.mol2" % fileNamePrefix)
					pass
		# babel_convert("%s.mol2" % fileNamePrefix, "mol2", "%s" % visuSdfFileName, "sdf", verbose = verbose)
			if not split:
				# os.system("rm  %s[.,_]*" % (filePrefix))
				pass
			
			babel_convert("%s.mol2" % filePrefix, "mol2", "%s.sdf" % filePrefix, "sdf", verbose = verbose)
			cmd = "cat %s.sdf >> %s" % (filePrefix, visuSdfFileName)
			os.system(cmd)

                	# if format == "mol2":
                        # 	coordsToMol2(self, [coords], ["no computed"], fileName)
                	# elif format == "pdb":	
                        # 	coordsToPdb(self, [coords], ["no computed"], fileName)
                	# else:
                        # 	coordsToSdf(self, [coords], ["no computed"], fileName)
			# coordsToSdf(self, [coords], ["no computed"], visuSdfFileName)
			self.eOK = True
			return 1
		
		# """
		# We need atomic types (MMFF is OK)
		# """
		if verbose:
			sys.stderr.write("to1confWWW: MMFF Types. \n")
                tmp = str(len(self.molecule.atoms)) + "\n" + str(int(energeticBarrer)) + "\n" + str(mcsteps) + "\n"
                for atom in self.molecule.atoms:
			try:
				vdwParams = VDW_PARAM[str(atom.mmffAtomType)]
			except:
				vdwParams = VDW_PARAM[str(getSimplestMMFFAtomType(atom))]
                        tmp += str(atom.x)+" "+str(atom.y)+" "+str(atom.z)+" "+str(getMMFFPartialAtomicCharge(atom))+" "+ \
                               str(vdwParams[0])+" "+str(vdwParams[1])+" "+str(vdwParams[2])+" "+str(vdwParams[3])+"\n"
		# """
		# We need to know what atoms depend on each dihedral modification
		# """
		if verbose:
			sys.stderr.write("to1confWWW: Atom dihedral dependencies. \n")
                i = 0
		for atom in self.molecule.atoms[:-1]:
                        atomsJoined = [atom]
                        for atJ in atom.oatoms:
                                if (atJ not in atomsJoined) and (atJ.index > i):
                                        atomsJoined.append(atJ)
                                for oat in atJ.oatoms:
                                        if (oat not in atomsJoined) and (oat.index > i):
                                                atomsJoined.append(oat)
                        for ring in atom.rings:
                                for at in ring.atoms:
                                        if (at not in atomsJoined) and (at.index > i):
                                                atomsJoined.append(at)
                        for j in range(i+1, len(self.molecule.atoms)):
                                if j != i+1:
                                        tmp += " "
                                if self.molecule.atoms[j] not in atomsJoined:
                                        tmp += "1"
                                else:
                                        tmp += "0"
                        tmp += "\n"
                        i += 1

		# """
		# We need to know what canonical conformations (number, values) for each dihedral
		# """
		if verbose:
			sys.stderr.write("to1confWWW: dihedral canonical values. \n")
		nbConfsPoss = 1
		pepLikes = []
		for index, diedre in enumerate(self.diedres):
			# We check for peptidic bonds
			# We could check for other such particular types
			pepLike = 0
			# if (diedre.atoms[2].rings == []) and (diedre.atoms[3].rings == []) and \
			#    (len(diedre.atoms[2].oatoms) == 3) and (len(diedre.atoms[3].oatoms) == 3):
			# 	if (diedre.atoms[2].symbol == 'C' and diedre.atoms[3].symbol == 'N') or \
			# 	   (diedre.atoms[2].symbol == 'N' and diedre.atoms[3].symbol == 'C'):
			# 		pepLike = 1
			if (diedre.atoms[1].rings == []) and (diedre.atoms[2].rings == []):
				cat = nat = None
				if ((diedre.atoms[1].symbol == "C") and (diedre.atoms[2].symbol == "N")):
					cat = diedre.atoms[1]
					nat = diedre.atoms[2]
				elif ((diedre.atoms[1].symbol == "N") and (diedre.atoms[2].symbol == "C")):
					cat = diedre.atoms[2]
					nat = diedre.atoms[1]
				if (cat is not None) and (nat is not None):	
					possibly = 0
					for oat in nat.oatoms:
						if oat.symbol == "H":
							possibly = 1
					if possibly:
						for bond in cat.bonds:
							if (bond.xatom(cat).symbol == "O") and (bond.bondtype == 2):
								pepLike = 1
			pepLikes.append(pepLike)
			if pepLike == 1:
				nbConfsPoss *= 2
			else:
				nbConfsPoss *= len(self.defaultDihedralValues[index])

		tmpp = str(len(self.diedres)) + "\n"
		for index, diedre in enumerate(self.diedres):
			tmpp += str(diedre.atoms[0].index)+" "+str(diedre.atoms[1].index)+" "+str(diedre.atoms[2].index)+" "+str(diedre.atoms[3].index)
			if verbose:
				# sys.stderr.write("diedre %d %s-%d-%s %s-%d-%s %s-%d-%s %s-%d-%s peplike\n" % (index, diedre.atoms[0].symbol, diedre.atoms[0].index, str(diedre.atoms[0].rings), diedre.atoms[1].symbol, diedre.atoms[1].index, str(diedre.atoms[1].rings), diedre.atoms[2].symbol, diedre.atoms[2].index, str(diedre.atoms[2].rings), diedre.atoms[3].symbol, diedre.atoms[3].index, str(diedre.atoms[3].rings)))
				sys.stderr.write("diedre %d %s-%d-%d %s-%d-%d %s-%d-%d %s-%d-%d peplike: %d\n" % (index, diedre.atoms[0].symbol, diedre.atoms[0].index, diedre.atoms[0].mmffAtomType, diedre.atoms[1].symbol, diedre.atoms[1].index, diedre.atoms[1].mmffAtomType, diedre.atoms[2].symbol, diedre.atoms[2].mmffAtomType, diedre.atoms[2].index, diedre.atoms[3].symbol, diedre.atoms[3].mmffAtomType, diedre.atoms[3].index,pepLikes[index]))
				sys.stderr.write("%s\n" % str((diedre.atoms[2].rings == []) and (diedre.atoms[3].rings == [])))
			# tmpp += str(diedre.atoms[1].index)+" "+str(diedre.atoms[2].index)
			if pepLikes[index] == 1:
				tmpp += " 2 0.0 0.0 180.0 0.0\n"
			else:
				tmpp += " "+str(len(self.defaultDihedralValues[index]))
				for j in range(len(self.defaultDihedralValues[index])):
					tmpp += " "+str(self.defaultDihedralValues[index][j][0])+" "+str(self.defaultDihedralValues[index][j][1])
				tmpp += "\n"
			mvL = diedre.movingList
			tmpp += str(len(mvL))
			for atom in mvL:
				tmpp += " "+str(atom.index)
			tmpp += "\n" + str(pepLikes[index]) + "\n"
		tmp += tmpp + str(NB_MAX_CONFS_TO_EXPLORE) + "\n"	

		# START This is for debug purpose
		coords = []
		for at in self.molecule.atoms:
			coords.append((at.x, at.y, at.z))
			# print at.x, at.y, at.z
			 
		if verbose:
			sys.stderr.write("to1confWWW: Output of .ori conformation. \n")
		coordsToMol2(self, [coords], ["no energy computed"], filePrefix+".ori.mol2")
		if verbose:
			sys.stderr.write("to1confWWW: Did output of .ori conformation. \n")

                # if format == "mol2":
                #         coordsToMol2(self, [coords], ["no energy computed"], filePrefix+".ori.mol2")
                # elif format == "pdb":
                #         coordsToPdb(self, [coords], ["no energy computed"], filePrefix+".ori.pdb")
                # else:
                #         coordsToSdf(self, [coords], ["no energy computed"], filePrefix+".ori.sdf")
		# END DEBUG ADDENDUM

		# RINGS
		rings = self.molecule.cycles
		tmp = tmp + str(len(rings))+"\n"
		# Output number of rings
		# then per ring: number of atoms, atom indices
		for ring in rings:
			ratoms = ring.atoms
			tmp = tmp + str(len(ratoms))+" "
			idex = []
			for atom in ratoms:
				idex.append(atom.index)
			idex.sort()
			idex = map(str,idex)
			tmp = tmp + " ".join(idex)
			tmp = tmp + "\n"
		# """
		# Now, we generate the conformation
		# """

		# """
		# We need a full mol2 with H
		# """
		if refFileName == None:
			# PT 13/01/2011
			if BABEL_VERSION < 220:
				babel_reprotonate_mol2_old(filePrefix,".ori.mol2", ".mol2", verbose=verbose)
			else:
				babel_reprotonate_mol2(filePrefix,".ori.mol2", ".babel.mol2", verbose=verbose)
		else:
			cmd = "cp %s  %s" % (refFileName, filePrefix+".babel.mol2")
			os.system(cmd)
			if verbose:
				sys.stderr.write("%s\n" % cmd)

		# Here we minimize input for monoconfv3 TO CHECK
		# We cannot minimize here since we need coordinates preserved 
		# to detect equivalence.
		# cmd = "cp %s.babel.mol2 %s.babel.ori.mol2" % (filePrefix, filePrefix)
		# os.system(cmd)
		# AMMOS_Minimize("%s.babel.mol2" % filePrefix)
		# cmd = "cp %s.babel_minimized.mol2 %s.babel.mol2" % (filePrefix, filePrefix)
		# if verbose:
		# 	sys.stderr.write("%s\n" % cmd)
		# # prefer this to cp since solves some format problems
		# babel_convert("%s.babel_minimized.mol2" % filePrefix, "mol2", "%s.babel.mol2" % filePrefix, "mol2", verbose = verbose)

		# os.system(cmd)
		
		# """
		# We also write the internal data
		# """
		# ftmp = open("TestData.data","w")
		try:
			if verbose:
				sys.stderr.write("data is: %s.ori.data\n" % (filePrefix))
			ftmp = open("%s.ori.data" % filePrefix, "w")
			ftmp.write(tmp)
			ftmp.close()
		except:
			if verbose:
				sys.stderr.write("could not write DATA file\n")
			# sys.stderr.write(tmp)
			sys.exit(0)
		# sys.exit(0)

		# """
		# Now we call the generator
		# """
		genMode = "-mono"
		# decide which strategy to adopt if multi
		if mode == "multi":
			# if nbConfsPoss > NB_MAX_CONFS_TO_EXPLORE:
			if nbConfsPoss > nbBestResults:
				genMode = "-rmulti %d" % int(nbBestResults)
#				fi, fo = os.popen2(MULTICONF_EXPLORE_SOME,"t")
#				fi.write(tmp+str(NB_MAX_CONFS_TO_EXPLORE)+"\n")
			else:
				genMode = "-fmulti %d " % int(nbBestResults)
#				fi, fo = os.popen2(MULTICONF_EXPLORE_ALL,"t")
#				fi.write(tmp)
		

		if split:
			oMol2FName = "%s.mol2" % filePrefix
		else:
			oMol2FName = fileNamePrefix
		cmd = "%s %s -idata %s.ori.data -imol2 %s.babel.mol2 -omol2 %s.mol2 -id %s -eIni %f" % (MONOCONFv2, genMode, filePrefix, filePrefix, filePrefix, self.id, eIni)
		if vibrate:
			cmd = "%s -vb " % cmd	

		if miniEach:
			cmd = "%s -mini %s " % (cmd, AMMOSMINIHOME)
		else:
			cmd = "%s -ammos %s " % (cmd, AMMOSMINIHOME)
		if clusterize:
			cmd = "%s -rmsd %s " % (cmd, str(clus_trhld))
		# if verbose:
		# 	cmd = "%s -v " % cmd
		if not verbose:
			cmd = "%s 2> /dev/null" % cmd
		if verbose:
			sys.stderr.write("%s\n" % cmd)
		# fi, fo = os.popen2(cmd,"t")
                # fi.write(tmp)
                # fi.close()
                # lines = fo.readlines()
                # fo.close()
		nConfGen = "0"

		# try:
		p = Popen([cmd], shell=True, # bufsize=bufsize,
			  stdin=PIPE, stdout=PIPE, close_fds=True)
		# (fi, fo) =  (p.stdin, p.stdout)
		# fi.write(tmp)
		# fi.close()

		# lines = []
		# while True:
		# 	line = fo.readline()
		# 	if not line:
		# 		break
		# 	lines.append(line)
		
		# 	sys.stderr.write("%s %s\n" % (self.id, lines[0]))
		lines = p.communicate(tmp)[0] # Will just return the number of confs generated.
		nConfGen      = lines.split()[0]
		bestEnergy    = lines.split()[1]
		bestNbEnergy  = lines.split()[2]
		bestBEnergy   = lines.split()[3]
		quickMinimize = False
		fullMinimize = False
		if autoMinimize:
			if float(bestEnergy) > float(eIni):
				if float(bestBEnergy) < float(eIni):
					quickMinimize = True
				else:
					fullMinimize = True
		if verbose:
			sys.stderr.write("%s %s %s\n" % (self.id, nConfGen, bestEnergy))
		# print "finished"
		p.wait()
		if p.returncode != 0:
			sys.stderr.write("%s conformer generation failed\n" % (self.id))
		# except:
		# 	sys.stderr.write("Failed to generate multiconf\n")

		# coords = []
                # for i in range(len(self.molecule.atoms)):
		# 	 atCoords = lines[i+2].split()
                #          coords.append((float(atCoords[0]), float(atCoords[1]), float(atCoords[2])))
                # mCoords = [coords]
		# nrjs = [lines[1].split()[0]]
		
                # if format == "mol2":
                #         coordsToMol2(self, mCoords, nrjs, fileName)
                # elif format == "pdb":
                #         coordsToPdb(self, mCoords, nrjs, fileName)
                # else:
                #         coordsToSdf(self, mCoords, nrjs, fileName)
		# coordsToSdf(self, mCoords, nrjs, visuSdfFileName)

		# Now we minimize if requested
		if verbose:
			sys.stderr.write( "minimize %s ammp_energy %s\n" % (str(minimize), str(ammp_energy)))
		if minimize or fullMinimize:
			# sys.stderr.write("WILL MINIMIZE\n")
			cmd = "cp %s.mol2 %s.raw.mol2" % (filePrefix, filePrefix)
			if verbose:
				sys.stderr.write("%s\n" % cmd)
			os.system(cmd)
			AMMOS_Minimize("%s.mol2" % filePrefix)
			cmd = "cp %s_minimized.mol2 %s.mol2" % (filePrefix, filePrefix)
			if verbose:
				sys.stderr.write("%s\n" % cmd)
			os.system(cmd)
			try:
				f = open("%s_energy.txt" % filePrefix)
				lines = f.readlines()
				f.close()
				bestEnergy = lines[1].split()[4]
			except:
				pass
		elif quickMinimize:
			# sys.stderr.write("WILL QUICK MINIMIZE\n")
			cmd = "cp %s.mol2 %s.raw.mol2" % (filePrefix, filePrefix)
			if verbose:
				sys.stderr.write("%s\n" % cmd)
			os.system(cmd)
			AMMOS_QuickMinimize("%s.mol2" % filePrefix)
			cmd = "cp %s_minimized.mol2 %s.mol2" % (filePrefix, filePrefix)
			if verbose:
				sys.stderr.write("%s\n" % cmd)
			os.system(cmd)
			try:
				f = open("%s_energy.txt" % filePrefix)
				lines = f.readlines()
				f.close()
				bestEnergy = lines[1].split()[4]
			except:
				pass
		elif ammp_energy:
			sys.stderr.write("AMMP energy calculation\n")
			AMMOS_Energy("%s.mol2" % filePrefix)
			cmd = "cp %s_energy.mol2 %s.mol2" % (filePrefix, filePrefix)
			if verbose:
				sys.stderr.write("%s\n" % cmd)
			os.system(cmd)

		if clusterize:
			cmd = "%s -imol2 %s.mol2 -omol2 %s-cluster.mol2 -rmsd %f" % (CLUSTCMD, filePrefix, filePrefix, clus_trhld)
			if verbose:
				sys.stderr.write("%s\n" % cmd)
			os.system(cmd)
			cmd = "cp %s-cluster.mol2 %s.mol2" % (filePrefix, filePrefix)
			os.system(cmd)

			
		# Now we need to convert into proper format
		if verbose:
			sys.stderr.write("Will check Energy is correct: %lf / %lf\n" % (float(bestEnergy) , eIni))
		if format == "mol2":
			if not split:
				if float(bestEnergy) > eIni:
					cmd = "cat %s.mol2 >> %s.mol2" % (filePrefix, badEnergyFileNamePrefix)
				else:
					cmd = "cat %s.mol2 >> %s.mol2" % (filePrefix, fileNamePrefix)
				os.system(cmd)
                elif format == "pdb":
			babel_convert("%s.mol2" % filePrefix, "mol2", "%s.pdb" % filePrefix, "pdb", verbose = verbose)
			if not split:
				if float(bestEnergy) > eIni:
					cmd = "cat %s.pdb >> %s.pdb" % (filePrefix, badEnergyFileNamePrefix)
				else:
					cmd = "cat %s.pdb >> %s.pdb" % (filePrefix, fileNamePrefix)
				os.system(cmd)
			else:
				# os.system("rm %s.mol2" % filePrefix)
				pass
                else:
                        babel_convert("%s.mol2" % filePrefix, "mol2", "%s.sdf" % filePrefix, "sdf", verbose = verbose)
			# Bugfix  for babel does not add $$$$ on file end if single compound !
			if mode == "mono":
				cmd = "echo '$$$$' >> %s.sdf" % filePrefix
				os.system(cmd)
			if not split:
				if float(bestEnergy) > eIni:
					cmd = "cat %s.sdf >> %s.sdf" % (filePrefix, badEnergyFileNamePrefix)
				else:
					cmd = "cat %s.sdf >> %s.sdf" % (filePrefix, fileNamePrefix)
				os.system(cmd)
			else:
				# os.system("rm %s.mol2" % fileNamePrefix)
				pass
		# babel_convert("%s.mol2" % fileNamePrefix, "mol2", "%s" % visuSdfFileName, "sdf", verbose = verbose)
		if not split and cleanup:
			# os.system("rm  %s[.,_]*" % (filePrefix))
			pass
		babel_convert("%s.mol2" % filePrefix, "mol2", "%s.sdf" % filePrefix, "sdf", verbose = verbose)
		cmd = "cat %s.sdf >> %s" % (filePrefix, visuSdfFileName)
		os.system(cmd)

		self.eOK = (float(bestEnergy) < eIni)
		# sys.stderr.write("Will return %d \n" % (int(nConfGen)))
		return int(nConfGen)


	def isAxialEquatorial(self, atom):
		"""
		Is the atom positoin axial/eqautorial ?
		"""
		nbCycles = 0
		indexCycle = -1
		for index, ring in enumerate(atom.rings):
			if atom in ring.atoms:
				nbCycles += 1
				indexCycle = index
		if nbCycles == 1 and len(atom.rings[indexCycle].atoms) > 4:
			#print atom, neighboursAtomsOutsideCycle(atom)
			neighbours = neighboursAtomsOutsideCycle(atom)
			if len(neighbours) == 2:
				chainsNeighbours = []
				for oatom in neighbours:
					chainsNeighbours.append(toCanonicalChain(oatom,oatom.findbond(atom)))
				if chainsNeighbours[0] != chainsNeighbours[1]:
					for oatom in neighbours:
						if oatom.symbol != 'H':
							return True, oatom.index
		return False, None

	def ambiguousStereoAtoms(self, verbose = 0):
		stereoAtoms = []
		# 1st pass: crude
		for atom in self.molecule.atoms:
			if (not atom.chiral_class) and isChiralAtom(atom):
				stereoAtoms.append(atom)
				if verbose:
					sys.stderr.write("Atom %d %s is chiral\n" % (atom.index, atom.symbol))
		# 2nd pass: detect positions induced by other steroSites (in rings)
		fRings = self.detectFusedRings()
		# Atom in ring, surrounded by chiral sites -> maybe chiral
		# TODO
		
		return stereoAtoms


	def toNsmiles_new(self,axialEquatorialPos=False, verbose = 0):
		"""
		toNsmiles: liste tous les smiles non ambigus a partir du smiles de la molecule
		axialEquatorialPos: True, sort aussi la combinatoire des positions axial/equatorial

		ex: self.toNsmiles(axialEquatorialPos=True)
		return a List of 2 lists:
		First is the list of the disambiguate smiles
		Second is the list of strings. Each is of size the number of atoms of the smiles.
		Symbols are: x, A, E. A for Axial, E for Equatorial.
		The order of the string matches that of the atoms of the smile.
		each AE string can be used against each smiles of the first list

		If the smiles is not ambiguous, then the smiles is returned.
		"""
		# check for chiral atoms (@, axial/equatorial)
		# Smiles:
		# [ Square brackets delimit individual atoms
		# - : single bond (omitted)
		# = : double bond
		# # : triple bond
		# . : non bond. WE SHOULD NOT ENCOUNTER IT HERE
		# ( : branches
		# 0..9 : ring closure
		# lower case: aromaticity
		# / /\: Z (//) E (/\) 
		# @ : chirality (@ anticlockwise, @@ clockwise, Looking FROM the 1st neighbor listed in the SMILES TO the chiral atom, 
		#     the other three neighbors appear anticlockwise or clockwise in the order listed. )
		# > : reactions. WE SHOULD NOT ENCOUNTER IT HERE

		# Only works from smiles
		if not self.smiles:
			return [ [], [] ]

		if verbose:
			sys.stderr.write("toNsmiles: Will identify @ and axial/equatorial sites\n")
		
		# @ / @@ stuff
		stereoAtoms = []
		# 1st pass: crude
		for atom in self.molecule.atoms:
			if (not atom.chiral_class) and isChiralAtom(atom):
				stereoAtoms.append(atom)
		# 2nd pass: detect positions induced by other steroSites (in rings)
		fRings = self.detectFusedRings()
		for aRing in fRings:
			pass
		# Axial / Equatorial stuff
		AEatoms = []
		for atom in self.molecule.atoms:
			if atom in stereoAtoms:
				continue
			isAxEq, index =  self.isAxialEquatorial(atom)
			if isAxEq:
				AEatoms.append(index)


		stereoIsos = []
		i = 0
		indexAtom = 0
		indexBond = 0
		openBracket = False

		# First atom
		atom = self.molecule.atoms[indexAtom]
		symbol = atom.symbol
		if (len(symbol) == 1) and (atom.aromatic):
			symbol = symbol.lower()
		# find occurrence of atom symbol in smiles
		# Should it be the exact atom ??
		# CHECK MAPPING BETWEEN SMILES AND molecule ATOM ORDER
		pos = self.smiles.find(symbol)
		while i < len(self.smiles):
			#print i
			# Square brackets delimit individual atoms. Cascading brackets should not exist
			if self.smiles[i] == '[':
				openBracket = True
			elif i == pos:
				# This is a position where chirality is unspecified
				if (not atom.chiral_class) and isChiralAtom(atom):
					if openBracket:
						# Should we really get there ?
						# This avoids cascading brackets
						# But it could be erroneous
						pass
					else:
						explH = ""
						#if len(atom.oatoms) == 2:
						#if (atom.explicit_hcount == 1) or (atom.imp_hcount == 1):
						if getNbOfHInOatoms(atom) == 1:
							explH = 'H'
						# Expand smiles to have explicit H
						# Here, we should check if [H] arn't already specified ...
						self.smiles = self.smiles[:i] + '[' + self.smiles[i:i+len(atom.symbol)] + explH + \
														']' + self.smiles[i+len(atom.symbol):]								
						i += 1
					stereoIsos.append(("RS",i+1))
				# if we want to, we check for axialEquatorial: atom belongs to a non aromatic ring
				elif (not atom.chiral_class) and (axialEquatorialPos) and (atom.rings) and (not atom.aromatic):
					nbCycles = 0
					indexCycle = -1
					for index, ring in enumerate(atom.rings):
						if atom in ring.atoms:
							nbCycles += 1
							indexCycle = index
					if nbCycles == 1 and len(atom.rings[indexCycle].atoms) > 4:
						#print atom, neighboursAtomsOutsideCycle(atom)
						neighbours = neighboursAtomsOutsideCycle(atom)
						if len(neighbours) == 2:
							chainsNeighbours = []
							for oatom in neighbours:
								chainsNeighbours.append(toCanonicalChain(oatom,oatom.findbond(atom)))
							if chainsNeighbours[0] != chainsNeighbours[1]:
								for oatom in neighbours:
									if oatom.symbol != 'H':
										axialEqu.append(oatom.index)
				indexAtom += 1
				if indexAtom < self.numberOfNonHAtoms:
				    #if indexAtom < len(self.molecule.atoms):
					atom = self.molecule.atoms[indexAtom]
					symbol = atom.symbol
					if (len(symbol) == 1) and (atom.aromatic):
						symbol = symbol.lower()
					pos = self.smiles.find(symbol,i+1)
				else:
					pos = -1
					break
				openBracket = False
			# Loop over next smiles position
			i += 1
			#print "pos", pos, symbol, i, atom
		if verbose:
			sys.stderr.write("toNsmiles: found %d @/@@ sites and %d axial/equatorial sites\n" % (len(stereoIsos), len(axialEqu)))

		# check for Z/E bonds
		if verbose:
			sys.stderr.write("toNsmiles: Will identify ZE sites\n")
			sys.stderr.write("toNsmiles: smiles is %s\n" % self.smiles)
		# We start from indexBond = 0
		while indexBond < (len(self.molecule.bonds) - len(self.molecule.atoms) + self.numberOfNonHAtoms):
			atoms = self.molecule.bonds[indexBond].atoms
			isChiral = False
			for bond in atoms[0].bonds:
				if (bond.symbol =='\\') or (bond.symbol =='/'):
					isChiral = True
			# We look for sites of unspecified chirality
			# What is bonds[].dbo ??
			# atoms are the atoms defining the bond of interest
			if (not isChiral) and ((not self.molecule.bonds[indexBond].dbo) and (isChiralBond(self.molecule.bonds[indexBond]))):
				if verbose:
					sys.stderr.write("considering chiral bond: %s %d - %s %d\n" % (atoms[0].symbol, atoms[0].index, atoms[1].symbol, atoms[1].index))
				# atoms: the atoms of the bond
				# we sort them according to their index
				atoms.sort(key = lambda x: x.index)
				
				# is the first atom of the bond the first atom ? no : ...
				if atoms[0].index != 0:

					# Not 1st position
					indexAtom = 0
					# find the position of indexatom in smiles
					atom = self.molecule.atoms[indexAtom]
					symbol = self.molecule.atoms[indexAtom].symbol
					if (len(symbol) == 1) and (atom.aromatic):
						symbol = symbol.lower()
					pos = self.smiles.find(symbol)
					indexBracket = 0
					openBracket = False
					if verbose:
						sys.stderr.write("atom %d at smiles position %d\n" % (indexAtom, pos) )

					# i is position in smiles string
					# pos is position in smiles of the indexAtom atom 
					# indexAtom is ??? (mapping smiles atoms with molecule atoms ?)
					i = 0
					while i < len(self.smiles):
						if self.smiles[i] == '[':
							openBracket = True
							indexBracket = i
						# We consider the next atom of correct symbol
						elif i == pos:
							# We are on the correct index of first bond atom
							if indexAtom == atoms[0].index:
								init = i
								if openBracket:
									init = indexBracket
								if verbose:
									sys.stderr.write("Found init atom at smiles %d\n" % i )
								# Now we shall search for second atom.
								# We start from next indexAtom
								indexAtom += 1
								symbol = self.molecule.atoms[indexAtom].symbol
								# if (len(symbol) == 1) and (atom.aromatic):
								if (len(symbol) == 1) and (self.molecule.atoms[indexAtom].aromatic):
									symbol = symbol.lower()
								pos = self.smiles.find(symbol,i+1)
								if verbose:
									sys.stderr.write("second atom search starts as atom %d %s at smiles position %d (search from %d)\n" % (indexAtom, symbol, pos, i+1) )

								openBracket = False
								while True:
									# Uncomment next line if bug here
									# print "COUCOU TRUE", self.smiles, i
									i += 1
									if i >= len(self.smiles):
										sys.stderr.write("toNsmiles: overflow on smiles walk\n")
										raise Exception("toNsmiles: overflow on smiles walk")
									if self.smiles[i] == '[':
										openBracket = True
										indexBracket = i
									elif i == pos:
										if indexAtom == atoms[1].index:
											if verbose:
												sys.stderr.write("Found end atom at smiles %d\n" % i )
											if openBracket:
												i += 1
												while self.smiles[i] != ']':
													i += 1
											i += 1
											if self.smiles[i] == '(':
												i += 1
											break
										indexAtom += 1
										atom = self.molecule.atoms[indexAtom]
										symbol = atom.symbol
										if (len(symbol) == 1) and (atom.aromatic):
											symbol = symbol.lower()
										#print "self.smiles.find(symbol,i+1)", symbol,i+1, "\n"
										pos = self.smiles.find(symbol,i+1)
										openBracket = False
								tmp = i
								fin = i
								#print "fin", i
								while self.smiles[i] != '=':
									i -= 1
								if verbose:
									sys.stderr.write("Quoting end at smiles %d\n" % i )
								
								while not self.smiles[i].isalpha():
									if self.smiles[i] == ')':
										break
									if self.smiles[i].isdigit():
										digit = self.smiles[i]
										i -= 1
										while self.smiles[i] != digit:
											i -= 1
										#init = i 
										fin = i
										break
									i -= 1	
								if fin < init:
									init = tmp
								stereoIsos.append(("ZE1",init))
								stereoIsos.append(("ZE2",fin))
								break
							if verbose:
								sys.stderr.write( "No on correct index loop\n")
							# We are not on the correct index
							openBracket = False
							indexAtom += 1
							if indexAtom < self.numberOfNonHAtoms:
							#if indexAtom < len(self.molecule.atoms):
								atom = self.molecule.atoms[indexAtom]
								symbol = atom.symbol
								if (len(symbol) == 1) and (atom.aromatic):
									symbol = symbol.lower()
								pos = self.smiles.find(symbol,i+1)
								if verbose:
									sys.stderr.write("atom %d at smiles position %d\n" % (indexAtom, pos) )
								#print "pos", pos, indexAtom
							else:
								# Last position
								pos = -1
								break
						i += 1	


				else:

					# print "en 1ere position"
					init = lg = len(self.molecule.atoms[0].symbol)
					envers = False
					if self.smiles[lg] == '(':
						init = lg + 1
						if self.smiles[lg+1] == '=':
							init += len(self.molecule.atoms[1].symbol) + 1
							if self.smiles[init] == '(':
								init += 1
							envers = True
							#print atoms[0].oatoms
							nextIndex = atoms[0].oatoms[1].index
							#if atoms[0].oatoms[0].symbol == 'H':
							#	nextIndex = atoms[0].oatoms[2].index
							#print "envers"
					elif self.smiles[lg].isdigit():
						if self.smiles[lg+1].isdigit():
							init = lg + 1
						else:
							init = lg
					if verbose:
						sys.stderr.write("Found special init atom at smiles %d\n" % init )
					stereoIsos.append(("ZE1",init))
					indexAtom = 0
					i = 0
					atom = self.molecule.atoms[indexAtom]
					symbol = atom.symbol
					if (len(symbol) == 1) and (atom.aromatic):
						symbol = symbol.lower()
					pos = self.smiles.find(symbol)
					indexBracket = 0
					openBracket = False	
					while i < len(self.smiles):
						if self.smiles[i] == '[':
							openBracket = True
							indexBracket = i
						elif i == pos:
							if (envers) and (indexAtom == nextIndex):
								if openBracket:
									i -= 1
								stereoIsos.append(("ZE2",i))
								if verbose:
									sys.stderr.write("Found special end atom at smiles %d\n" % i )
								break
							elif (not envers) and (indexAtom == atoms[1].index):
								if openBracket:
									i += 1
									while self.smiles[i] != ']':
										i += 1
								i += 1
								if self.smiles[i] == '(':
									i += 1
								stereoIsos.append(("ZE2",i))
								if verbose:
									sys.stderr.write("Found special end atom at smiles %d\n" % i )
								break
							openBracket = False
							indexAtom += 1
							if indexAtom < self.numberOfNonHAtoms:
								atom = self.molecule.atoms[indexAtom]
								symbol = atom.symbol
								if (len(symbol) == 1) and (atom.aromatic):
									symbol = symbol.lower()
								pos = self.smiles.find(symbol,i+1)
							else:
								pos = -1
								break			
						i += 1
					# print "FIN en  1ere position"
			indexBond += 1

		if verbose:
			sys.stderr.write("toNsmiles: found %d @/@@ and ZE sites\n" % (len(stereoIsos)) )

		#print self.smiles	
		stereoIsos.sort(key = lambda x: x[1])
		if verbose:
			for s in stereoIsos:
				print s
				# sys.stderr.write("toNsmiles: %s\n" % (s) )
		print "toNsmiles stereoIsos"
		#print stereoIsos
		lg = len(stereoIsos)
		deleted = 0
		for i in range(lg-1):
			if stereoIsos[i-deleted][1] == stereoIsos[i+1-deleted][1]:
				#stereoIsos[i][0] = "ZE2"
				stereoIsos[i+1-deleted] = ("ZEX", stereoIsos[i+1-deleted][1])
				del stereoIsos[i-deleted]			
				deleted += 1
		#print stereoIsos
		if verbose:
			for s in stereoIsos:
				print s

		if verbose:
			sys.stderr.write("toNsmiles: will generate unambiguous smiles for %d enantiomer sites\n" % len(stereoIsos))
		
		nSmiles = []
		if stereoIsos:
			backS = "\\"
			stereos = {"RS":('@','@@'), "ZE1":('/','/'), "ZE2":('/',"\\"), "ZEX":('/',"\\")}
			nb = 0
			for i in range(len(stereoIsos)):
				if stereoIsos[i][0] in ["RS","ZE1","ZEX"]:
					nb += 1
			nSmiles = [""]*pow(2,nb)
			if verbose:
				sys.stderr.write("toNsmiles: will generate %d unambiguous smiles\n" % len(nSmiles))
			
			indexStereo = i = nbZE2 = 0
			while i < len(self.smiles):
				# whichStereo is isomeric state (over two)
				whichStereo = 0
				if stereoIsos[indexStereo][0] == "ZEX":
				#if stereoIsos[indexStereo][0] == "ZE2":
					nbZE2 += 1
				for j in xrange(len(nSmiles)):
					nSmiles[j] = nSmiles[j] + self.smiles[i:stereoIsos[indexStereo][1]]
					#if stereoIsos[indexStereo][0] == "ZE1":
						#nSmiles[j] = nSmiles[j] + "/" + self.smiles[stereoIsos[indexStereo][1]:stereoIsos[indexStereo][2]]
					#print j, deleted, len(nSmiles), pow(2,indexStereo+1)
					print len(nSmiles), (len(nSmiles)*pow(2,deleted))/pow(2,indexStereo+1-nbZE2), j
					if (indexStereo+1 == len(stereoIsos)) or \
					   ((j != 0) and ((len(nSmiles)*pow(2,deleted))/pow(2,indexStereo+1-nbZE2) != 0) and (j%((len(nSmiles)*pow(2,deleted))/pow(2,indexStereo+1-nbZE2)) == 0)):
						whichStereo = (whichStereo+1)%2
					nSmiles[j] = nSmiles[j] + stereos[stereoIsos[indexStereo][0]][whichStereo]
				i = stereoIsos[indexStereo][len(stereoIsos[indexStereo])-1]
				indexStereo += 1
				if indexStereo == len(stereoIsos):
					for j in xrange(len(nSmiles)):
						nSmiles[j] = nSmiles[j] + self.smiles[i:]
					break
			if verbose:
				for i in range(len(nSmiles)):
					sys.stderr.write("%s\n" % nSmiles[i])
			for i in range(len(nSmiles)):
				nSmiles[i].replace(' ','')
				#nSmiles[i].replace("\\\\","\\")
			# print nSmiles
			for i in range(len(nSmiles)):
				nSmiles[i].replace('//','/')
				nSmiles[i].replace('/\\','/')
				nSmiles[i].replace('\\/','\\')
				nSmiles[i].replace('\\\\','\\')
			for i in range(len(nSmiles)):
				nSmiles[i].replace('//','/')
				nSmiles[i].replace('/\\','/')
				nSmiles[i].replace('\\/','\\')
				nSmiles[i].replace('\\\\','\\')	
		else:
			nSmiles = [self.smiles]
		if verbose:
			sys.stderr.write("toNsmiles: Will generate nAxEq\n")

		nAxEq = []
		if axialEquatorialPos and axialEqu:
			# print "axialEqu",axialEqu
			nAxEq = [""]*pow(2,len(axialEqu))
			demi = len(nAxEq)
			i = 0
			axEq = ("A", "E")
			while i < self.numberOfNonHAtoms:
				whichDescrip = 0
				if i in axialEqu:
					demi /= 2
				for j in range(len(nAxEq)):
					if i in axialEqu:
						if (j != 0) and (j%(demi) == 0):
							whichDescrip = (whichDescrip+1)%2
						nAxEq[j] += axEq[whichDescrip]
					else:
						nAxEq[j] += "x"
				i += 1
		else:
			nAxEq = [ "x" * self.numberOfNonHAtoms ]
		if verbose:
			sys.stderr.write("toNsmiles: returning ...\n" )
		return [nSmiles, nAxEq]


	def toNsmiles(self, axialEquatorialPos=False, verbose = 0):
		"""
		toNsmiles: liste tous les smiles non ambigus a partir du smiles de la molecule
		axialEquatorialPos: True, sort aussi la combinatoire des positions axial/equatorial

		ex: self.toNsmiles(axialEquatorialPos=True)
		return a List of 2 lists:
		First is the list of the disambiguate smiles
		Second is the list of strings. Each is of size the number of atoms of the smiles.
		Symbols are: x, A, E. A for Axial, E for Equatorial.
		The order of the string matches that of the atoms of the smile.
		each AE string can be used against each smiles of the first list

		If the smiles is not ambiguous, then the smiles is returned.
		"""
		if not self.smiles:
			return [ [], [] ]

		axialEqu = []
		stereoIsos = []
		unstereoIsos = [] # Sites for which we remove iso info
		i = 0
		indexAtom = 0
		indexBond = 0
		openBracket = False

		# First atom
		atom = self.molecule.atoms[indexAtom]
		symbol = atom.symbol
		if (len(symbol) == 1) and (atom.aromatic):
			symbol = symbol.lower()
		# find occurrence of atom symbol in smiles
		# Should it be the exact atom ??
		# CHECK MAPPING BETWEEN SMILES AND molecule ATOM ORDER
		pos = self.smiles.find(symbol)
		#print self.smiles
		#print "pos", pos, symbol

		# check for chiral atoms (@, axial/equatorial)
		# Smiles:
		# [ Square brackets delimit individual atoms
		# - : single bond (omitted)
		# = : double bond
		# # : triple bond
		# . : non bond. WE SHOULD NOT ENCOUNTER IT HERE
		# ( : branches
		# 0..9 : ring closure
		# lower case: aromaticity
		# / /\: Z (//) E (/\) 
		# @ : chirality (@ anticlockwise, @@ clockwise, Looking FROM the 1st neighbor listed in the SMILES TO the chiral atom, 
		#     the other three neighbors appear anticlockwise or clockwise in the order listed. )
		# > : reactions. WE SHOULD NOT ENCOUNTER IT HERE

		if verbose:
			sys.stderr.write("toNsmiles: Will identify @ and axial/equatorial sites\n")
		while i < len(self.smiles):
			#print i
			# Square brackets delimit individual atoms. Cascading brackets should not exist
			if self.smiles[i] == '[':
				openBracket = True
			elif i == pos:
				if (atom.chiral_class) and isRingedRing(atom, verbose = verbose):
					atom.chiral_class = None
					unstereoIsos.append((atom, i+1))
					if verbose:
						sys.stderr.write("REMOVING CHIRAL INFORMATION FOR ATOM %s %d\n" % (atom.symbol, atom.index))
					if self.smiles[i+1] == "@":
						self.smiles = self.smiles[:i+1] + self.smiles[i+2:]	
					if self.smiles[i+1] == "@":  # For @@
						self.smiles = self.smiles[:i+1] + self.smiles[i+2:]	
					if verbose:
						sys.stderr.write("%s\n" % self.smiles)

					#print self.molecule.cansmiles()
				# This is a position where chirality is unspecified
				if (not atom.chiral_class) and isChiralAtom(atom, verbose = verbose):
					if openBracket:
						# Should we really get there ?
						pass
					else:
						explH = ""
						#if len(atom.oatoms) == 2:
						#if (atom.explicit_hcount == 1) or (atom.imp_hcount == 1):
						if getNbOfHInOatoms(atom) == 1:
							explH = 'H'
						# Expand smiles to have explicit H
						# Here, we should check if [H] arn't already specified ...
						self.smiles = self.smiles[:i] + '[' + self.smiles[i:i+len(atom.symbol)] + explH + \
														']' + self.smiles[i+len(atom.symbol):]								
						i += 1
					stereoIsos.append(("RS",i+1))
				# if we want to, we check for axialEquatorial: atom belongs to a non aromatic ring
				elif (not atom.chiral_class) and (axialEquatorialPos) and (atom.rings) and (not atom.aromatic):
					nbCycles = 0
					indexCycle = -1
					for index, ring in enumerate(atom.rings):
						if atom in ring.atoms:
							nbCycles += 1
							indexCycle = index
					if nbCycles == 1 and len(atom.rings[indexCycle].atoms) > 4:
						#print atom, neighboursAtomsOutsideCycle(atom)
						neighbours = neighboursAtomsOutsideCycle(atom)
						if len(neighbours) == 2:
							chainsNeighbours = []
							for oatom in neighbours:
								chainsNeighbours.append(toCanonicalChain(oatom,oatom.findbond(atom)))
							if chainsNeighbours[0] != chainsNeighbours[1]:
								for oatom in neighbours:
									if oatom.symbol != 'H':
										axialEqu.append(oatom.index)
				indexAtom += 1
				if indexAtom < self.numberOfNonHAtoms:
				    #if indexAtom < len(self.molecule.atoms):
					atom = self.molecule.atoms[indexAtom]
					symbol = atom.symbol
					if (len(symbol) == 1) and (atom.aromatic):
						symbol = symbol.lower()
					pos = self.smiles.find(symbol,i+1)
				else:
					pos = -1
					break
				openBracket = False
			# Loop over next smiles position
			i += 1
			#print "pos", pos, symbol, i, atom
		if verbose:
			sys.stderr.write("toNsmiles: found %d @/@@ sites and %d axial/equatorial sites\n" % (len(stereoIsos), len(axialEqu)))

		# TO DEBUG @ WARNING: DANGEROUS:
		# axialEqu = []

		# check for Z/E bonds
		if verbose:
			sys.stderr.write("toNsmiles: Will identify ZE sites\n")
			sys.stderr.write("toNsmiles: smiles is %s\n" % self.smiles)
		# We start from indexBond = 0
		while indexBond < (len(self.molecule.bonds) - len(self.molecule.atoms) + self.numberOfNonHAtoms):
			atoms = self.molecule.bonds[indexBond].atoms
			isChiral = False
			for bond in atoms[0].bonds:
				if (bond.symbol =='\\') or (bond.symbol =='/'):
					isChiral = True
			# We look for sites of unspecified chirality
			# What is bonds[].dbo ??
			# atoms are the atoms defining the bond of interest
			if (not isChiral) and ((not self.molecule.bonds[indexBond].dbo) and (isChiralBond(self.molecule.bonds[indexBond]))):
				if verbose:
					sys.stderr.write("considering chiral bond: %s %d - %s %d\n" % (atoms[0].symbol, atoms[0].index, atoms[1].symbol, atoms[1].index))
				# atoms: the atoms of the bond
				# we sort them according to their index
				atoms.sort(key = lambda x: x.index)
				
				# is the first atom of the bond the first atom ? no : ...
				if atoms[0].index != 0:

					# Not 1st position
					indexAtom = 0
					# find the position of indexatom in smiles
					atom = self.molecule.atoms[indexAtom]
					symbol = self.molecule.atoms[indexAtom].symbol
					if (len(symbol) == 1) and (atom.aromatic):
						symbol = symbol.lower()
					pos = self.smiles.find(symbol)
					indexBracket = 0
					openBracket = False
					if verbose:
						sys.stderr.write("atom %d at smiles position %d\n" % (indexAtom, pos) )

					# i is position in smiles string
					# pos is position in smiles of the indexAtom atom 
					# indexAtom is ??? (mapping smiles atoms with molecule atoms ?)
					i = 0
					while i < len(self.smiles):
						if self.smiles[i] == '[':
							openBracket = True
							indexBracket = i
						# We consider the next atom of correct symbol
						elif i == pos:
							# We are on the correct index of first bond atom
							if indexAtom == atoms[0].index:
								init = i
								if openBracket:
									init = indexBracket
								if verbose:
									sys.stderr.write("Found init atom at smiles %d\n" % i )
								# Now we shall search for second atom.
								# We start from next indexAtom
								indexAtom += 1
								symbol = self.molecule.atoms[indexAtom].symbol
								# if (len(symbol) == 1) and (atom.aromatic):
								if (len(symbol) == 1) and (self.molecule.atoms[indexAtom].aromatic):
									symbol = symbol.lower()
								pos = self.smiles.find(symbol,i+1)
								if verbose:
									sys.stderr.write("second atom search starts as atom %d %s at smiles position %d (search from %d)\n" % (indexAtom, symbol, pos, i+1) )

								openBracket = False
								while True:
									# Uncomment next line if bug here
									# print "COUCOU TRUE", self.smiles, i
									i += 1
									if i >= len(self.smiles):
										sys.stderr.write("toNsmiles: overflow on smiles walk\n")
										raise Exception("toNsmiles: overflow on smiles walk")
									if self.smiles[i] == '[':
										openBracket = True
										indexBracket = i
									elif i == pos:
										if indexAtom == atoms[1].index:
											if verbose:
												sys.stderr.write("Found end atom at smiles %d\n" % i )
											if openBracket:
												i += 1
												while self.smiles[i] != ']':
													i += 1
											i += 1
											if self.smiles[i] == '(':
												i += 1
											break
										indexAtom += 1
										atom = self.molecule.atoms[indexAtom]
										symbol = atom.symbol
										if (len(symbol) == 1) and (atom.aromatic):
											symbol = symbol.lower()
										#print "self.smiles.find(symbol,i+1)", symbol,i+1, "\n"
										pos = self.smiles.find(symbol,i+1)
										openBracket = False
								tmp = i
								fin = i
								#print "fin", i
								while self.smiles[i] != '=':
									i -= 1
								if verbose:
									sys.stderr.write("Quoting end at smiles %d\n" % i )
								
								while not self.smiles[i].isalpha():
									if self.smiles[i] == ')':
										break
									if self.smiles[i].isdigit():
										digit = self.smiles[i]
										i -= 1
										while self.smiles[i] != digit:
											i -= 1
										#init = i 
										fin = i
										break
									i -= 1	
								if fin < init:
									init = tmp
								stereoIsos.append(("ZE1",init))
								stereoIsos.append(("ZE2",fin))
								break
							if verbose:
								sys.stderr.write( "No on correct index loop\n")
							# We are not on the correct index
							openBracket = False
							indexAtom += 1
							if indexAtom < self.numberOfNonHAtoms:
							#if indexAtom < len(self.molecule.atoms):
								atom = self.molecule.atoms[indexAtom]
								symbol = atom.symbol
								if (len(symbol) == 1) and (atom.aromatic):
									symbol = symbol.lower()
								pos = self.smiles.find(symbol,i+1)
								if verbose:
									sys.stderr.write("atom %d at smiles position %d\n" % (indexAtom, pos) )
								#print "pos", pos, indexAtom
							else:
								# Last position
								pos = -1
								break
						i += 1	


				else:

					# print "en 1ere position"
					init = lg = len(self.molecule.atoms[0].symbol)
					envers = False
					if self.smiles[lg] == '(':
						init = lg + 1
						if self.smiles[lg+1] == '=':
							init += len(self.molecule.atoms[1].symbol) + 1
							if self.smiles[init] == '(':
								init += 1
							envers = True
							#print atoms[0].oatoms
							nextIndex = atoms[0].oatoms[1].index
							#if atoms[0].oatoms[0].symbol == 'H':
							#	nextIndex = atoms[0].oatoms[2].index
							#print "envers"
					elif self.smiles[lg].isdigit():
						if self.smiles[lg+1].isdigit():
							init = lg + 1
						else:
							init = lg
					if verbose:
						sys.stderr.write("Found special init atom at smiles %d\n" % init )
					stereoIsos.append(("ZE1",init))
					indexAtom = 0
					i = 0
					atom = self.molecule.atoms[indexAtom]
					symbol = atom.symbol
					if (len(symbol) == 1) and (atom.aromatic):
						symbol = symbol.lower()
					pos = self.smiles.find(symbol)
					indexBracket = 0
					openBracket = False	
					while i < len(self.smiles):
						if self.smiles[i] == '[':
							openBracket = True
							indexBracket = i
						elif i == pos:
							if (envers) and (indexAtom == nextIndex):
								if openBracket:
									i -= 1
								stereoIsos.append(("ZE2",i))
								if verbose:
									sys.stderr.write("Found special end atom at smiles %d\n" % i )
								break
							elif (not envers) and (indexAtom == atoms[1].index):
								if openBracket:
									i += 1
									while self.smiles[i] != ']':
										i += 1
								i += 1
								if self.smiles[i] == '(':
									i += 1
								stereoIsos.append(("ZE2",i))
								if verbose:
									sys.stderr.write("Found special end atom at smiles %d\n" % i )
								break
							openBracket = False
							indexAtom += 1
							if indexAtom < self.numberOfNonHAtoms:
								atom = self.molecule.atoms[indexAtom]
								symbol = atom.symbol
								if (len(symbol) == 1) and (atom.aromatic):
									symbol = symbol.lower()
								pos = self.smiles.find(symbol,i+1)
							else:
								pos = -1
								break			
						i += 1
					# print "FIN en  1ere position"
			indexBond += 1

		if verbose:
			sys.stderr.write("toNsmiles: found %d @/@@ and ZE sites\n" % (len(stereoIsos)) )

		#print self.smiles	
		stereoIsos.sort(key = lambda x: x[1])
		# if verbose:
		# 	for s in stereoIsos:
		# 		print s
				# sys.stderr.write("toNsmiles: %s\n" % (s) )
		# print "toNsmiles stereoIsos"
		#print stereoIsos
		lg = len(stereoIsos)
		deleted = 0
		for i in range(lg-1):
			if stereoIsos[i-deleted][1] == stereoIsos[i+1-deleted][1]:
				#stereoIsos[i][0] = "ZE2"
				stereoIsos[i+1-deleted] = ("ZEX", stereoIsos[i+1-deleted][1])
				del stereoIsos[i-deleted]			
				deleted += 1
		#print stereoIsos
		# if verbose:
		# 	for s in stereoIsos:
		# 		print s

		if verbose:
			sys.stderr.write("toNsmiles: will generate unambiguous smiles for %d enantiomer sites\n" % len(stereoIsos))
		
		nSmiles = []
		if stereoIsos:
			backS = "\\"
			stereos = {"RS":('@','@@'), "ZE1":('/','/'), "ZE2":('/',"\\"), "ZEX":('/',"\\")}
			nb = 0
			for i in range(len(stereoIsos)):
				if stereoIsos[i][0] in ["RS","ZE1","ZEX"]:
					nb += 1
			nSmiles = [""]*pow(2,nb)
			if verbose:
				sys.stderr.write("toNsmiles: will generate %d unambiguous smiles\n" % len(nSmiles))
			
			indexStereo = i = nbZE2 = 0
			while i < len(self.smiles):
				# whichStereo is isomeric state (over two)
				whichStereo = 0
				if stereoIsos[indexStereo][0] == "ZEX":
				#if stereoIsos[indexStereo][0] == "ZE2":
					nbZE2 += 1
				for j in xrange(len(nSmiles)):
					nSmiles[j] = nSmiles[j] + self.smiles[i:stereoIsos[indexStereo][1]]
					#if stereoIsos[indexStereo][0] == "ZE1":
						#nSmiles[j] = nSmiles[j] + "/" + self.smiles[stereoIsos[indexStereo][1]:stereoIsos[indexStereo][2]]
					#print j, deleted, len(nSmiles), pow(2,indexStereo+1)
					# print len(nSmiles), (len(nSmiles)*pow(2,deleted))/pow(2,indexStereo+1-nbZE2), j
					if (indexStereo+1 == len(stereoIsos)) or \
					   ((j != 0) and ((len(nSmiles)*pow(2,deleted))/pow(2,indexStereo+1-nbZE2) != 0) and (j%((len(nSmiles)*pow(2,deleted))/pow(2,indexStereo+1-nbZE2)) == 0)):
						whichStereo = (whichStereo+1)%2
					nSmiles[j] = nSmiles[j] + stereos[stereoIsos[indexStereo][0]][whichStereo]
				i = stereoIsos[indexStereo][len(stereoIsos[indexStereo])-1]
				indexStereo += 1
				if indexStereo == len(stereoIsos):
					for j in xrange(len(nSmiles)):
						nSmiles[j] = nSmiles[j] + self.smiles[i:]
					break
			if verbose:
				for i in range(len(nSmiles)):
					sys.stderr.write("%s\n" % nSmiles[i])
			for i in range(len(nSmiles)):
				nSmiles[i].replace(' ','')
				#nSmiles[i].replace("\\\\","\\")
			# print nSmiles
			for i in range(len(nSmiles)):
				nSmiles[i].replace('//','/')
				nSmiles[i].replace('/\\','/')
				nSmiles[i].replace('\\/','\\')
				nSmiles[i].replace('\\\\','\\')
			for i in range(len(nSmiles)):
				nSmiles[i].replace('//','/')
				nSmiles[i].replace('/\\','/')
				nSmiles[i].replace('\\/','\\')
				nSmiles[i].replace('\\\\','\\')	
		else:
			nSmiles = [self.smiles]
		if verbose:
			sys.stderr.write("toNsmiles: Will generate nAxEq\n")

		nAxEq = []
		if axialEquatorialPos and axialEqu:
			# print "axialEqu",axialEqu
			nAxEq = [""]*pow(2,len(axialEqu))
			demi = len(nAxEq)
			i = 0
			axEq = ("A", "E")
			while i < self.numberOfNonHAtoms:
				whichDescrip = 0
				if i in axialEqu:
					demi /= 2
				for j in range(len(nAxEq)):
					if i in axialEqu:
						if (j != 0) and (j%(demi) == 0):
							whichDescrip = (whichDescrip+1)%2
						nAxEq[j] += axEq[whichDescrip]
					else:
						nAxEq[j] += "x"
				i += 1
		else:
			nAxEq = [ "x" * self.numberOfNonHAtoms ]
		if verbose:
			sys.stderr.write("toNsmiles: returning ...\n" )
		return [nSmiles, nAxEq, len(unstereoIsos)]

		
	
		
		
		
		
		
		
		
		
		
		

