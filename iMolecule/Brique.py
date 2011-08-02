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

import re
from Exceptions import CycleUnknown
from frowns import Smiles, Cycle
from frowns import MDL
import os
import sys
from frowns.Canonicalization import Disambiguate

from Config import *

class Brique :
    """
    une brique est un fragment d'une molecule
    """
    def __init__(self,m,c) :
        self.molecule = m
        self.connectInfo = c
        self.handle = id(self)  # Unique Id, P. Tuffery 2009

#     def __repr__(self):
#         return str(self.handle)

    def cansmiles(self) :
        """
        affiche en canonical smiles une brique
        """
        # sys.stderr.write("Brique.cansmiles: will call molecule.cansmiles() ...\n")
        a = self.molecule.cansmiles()

        m = Smiles.smilin("".join(re.split("[\[\]+-]","".join(re.split("H[0-9]*",a)))))
        return m.cansmiles()

    def write(self,file="test2",mode="w") :
        """
        ecrire une brique en sdf
        """
        ### HEADER
        result = [self.cansmiles()]
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
            result.append("%3i%3i%3i%3i%3i%3i"%(bond.atoms[0].index+1, bond.atoms[1].index+1, bond.bondorder,0, 0, 0))
        ### FIN DU FICHIER
        result.append("M  END")
        result.append("$$$$")
        
        f = open(file,mode)
        f.write("\n".join(result))
        f.write("\n")
        f.close()
        
    def load(self, workingPath, verbose = False) :
        """
        charger le fichier sdf de la brique
        """
        if verbose:
            sys.stderr.write("Brique.load: %s\n" % LIBRARY_PATH)
        if not LIBRARY_PATH :
            sys.stderr.wirte( "variable path de la classe Brique non initialisee\n")
            return 0
        try:
            listFile = os.listdir(LIBRARY_PATH) # liste des fichier contenu ds 'PAHT' -> listFile
        except:
            sys.stderr.write("Brique.load: Failed to access ring library at %s\n" % LIBRARY_PATH)
        # sys.exit(0)
        name = self.cansmiles()

        if verbose:
            sys.stderr.write("Brique.load: Attempting %s ...\n" % str(name))
            # print "Brique.load:",name
        if workingPath[-1] != "/":
            workingPath += "/"
    
        fileSmi = open(workingPath + "brique.smi", "w")
        fileSmi.write(name + "\n")
        fileSmi.close()
        if verbose:
            sys.stderr.write("Brique.load: babel sdf conversion ...\n")
        os.system(BABEL_PATH + " -ismi " + workingPath + "brique.smi -osdf " + workingPath + "brique.sdf 2> /dev/null")
        os.system("echo '$$$$' >> " + workingPath + "brique.sdf")
        # sys.stderr.write("Brique.load: Will read %s/brique.sdf\n" % (workingPath))
        fileSdf = open(workingPath + "brique.sdf")
        reader = MDL.sdin(fileSdf)
        name2 = reader.next()[0].cansmiles()
        # sys.stderr.write("Brique.load: read sdf as %s\n" %  name2)
        fileSdf.close()
            # print "brique", name, name2
        # if name2 in listFile: # si le fichier est ds la liste
        #     return LIBRARY_PATH+name2 # on renvoie le chemin complet du fichier
        # elif name in listFile:
        #     return LIBRARY_PATH+name
        if name in listFile:
            return LIBRARY_PATH+name
        else:
            # 11 fevrier 2010: We attempt to generate using AMMOS, not to use babel conversion for name.
            name2 = name
            if verbose:
                sys.stderr.write("Brique.load: Building %s as %s ...\n" % ( str(name), str(name2)))
            try:
                # Generer du mol2.
                uName = name2.replace("(","\(").replace(")","\)").replace("[","\[").replace("]","\]")
                # f  =open("%s/%s.smi" % (LIBRARY_PATH, name2), "w")
                sys.stderr.write("Creating %s%s.smi\n" % (LIBRARY_PATH, uName))
                f = open("%s%s.smi" % (LIBRARY_PATH, name2), "w") # in python do not escape names !!!
                f.write("%s\n" %  name2)
                f.close()
                cmd = "%s -h -ismi %s/%s.smi -omol2 %s/%s.mol2" % (BABEL_PATH, LIBRARY_PATH, uName, LIBRARY_PATH, uName)
                # cmd = "%s -h -ismi %s/%s.smi -omol2 %s/%s.mol2" % (BABEL_PATH, LIBRARY_PATH, "itest", LIBRARY_PATH, "itest")
                sys.stderr.write("%s\n" % cmd)
                os.system(cmd)

                # Generer le fichier input.param
                # path_of_DG-AMMOS= /Users/DG-AMMOS
                # bank= input_dataset.mol2
                f=  open("%s/%s.dgammos" % (LIBRARY_PATH, name2), "w")
                # f=  open("%s.dgammos" % ("itest"), "w")
                f.write("path_of_DG-AMMOS= %s\n" % DGAMMOSHOME)
                f.write("bank= %s%s.mol2\n" % (LIBRARY_PATH, name2))
                # f.write("bank= %s%s.mol2\n" % (LIBRARY_PATH, "itest"))
                f.close()
                sys.stderr.write("Generated %s/%s.dgammos\n" % (LIBRARY_PATH, uName))

                sys.stderr.write("Generating %s/%s.mol2\n" % (LIBRARY_PATH, uName))
                cmd = "%s %s/%s.dgammos" % (AMMOSBUILD, LIBRARY_PATH, uName)
                # cmd = "%s %s.dgammos" % (AMMOSBUILD, uName)
                sys.stderr.write("DG_AMMOS cmd:  %s\n" % (cmd))
                os.system(cmd)
                # sys.exit(0)

                # Convertir en sdf, oter les H
                # cmd = "%s -d -imol2 %s/%s_Built.mol2 -osdf %s/%s.sdf 2> %s/%s.log" % (BABEL_PATH, LIBRARY_PATH, name2, LIBRARY_PATH, name2, LIBRARY_PATH, name2)
                cmd = "%s -d -imol2 %s/%s_Built.mol2 -osdf %s/%s 2> %s/%s.log" % (BABEL_PATH, LIBRARY_PATH, uName, LIBRARY_PATH, uName, LIBRARY_PATH, uName)
                os.system(cmd)
                os.system("echo '$$$$' >> %s/%s" % (LIBRARY_PATH, uName))
                # sys.exit(0)
                # self.load(workingPath)
                return LIBRARY_PATH+name2
                # SI AMMOS N'A PAS PUS GENERER LE CYCLE: le mol2 sera dans %s/%s_BadMolecules.mol2
                # os.system(cmd)
            except:
                return 0
            return 0
            
	
    def map(self, workingPath, verbose = False):
        if verbose:
            sys.stderr.write("Brique.map: Will load %s ... wp %s\n" % (self.cansmiles(), workingPath))
        file = self.load(workingPath, verbose = verbose)
        if verbose:
            sys.stderr.write("Brique.map: %s Loaded ...\n" % workingPath)
        if not file :
            if verbose:
                sys.stderr.write("Brique.map: ring not found ...\n")
            # Here we should use AMMOS to generate ring conformation
            raise CycleUnknown, "cycle " + self.cansmiles() + " is not present in the ring library"
        # sys.stderr.write("Brique.map: will read %s ...\n" % file)
        f = open(file)
        reader = MDL.sdin(f)
        # m = reader.next()[0]
        m, err, txt = reader.next()
        
        # sys.stderr.write("Brique.map: MDL reader got %s\n" % err)
        f.close()
        # sys.stderr.write("Brique.map: read as SDF ...\n")
        Disambiguate.FreedDisambiguate(m)
        atomsSDF = m.atoms[:]
        atomsSDF.sort(key = lambda x: x.symorder)
        # sys.stderr.write("Brique.map: ordered by symorder ...\n")
        try:
            Disambiguate.FreedDisambiguate(self.molecule)
            atomsBrick = self.molecule.atoms[:]
            atomsBrick.sort(key = lambda x: x.symorder)
            # sys.stderr.write("Brique.map: Loaded %d atoms ...\n" % len(atomsSDF))
            for i in range(len(atomsBrick)):
                for oatB in atomsBrick[i].oatoms:
                    if oatB.symorder not in [oatS.symorder for oatS in atomsSDF[i].oatoms]:
                        for j in range(i):
                            atomsBrick[j].x = None
                            atomsBrick[j].y = None
                            atomsBrick[j].z = None
                            self.map2(workingPath)
                            return			
                atomsBrick[i].x = atomsSDF[i].x
                atomsBrick[i].y = atomsSDF[i].y
                atomsBrick[i].z = atomsSDF[i].z
        except:
            # sys.stderr.write("Brique.map: Something went wrong. Attempting map2\n")
            self.map2(workingPath)
            return

    def map2(self, workingPath) :
        """
        #correspondance entre atomes de la brique et du fichier sdf
        """
        file = self.load(workingPath)
        if not file :
            raise CycleUnknown, self.cansmiles()
        f = open(file)
        reader = MDL.sdin(f)
        m = reader.next()[0]
		# mx = reader.next() # ou while reader.next()
		# m = mx[0]
        f.close()

		# output : 
        
        i = j = tentativeRate = 0
        cyclesChecked = []
        cyclesMapped = []
        #print "nb cycles", len(self.molecule.cycles)
        while tentativeRate < len(self.molecule.cycles): # boucle tant que le nb de tentative rate de mapping n'est pas egale au nb de cycle ds la molecule
            #print j,i,len(self.molecule.cycles[j]),len(m.cycles[i])
            if len(self.molecule.cycles[j]) == len(m.cycles[i]) : # si les 2 cycles ont le meme nb d'atomes, on peut essayer de les mapper
                if mapCycles(self.molecule.cycles[j],m.cycles[i]) :
                    #print "mapCycles(self.molecule.cycles[j],m.cycles[i])"
                    cyclesChecked.append(i)
                    cyclesMapped.append(j)
                    j += 1
                    tentativeRate = 0
                    break
                else :
                    tentativeRate += 1
            else :
                tentativeRate += 1
            i += 1


        while len(cyclesMapped) < len(self.molecule.cycles) and tentativeRate < len(self.molecule.cycles)-len(cyclesChecked) :
            if i in cyclesChecked :
                i += 1
                continue
            if i == len(m.cycles) :
                i = 0
                continue
            if j in cyclesMapped :
                j += 1
                continue
            if j == len(self.molecule.cycles) :
                j = 0
                continue
            coord = 0
            for atom in self.molecule.cycles[j].atoms :
                if atom.x != None :
                    coord = 1
                    break
            if not coord :
                j += 1
                continue

            if len(self.molecule.cycles[j]) == len(m.cycles[i]) :
                if mapCycles2(self.molecule.cycles[j],m.cycles[i]) :
                    cyclesChecked.append(i)
                    cyclesMapped.append(j)
                    j += 1
                    tentativeRate = 0
                else :
                    tentativeRate += 1
            else :
                tentativeRate += 1
            i += 1
        
        return tentativeRate == 0

def mapCycles2(c1,c2):
    i = j = tentativeRate = 0
    atomsChecked = []
    while j < len(c1.atoms) and tentativeRate < len(c1.atoms)-len(atomsChecked) :
        if i in atomsChecked :
            i += 1
            continue
        if i == len(c2.atoms) :
            i = 0
            continue

        if c1.atoms[j].x != None :
            ok = 0
            for k in range(len(c2.atoms)) :
                if [c2.atoms[k].x,c2.atoms[k].y,c2.atoms[k].z] == [c1.atoms[j].x,c1.atoms[j].y,c1.atoms[j].z] :
                    atomsChecked.append(k)
                    ok = 1
                    break
            if not ok :
                tentativeRate = len(c1.atoms)
            j += 1
            continue

        if mapAtoms3(c1.atoms[j],c2.atoms[i],c1,c2,j,i) or mapAtoms4(c1.atoms[j],c2.atoms[i],c1,c2,j,i) :
            atomsChecked.append(i)
            c1.atoms[j].x = c2.atoms[i].x
            c1.atoms[j].y = c2.atoms[i].y
            c1.atoms[j].z = c2.atoms[i].z
            tentativeRate = 0
            j += 1
        else :
            tentativeRate += 1
        i += 1
    return j==len(c1.atoms)
    
def mapCycles(c1,c2):
    """
    verif si 2 cycles sont equivalents
    """
    i = j = tentativeRate = 0
    atomsChecked = []
    
    while j < len(c1.atoms) and tentativeRate < len(c1.atoms)-len(atomsChecked) :
        if i in atomsChecked :
            i += 1
            continue
        if i == len(c2.atoms) :
            i = 0
            continue
        
        if c1.atoms[j].x :
            ok = 0
            for k in range(len(c2.atoms)) :
                if [c2.atoms[k].x,c2.atoms[k].y,c2.atoms[k].z] == [c1.atoms[j].x,c1.atoms[j].y,c1.atoms[j].z] :
                    atomsChecked.append(k)
                    ok = 1
                    break
            if not ok :
                tentativeRate = len(c1.atoms)
            j += 1
            continue
        
        if mapAtoms(c1.atoms[j],c2.atoms[i],c1,c2,j,i) or mapAtoms2(c1.atoms[j],c2.atoms[i],c1,c2,j,i) :
	    atomsChecked.append(i)
            c1.atoms[j].x = c2.atoms[i].x
            c1.atoms[j].y = c2.atoms[i].y
            c1.atoms[j].z = c2.atoms[i].z
            tentativeRate = 0
            j += 1
        else :
            tentativeRate += 1
        i += 1
    return j==len(c1.atoms)

def mapAtoms3(a1,a2,c1,c2,index1,index2) :
    for k in range(1,len(c1.atoms)) :
        i = index1 + k
        if i >= len(c1.atoms) :
            i = i - len(c1.atoms)
                
        j = index2 + k
        if j >= len(c2.atoms) :
            j = j - len(c2.atoms)
        if c1.atoms[i].x != None :
            if [c1.atoms[i].x, c1.atoms[i].y, c1.atoms[i].z] == [c2.atoms[j].x, c2.atoms[j].y, c2.atoms[j].z] :
                continue
            else :
                return 0
    return 1


def mapAtoms(a1,a2,c1,c2,index1,index2) :
    """
    verifie si 2 atomes sont equivalents
    """
    if a1.number == a2.number :
        for k in range(1,len(c1.atoms)) :
            i = index1 + k
            if i >= len(c1.atoms) :
                i = i - len(c1.atoms)
                
            j = index2 + k
            if j >= len(c2.atoms) :
                j = j - len(c2.atoms)

            if c1.atoms[i].number == c2.atoms[j].number :
                
                if c1.atoms[i].x != None and c1.atoms[i].y != None and c1.atoms[i].z != None :
                    if [c1.atoms[i].x, c1.atoms[i].y, c1.atoms[i].z] == [c2.atoms[j].x, c2.atoms[j].y, c2.atoms[j].z] :
                        continue
                    else :
                        return 0
                #print "len rings", len(c1.atoms[i].rings), len(c2.atoms[j].rings)
                if len(c1.atoms[i].rings) == len(c2.atoms[j].rings) :
                    if len(c1.atoms[i].rings) > 1 :
                        if mapCyclesVoisins(c1.atoms[i],c2.atoms[j]) :
                            if mapVoisinsCyclesVoisins(c1.atoms[i],c2.atoms[j],c1,c2) :
                                #print "mapCyclesVoisins(c1.atoms[i],c2.atoms[j]) && mapVoisinsCyclesVoisins(c1.atoms[i],c2.atoms[j],c1,c2)"
				continue
                            else :
                                return 0
                        else :
                            return 0
                    else :
                        continue
                else :
                    return 0
                continue
            else :
                return 0
        return 1
    else :
        return 0

def mapAtoms4(a1,a2,c1,c2,index1,index2) :
    for k in range(1,len(c1.atoms)) :
        j = index2 - k
        if j < 0 :
            j = len(c2.atoms) + j

        i = index1 + k
        if i >= len(c1.atoms) :
            i = i - len(c1.atoms)

        if c1.atoms[i].x != None :
            if [c1.atoms[i].x, c1.atoms[i].y, c1.atoms[i].z] == [c2.atoms[j].x, c2.atoms[j].y, c2.atoms[j].z] :
                continue
            else :
                return 0
    return 1


def mapAtoms2(a1,a2,c1,c2,index1,index2) :
    if a1.number == a2.number :
        for k in range(1,len(c1.atoms)) :
            j = index2 - k
            if j < 0 :
                j = len(c2.atoms) + j

            i = index1 + k
            if i >= len(c1.atoms) :
                i = i - len(c1.atoms)

            if c1.atoms[i].number == c2.atoms[j].number :
                if c1.atoms[i].x != None and c1.atoms[i].y != None and c1.atoms[i].z != None :
                    if [c1.atoms[i].x, c1.atoms[i].y, c1.atoms[i].z] == [c2.atoms[j].x, c2.atoms[j].y, c2.atoms[j].z] :
                        continue
                    else :
                        return 0
                
                if len(c1.atoms[i].rings) == len(c2.atoms[j].rings) :
                    if len(c1.atoms[i].rings) > 1 :
                        if mapCyclesVoisins(c1.atoms[i],c2.atoms[j]) :
                            if mapVoisinsCyclesVoisins(c1.atoms[i],c2.atoms[j],c1,c2) :
                                continue
                            else :
                                return 0
                        else :
                            return 0
                    else :
                        continue
                else :
                    return 0
                continue
            else :
                return 0
        return 1
    else :
        return 0
                    

def mapCyclesVoisins(a1,a2) :
    i = j = tentativeRate = 0
    cyclesChecked = []
    while j < len(a1.rings) and tentativeRate < len(a1.rings)-len(cyclesChecked) :
        if i in cyclesChecked :
            i += 1
            continue
        if i == len(a1.rings) :
            i = 0
            continue
        if len(a1.rings[j].atoms) == len(a2.rings[i].atoms) :
            cyclesChecked.append(i)
            j += 1
            tentativeRate = 0
        else :
            tentativeRate += 1
        i += 1
    return j == len(a1.rings)

def mapVoisinsCyclesVoisins(a1,a2,c1,c2):
    for v in a1.oatoms :
        if c1 not in v.rings :
            v1 = v
            break
    else :
        return 0

    for w in a2.oatoms :
        if c2 not in w.rings :
            v2 = w
            break
    else :
        return 0

    if v1.number != v2.number :
        return 0

    if len(v1.rings) == len(v2.rings) :
        if len(v1.rings) > 1 :
            if mapCyclesVoisins(v1,v2) :
                return 1
            else :
                return 0
        else :
            return 1
    else :
        return 0

