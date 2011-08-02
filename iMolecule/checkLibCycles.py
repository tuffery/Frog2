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
sys.path.append("/home/tintin/tuffery/wrk/Frog/Frog/")
from SetMol import *
import os, tempfile

#The repertory's path of the cycles already created.
#If the option -m is used in the call of createLibCycles.sh, this path is not used
PATH_LIB_EXISTANT_CYCLES = "/media/disk/tuffery/wrk/Ambinter/librairie"

#The repertory's path of the cycles to be created. 
PATH_LIB_NEW_CYCLES = "/media/disk/tuffery/wrk/Ambinter/librairie2009/"

#The path of babel command
PATH_BABEL = "/home/tintin/tuffery/wrk/Frog/Frog/bin/babel"


def cycle(mol):
    """
    Only check if cycle there, do not write any mol2
    """
    try:
        global _nbCyclesEx
        global _listFiles
        """
        HERE IS A BUG !!!
        tmp file could have been overwritten ...
        """
        f = "temp.mol2"
        sdfFile = "temp.sdf"
        fd= tempfile.NamedTemporaryFile(suffix=".mol2",prefix="tmp",dir="./")
        sdfFile = fd.name
        # print sdfFile
        mol2File = open(f,"w")
        mol2File.write(mol)
        mol2File.close()
        os.system(PATH_BABEL + " -imol2 " + f + " -osdf " + sdfFile + " 2> babel_print.txt")
        os.system("echo '$$$$' >> " + sdfFile)
        # m = iMolecule(sdf="temp.sdf", extractCoordsFromSdf=True)
        m = iMolecule(sdf=sdfFile, extractCoordsFromSdf=True)
        fd.close()
        #try:
        #    m = 
        #    ff = open(sdfFile)
        #    reader = MDL.sdin(ff)
        #    m = reader.next()
        #    whil m:
        #        monoSdfFile = open("monoSdfTemp.sdf","w")
        #        monoSdfFile.write(m[2])
        #        monoSdfFile.close()
        #
        #
        #m = readSDF(sdfFile)
        if m:
            #for b in m[0].brickgraph :
            for b in m.brickgraph :
                cansmiles = b[0].cansmiles()
                if (sys.argv[1] == "False") and (cansmiles in _listFiles):
                    print cansmiles, "already done"
                    continue
                else:
                    print cansmiles, "not in library!!"
                    # b[0].write(file = PATH_LIB_NEW_CYCLES + cansmiles, mode="a")
                    _nbCyclesEx += 1
#                     if sys.argv[1] == "False":
#                         _listFiles.append(cansmiles)
    except:
            mode = "w"
            if os.path.isfile("./errors.txt"):
                mode = "a"
            error = open("./errors.txt",mode)
            error.write(mol[1])
            error.close()
            mol2File.close()
    sys.stdout.flush()

_listFiles = []
if sys.argv[1] == "False":
	_listFiles = os.listdir(PATH_LIB_EXISTANT_CYCLES)
        # print _listFiles
mol2banque = open(sys.argv[2],"r")
line = mol2banque.readline()
mol = ""
nbMolRead = 0
_nbCyclesEx = 0
while line != "":
	if mol and ("@<TRIPOS>MOLECULE" in line) :
		nbMolRead += 1
		cycle(mol)
		mol = ""+line
		line = mol2banque.readline()
	else:
		mol += line
		line = mol2banque.readline()
else:
	nbMolRead += 1
	cycle(mol)
print nbMolRead,"  molecules read"
print _nbCyclesEx," cycles extracted"
mol2banque.close()

