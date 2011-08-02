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


def cycle(mol, verbose = 0):
    cname = mol.split("\n")[1]
    try:
        global _nbCyclesEx
        global _listFiles
        # mol2File = "temp.mol2"
        sys.stdout.write("%s\n" % cname)
        mol2FileName = "%s.mol2" % cname
        sdfFileName = "%s.sdf" % cname
        
#        tempfile.NamedTemporaryFile(suffix=".mol2",prefix="tmp",dir="./")
#        mol2FileName = mol2File.name
        # sdfFile = "temp.sdf"
#        sdfFile = tempfile.NamedTemporaryFile(suffix=".sdf",prefix="tmp",dir="./")
#        sdfFileName = sdfFile.name
#         print mol2FileName, sdfFileName
        # mol2File = open(f,"w")
        mol2File = open(mol2FileName,"w")
        mol2File.write("%s" % mol)
        mol2File.close()
        # mol2File.close()
        cmd = "%s -imol2 %s -osdf %s 2> babel_print.txt" % (PATH_BABEL, mol2FileName, sdfFileName)
#         print cmd
        # os.system(PATH_BABEL + " -imol2 " + f + " -osdf " + sdfFile + " 2> babel_print.txt")
        os.system(cmd)
        os.system("echo '$$$$' >> " + sdfFileName)
        os.remove("%s" % mol2FileName)
        # m = iMolecule(sdf="temp.sdf", extractCoordsFromSdf=True)
    except:
            mode = "w"
            if os.path.isfile("./errors.txt"):
                mode = "a"
            error = open("./errors.txt",mode)
            error.write("Failed to convert mol2 to sdf\n")
            error.close()
            return
    try:
        m = iMolecule(sdf=sdfFileName, extractCoordsFromSdf=True)
        if m:
            #for b in m[0].brickgraph :
#             print "briques: %d\n" % len(m.brickgraph)
            # print m.brickgraph
            todo = {}
            for b in m.brickgraph :
                todo[b[0].handle] = 1
            for b in m.brickgraph :
                if not todo[b[0].handle]:
                    continue
                cansmiles = b[0].cansmiles()
                todo[b[0].handle] = 0
                if (sys.argv[1] == "False") and (cansmiles in _listFiles):
                    if verbose:
                        sys.stdout.write("Skip %s : already done ...\n" % (cansmiles))
                    continue
                else:
                    sys.stdout.write("Creating %s \n" % (cansmiles))
                    b[0].write(file = PATH_LIB_NEW_CYCLES + cansmiles, mode="a")
                    _nbCyclesEx += 1
#                     if sys.argv[1] == "False":
#                         _listFiles.append(cansmiles)
    except:
            mode = "w"
            if os.path.isfile("./errors.txt"):
                mode = "a"
            error = open("./errors.txt",mode)
            error.write(mol.split("\n")[1])
            error.close()
            mol2File.close()
    os.remove("%s" % sdfFileName)

"""
Cleanup by P. Tuffery, 2009:
- Will not scan same cycle more than once 
- Will not pass on cycle occurrences within same library

TODO: 
Should manage a count for occurrences, in order to pass from a minimal number of occurrences
"""
if __name__ == "__main__":
    _listFiles = []
    if sys.argv[1] == "False":
            _listFiles = os.listdir(PATH_LIB_EXISTANT_CYCLES)
            print "Read existing entries (%d)" % len(_listFiles)
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
                    # break
            else:
                    mol += line
                    line = mol2banque.readline()
    else:
        pass
            # nbMolRead += 1
            # cycle(mol)
    print nbMolRead,"  molecules read"
    print _nbCyclesEx," cycles extracted"
    mol2banque.close()

