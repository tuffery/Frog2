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

from SetMol import *
import os

def cycle(mol,sdfFile,listFile,f) :
    try :
        mol2File = open(f,"w")
        mol2File.write(mol)
        mol2File.close()
        os.system("./bin/babel -imol2 "+f+" -osdf "+sdfFile)
        m = readSDF(sdfFile)
        if m :
            for b in m[0].brickgraph :
                if b[0].cansmiles() in listFile :
                    continue
                else :
                    b[0].write(file=PATH+b[0].cansmiles())
                    listFile = os.listdir(PATH)
    except :
        print "error"
        return
        
f = "temp.mol2"
sdfFile = "temp.sdf"
listFile = os.listdir(PATH)
banque = ["/home/prot1/gomes/python/banque_sdf/nci/NCI_aug00_filtered_3D.mol2.gz",
          "/home/prot1/gomes/python/banque_sdf/asinex/AsinexGoldCollection_Feb2004_filtered_3D.mol2.gz",
          "/home/prot1/gomes/python/banque_sdf/asinex/AsinexPlatinumCollection_Feb2004_filtered_3D.mol2.gz",
          "/home/prot1/gomes/python/banque_sdf/chembridge/diverset_march2004_filtered_3D.mol2.gz",
          "/home/prot1/gomes/python/banque_sdf/chembridge/express_pick_march2004_part1_filtered3D.mol2.gz",
          "/home/prot1/gomes/python/banque_sdf/chembridge/express_pick_march2004_part2_filtered3D.mol2.gz",
          "/home/prot1/gomes/python/banque_sdf/chembridge/hit2lleadscreening_march2004_part1_filtered3D.mol2.gz",
          ]
for file in banque :
    print "decopression",file,"en cours"
    os.system("gunzip "+file)
    mol = ""
    mol2banque = open(file[:-3])
    line = mol2banque.readline()
    i = 0
    
    while line != "" :
        if mol and ("@<TRIPOS>MOLECULE" in line) :
            i += 1
            cycle(mol,sdfFile,listFile,f)
            print i
            mol = ""+line
            line = mol2banque.readline()
        else :
            mol += line
            line = mol2banque.readline()
    else :
        cycle(mol,sdfFile,listFile,f)
    mol2banque.close()
    os.remove(file[:-3])
    res = open("/home/prot1/gomes/python/banque_sdf/banques_faites","a")
    res.write(file[:-3]+"\n")
    res.close()
    print file[:-3],"OK !!!!"

