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

from iMolecule import *
from frowns import MDL
import os

def readSmiles(file):
    """
    lit un fichier SMILES et renvoi les infos (liste des molecules, listes des noms, listes des id)
    """
    f=open(file)
    lines = f.readlines() # lit le contenu du fichier, chaque ligne est une chaine de caractere
    f.close()
    smilesError = "error.smi"
    sdfError = "error.sdf"
    mol=[]
    for line in lines : # pour chacune des lignes du fichier
        data = line.split() # decoupe la ligne en plusieurs chaines de caracteres pr pouvoir acceder a chacune des informations conenues ds la ligne
        if len(data)>1 :
            try :
                mol.append(iMolecule(smiles=data[0],id=data[1]))
            except :
                ferror=open(smilesError,"w")
                ferror.write(data[0]+"\t"+data[1]+"\n")
                ferror.close()
                os.system("babel -ismi "+smilesError+" -osdf "+sdfError)
                mol.extend(readSDF(sdfError))
        else :
            try :
                mol.append(iMolecule(smiles=data[0]))
            except :
                ferror=open(smilesError,"w")
                ferror.write(data[0]+"\n")
                ferror.close()
                os.system("babel -ismi "+smilesError+" -osdf "+sdfError)
                mol.extend(readSDF(sdfError))       
    return mol


def readSDF(file) :
    """
    idem ci-dessus pr fichier SDF
    """
    f = open(file)
    reader = MDL.sdin(f)
    mo = []
    try :
        m = reader.next()[0]
        while m :
            if m.fields.keys():
                for key in m.fields.keys():
                    if "ID" in key :
                        mo.append(iMolecule(sd=m,id=m.fields[key]))
                        break
            else :
                mo.append(iMolecule(sd=m,id=m.name))
            m = reader.next()[0]
            
    except :
        f.close()
        return mo

