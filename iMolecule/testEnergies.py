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


mOri = iMolecule("c1(C(=NO)Cc2ccccc2)c([n]([c]3[c]1cccc3)C)C")
mOri.graph2molecule()

fileSdf = open("testEnergyThiago1_3D_multiconf_bruno.sdf")
reader = MDL.sdin(fileSdf)
for i in range(21):
    m = reader.next()[0]

    for at in m.atoms:
        at.mmffAtomType = getMMFFAtomType(at)

    tmp = str(len(m.atoms)) + "\n500\n100\n"  
                
    for atom in m.atoms:
        vdwParams = VDW_PARAM[str(atom.mmffAtomType)]
        tmp += str(atom.x)+" "+str(atom.y)+" "+str(atom.z)+" "+str(getMMFFPartialAtomicCharge(atom))+" "+ \
               str(vdwParams[0])+" "+str(vdwParams[1])+" "+str(vdwParams[2])+" "+str(vdwParams[3])+"\n"
    i = 0
    for atom in m.atoms[:-1]:
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
        for j in range(i+1, len(m.atoms)):
            if j != i+1:
                tmp += " "
            if m.atoms[j] not in atomsJoined:
                tmp += "1"
            else:
                tmp += "0"
        tmp += "\n"
        i += 1
    tmp += "0\n"
    print tmp
    fi, fo = os.popen2("./monoconf","t")
    fi.write(tmp)
    fi.close()
    lines = fo.readlines()
    fo.close()
    print lines

