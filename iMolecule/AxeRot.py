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
import popen2

from Config import *
method = "AxeRot"

class AxeRot :

    def __init__(self, at1, at2, angle) :
        self.at1 = at1
        self.at2 = at2
        self.angle = angle
        self.M1 = ""
        self.M2 = ""
        self.M3 = ""
        self.M4 = ""

    def axeRot(self) :
        oristdout = sys.stdout
        cmd = ANCILLARY_PATH+"/"+method
        fin, sys.stdout = popen2.popen2(cmd)

        print self.at1.x, self.at1.y, self.at1.z
        print self.at2.x, self.at2.y, self.at2.z
        print self.angle
        sys.stdout.close()
        sys.stdout = oristdout

        self.M1 = fin.readline()
        self.M2 = fin.readline()
        self.M3 = fin.readline()
        self.M4 = fin.readline()

def TM(atomList, M1, M2, M3, M4) :
    oristdout = sys.stdout
    fin, sys.stdout = popen2.popen2(ANCILLARY_PATH+"/TM")
    print len(atomList)
    for a in atomList :
        print a.x, a.y, a.z
    print M1
    print M2
    print M3
    print M4
    sys.stdout.close()
    sys.stdout = oristdout
    
    for a in atomList :
        lrs = fin.readline()
        #rs = lrs.split()
        rs = [lrs[0:8].strip(), lrs[8:16].strip(), lrs[16:24].strip()]
        a.x = float(rs[0])
        a.y = float(rs[1])
        a.z = float(rs[2])

