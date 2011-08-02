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

import math

class Diedre :

    def __init__(self,atoms,movingList) :
        self.atoms = atoms
        self.value = self.computeValue()
        self.movingList = movingList

    def computeValue(self) :
        """
        calcul de la valeur de l'angle diedre
        """
        a = self.atoms[0]
        b = self.atoms[1]
        c = self.atoms[2]
        d = self.atoms[3]
        #print self.atoms
        ab = [0,0,0]
        bc = [0,0,0]
        cd = [0,0,0]
        
        ab[0] = b.x - a.x
        ab[1] = b.y - a.y
        ab[2] = b.z - a.z
        bc[0] = c.x - b.x
        bc[1] = c.y - b.y
        bc[2] = c.z - b.z
        cd[0] = d.x - c.x
        cd[1] = d.y - c.y
        cd[2] = d.z - c.z
        
        d012 = ab[0]*bc[0] + ab[1]*bc[1] + ab[2]*bc[2]
        d123 = cd[0]*bc[0] + cd[1]*bc[1] + cd[2]*bc[2]
        d0123 = ab[0]*cd[0] + ab[1]*cd[1] + ab[2]*cd[2]
        
        d01  = ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]
        d12  = bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2]
        d23  = cd[0]*cd[0] + cd[1]*cd[1] + cd[2]*cd[2]
        
        num = d012 * d123 - d12 * d0123
        den = (d01*d12 - d012*d012)*(d12*d23 - d123*d123)
        #for at in self.atoms:
        #    print at.x, at.y,at.z
        #print d01*d12, d012*d012,d12*d23,d123*d123
        arcos = num / math.sqrt(den)
        
        if arcos > 1. :
            arcos = 1.
            
        if arcos < -1. :
            arcos = -1
            
        RS = math.acos(arcos)
        
        RS1 = cd[0] * (ab[1] * bc[2] - ab[2] * bc[1]) + cd[1] * (bc[0] * ab[2] - ab[0] * bc[2]) + cd[2] * (ab[0] * bc[1] - ab[1] * bc[0])

        if RS1 > 0. :
            return RS
        else :
            return -RS

    def setValue(self,v) :
        self.value = v

