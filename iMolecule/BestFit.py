#!/usr/local/bin/python
# -*- coding: utf-8 -*-
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

#
# HMM-SA tools upon the PDB class
#
# Written (2001-2004) by P. Tuffery, INSERM, France
#
# No warranty of any kind is provided
#
"""
La classe SACodec sert a encoder dans HMM-SA
Les quelques globales ci-dessous peuvent etre editees
selon la machine
"""

import string
import sys
import os
import copy
import math
import gzip
import types
import popen2
"""
from Fastav2 import *
from PDB6_5  import *
from PDBSA   import *

from SAVars import *
"""


from Config import *

#
# This is a simple wrapper of QBestFit
# but, it (will) can iterate to refine the fit.
#
class bestFit :#(PDBSA, PDB):

    def __init__(self, ncrd = 0, qxyz = [], txyz = [], method = "QBestFit", binPath = ANCILLARY_PATH, niter = 1, maxRMSd = 3.):

        self.binPath  = binPath
        self.method   = method

        self.niter    = niter
        self.maxRMSd  = maxRMSd
        self.orirmsd  = 0.
        self.bestrmsd = 0.
        self.ncrd     = ncrd
        self.qxyz     = qxyz
        self.txyz     = txyz
        self.ocrd     = []
        self.M1       = ""
        self.M2       = ""
        self.M3       = ""
        self.M4       = ""

    #
    # Call to best fit calculation
    # The two sets of corresponding crds are specified.
    # crd1 and crd2 are ncrd lines of 3 crds (x,y,z) for each point.
    # This will return : ori and best fit RMSd, the transformation matrix,
    # and if required, the crd1 crds superimpoed onto crd2.
    #
    def bestFit(self, getFit = 0):
        oristdout = sys.stdout
        cmd = self.binPath+"/"+self.method+" "
        if getFit:
            cmd += "-bfxyz "
        fin, sys.stdout = popen2.popen2(cmd)

        # We send data to method
        print self.ncrd
        for i in self.qxyz:
            print i[0],i[1],i[2]
        for i in self.txyz:
            print i[0],i[1],i[2]
        sys.stdout.close()
        sys.stdout = oristdout

        # we get the results
        self.orirmsd  = fin.readline(),
        self.bestrmsd = fin.readline(),
        self.M1 = fin.readline()
        self.M2 = fin.readline()
        self.M3 = fin.readline()
        self.M4 = fin.readline()
        
        self.ocrd = []
        if getFit:
            for i in range(0,self.ncrd):
                self.ocrd.append(fin.readline()[:-1])
        fin.close()

        # return the result
        return string.split(self.orirmsd[0])[1], string.split(self.bestrmsd[0])[1], self.M1, self.M2, self.M3, self.M4, self.ocrd

    #
    # transform the crds of one PDB instance
    #
def TM(x, ffrom, tto, TM, BINPATH=ANCILLARY_PATH):
    oristdout = sys.stdout
    fin, sys.stdout = popen2.popen2(BINPATH+"/TM")
    z = atmList(x[ffrom:tto].atms)
    crds = z.crds()
    # sys.stderr.write("TM: crds %s\n" % str(crds) )
    print len(z)
    for i in crds:
        print i
    print TM[0]
    print TM[1]
    print TM[2]
    print TM[3]
    sys.stdout.close()
    sys.stdout = oristdout
    
    orirmsd = string.split(orirmsd[0])[2]
    rmsd = string.split(rmsd[0])[1]
    print "REMARK  ",sys.argv[1]," (",len(x),")", sys.argv[2]," (",len(y),")"
    print "REMARK  ",from1, to1, from2, to2,"match len :",len(xs)
    print "REMARK   Initial RMSd  : ",orirmsd
    print "REMARK   Best fit RMSd : ",rmsd
    
    rs = []
    for i in range(0,len(z)):
        lrs= fin.readline()
        # sys.stderr.write("%s" % lrs)
        rs.append(z.list[i][0:30]+lrs[:-1]+z.list[i][54:-1])
    return rs

	
def TM2(atomList, M1, M2, M3, M4, BINPATH=ANCILLARY_PATH):
    oristdout = sys.stdout
    fin, sys.stdout = popen2.popen2(BINPATH+"/TM")
    print len(atomList)
    for i in range(len(atomList)):
        print atomList[i].x,atomList[i].y,atomList[i].z
    print M1
    print M2
    print M3
    print M4
    sys.stdout.close()
    sys.stdout = oristdout
    
##     orirmsd = string.split(orirmsd[0])[2]
##     rmsd = string.split(rmsd[0])[1]
    #rs = []
    for i in range(0,len(atomList)):
        lrs= fin.readline()
        rs = lrs.split()
        atomList[i].x = float(rs[0])
        atomList[i].y = float(rs[1])
        atomList[i].z = float(rs[2])
        """
        rs.append(lrs[:-1])
        # print rs
        # sys.exit(0)
    return rs
    """

