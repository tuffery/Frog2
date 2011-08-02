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

"""
These procedures are to interface Frog with its specific version of AMMOS
"""

import sys
import os
from Config import *

def AMMOS_Ring_Generate(name2, LIBRARY_PATH, BABEL_PATH, verbose = 0):
    
    # Generer du mol2.
    uName = name2.replace("(","\(").replace(")","\)").replace("[","\[").replace("]","\]")
    # f  =open("%s/%s.smi" % (LIBRARY_PATH, name2), "w")
    f  =open("%s/%s.smi" % (LIBRARY_PATH, uName), "w")
    f.write("%s\n" %  name2)
    f.close()
    cmd = "%s -h -ismi %s/%s.smi -omol2 %s/%s.mol2" % (BABEL_PATH, LIBRARY_PATH, uName, LIBRARY_PATH, uName)
    # cmd = "%s -h -ismi %s/%s.smi -omol2 %s/%s.mol2" % (BABEL_PATH, LIBRARY_PATH, "itest", LIBRARY_PATH, "itest")
    if verbose:
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
    if verbose:
        sys.stderr.write("Generated %s/%s.dgammos\n" % (LIBRARY_PATH, uName))
        sys.stderr.write("Generating %s/%s.mol2\n" % (LIBRARY_PATH, uName))
    cmd = "%s %s/%s.dgammos >& /dev/null" % (AMMOSBUILD, LIBRARY_PATH, uName)
    # cmd = "%s %s.dgammos" % (AMMOSBUILD, uName)
    sys.stderr.write("DG_AMMOS cmd:  %s\n" % (cmd))
    if verbose:
        os.system(cmd)
    # sys.exit(0)

    # Convertir en sdf, oter les H
    # cmd = "%s -d -imol2 %s/%s_Built.mol2 -osdf %s/%s.sdf 2> %s/%s.log" % (BABEL_PATH, LIBRARY_PATH, name2, LIBRARY_PATH, name2, LIBRARY_PATH, name2)
    cmd = "%s -d -imol2 %s/%s_Built.mol2 -osdf %s/%s 2> %s/%s.log" % (BABEL_PATH, LIBRARY_PATH, uName, LIBRARY_PATH, uName, LIBRARY_PATH, uName)
    os.system(cmd)

def AMMOS_Generate(name2, LIBRARY_PATH, BABEL_PATH, verbose = 0):
    """
    Generic function to generate compound from scratch
    """
    
    # Generer du mol2.
    uName = name2.replace("(","\(").replace(")","\)").replace("[","\[").replace("]","\]")
    # f  =open("%s/%s.smi" % (LIBRARY_PATH, name2), "w")
    f  =open("%s/%s.smi" % (LIBRARY_PATH, uName), "w")
    f.write("%s\n" %  name2)
    f.close()
    cmd = "%s -h -ismi %s/%s.smi -omol2 %s/%s.mol2" % (BABEL_PATH, LIBRARY_PATH, uName, LIBRARY_PATH, uName)
    # cmd = "%s -h -ismi %s/%s.smi -omol2 %s/%s.mol2" % (BABEL_PATH, LIBRARY_PATH, "itest", LIBRARY_PATH, "itest")
    if verbose:
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
    if verbose:
        sys.stderr.write("Generated %s/%s.dgammos\n" % (LIBRARY_PATH, uName))
        sys.stderr.write("Generating %s/%s.mol2\n" % (LIBRARY_PATH, uName))
    cmd = "%s %s/%s.dgammos >& /dev/null" % (AMMOSBUILD, LIBRARY_PATH, uName)
    # cmd = "%s %s.dgammos" % (AMMOSBUILD, uName)
    if verbose:
        sys.stderr.write("DG_AMMOS cmd:  %s\n" % (cmd))
    os.system(cmd)
    # sys.exit(0)

    # Convertir en sdf, oter les H
    # cmd = "%s -d -imol2 %s/%s_Built.mol2 -osdf %s/%s.sdf 2> %s/%s.log" % (BABEL_PATH, LIBRARY_PATH, name2, LIBRARY_PATH, name2, LIBRARY_PATH, name2)
    cmd = "%s -d -imol2 %s/%s_Built.mol2 -osdf %s/%s 2> %s/%s.log" % (BABEL_PATH, LIBRARY_PATH, uName, LIBRARY_PATH, uName, LIBRARY_PATH, uName)
    os.system(cmd)

def AMMOS_Minimize(fName, verbose = 0):
    """
    Here, we suppose fName is already mol2 format
    @param fName: input file (must end without .mol2, BUT fName.mol2 musst exist)

    @return fName_minimized.mol2
    """
    
    # Generer du mol2.
    uName = fName.replace("(","\(").replace(")","\)").replace("[","\[").replace("]","\]")
    if verbose:
        sys.stderr.write("AMMOS_Minimize: fName %s uName %s\n" % (fName, uName))
    # cmd = "%s -h -i%s %s -omol2 %s.mol2" % (BABEL_PATH, o3dFileFormat, uName, uName)
    # sys.stderr.write("%s\n" % cmd)
    # status = os.system(cmd)

    # Generer le fichier input.param
    # path_of_DG-AMMOS= /Users/DG-AMMOS
    # bank= input_dataset.mol2
    f=  open("%s.ammos" % (fName), "w")
    f.write("path_of_AMMOS_SmallMol= %s\n" % AMMOSMINIHOME)
    f.write("bank= %s\n" % (fName))
    f.close()
    if verbose:
        sys.stderr.write("AMMOS_Minimize: Generated %s.ammos\n" % (uName))
        sys.stderr.write("AMMOS_Minimize: Generating %s_minimized.mol2\n" % (uName))
    cmd = "%s %s.ammos >& /dev/null" % (AMMOSMINI, uName)
    if verbose:
        sys.stderr.write("AMMOS cmd:  %s\n" % (cmd))
    status = os.system(cmd)
    # sys.exit(0)

    # Convertir au format final, sans oter les H
    # cmd = "%s -imol2 %s_minimized.mol2 -o%s %s 2> %s.log" % (BABEL_PATH, uName, o3dFileFormat, uName, uName)
    # os.system(cmd)

def AMMOS_QuickMinimize(fName, verbose = 0):
    """
    Here, we suppose fName is already mol2 format
    @param fName: input file (must end without .mol2, BUT fName.mol2 musst exist)

    @return fName_minimized.mol2
    """
    
    # Generer du mol2.
    uName = fName.replace("(","\(").replace(")","\)").replace("[","\[").replace("]","\]")
    if verbose:
        sys.stderr.write("AMMOS_Minimize: fName %s uName %s\n" % (fName, uName))
    # cmd = "%s -h -i%s %s -omol2 %s.mol2" % (BABEL_PATH, o3dFileFormat, uName, uName)
    # sys.stderr.write("%s\n" % cmd)
    # status = os.system(cmd)

    # Generer le fichier input.param
    # path_of_DG-AMMOS= /Users/DG-AMMOS
    # bank= input_dataset.mol2
    f=  open("%s.ammos" % (fName), "w")
    f.write("path_of_AMMOS_SmallMol= %s\n" % AMMOSMINIHOME)
    f.write("bank= %s\n" % (fName))
    f.close()
    if verbose:
        sys.stderr.write("AMMOS_Minimize: Generated %s.ammos\n" % (uName))
        sys.stderr.write("AMMOS_Minimize: Generating %s_minimized.mol2\n" % (uName))
    cmd = "%s %s.ammos >& /dev/null" % (AMMOSQUICKMINI, uName)
    if verbose:
        sys.stderr.write("AMMOS cmd:  %s\n" % (cmd))
    status = os.system(cmd)
    # sys.exit(0)

    # Convertir au format final, sans oter les H
    # cmd = "%s -imol2 %s_minimized.mol2 -o%s %s 2> %s.log" % (BABEL_PATH, uName, o3dFileFormat, uName, uName)
    # os.system(cmd)

def AMMOS_Energy(fName, verbose = 0):
    """
    Here, we suppose fName is already mol2 format
    @param fName: input file (must end without .mol2, BUT fName.mol2 musst exist)

    @return fName_energy.mol2
    """
    
    # Generer du mol2.
    uName = fName.replace("(","\(").replace(")","\)").replace("[","\[").replace("]","\]")
    if verbose:
        sys.stderr.write("AMMOS_Energy: fName %s uName %s\n" % (fName, uName))
    # cmd = "%s -h -i%s %s -omol2 %s.mol2" % (BABEL_PATH, o3dFileFormat, uName, uName)
    # sys.stderr.write("%s\n" % cmd)
    # status = os.system(cmd)

    # Generer le fichier input.param
    # path_of_DG-AMMOS= /Users/DG-AMMOS
    # bank= input_dataset.mol2
    f=  open("%s.eammos" % (fName), "w")
    f.write("path_of_AMMOS_SmallMol= %s\n" % AMMOSMINIHOME)
    f.write("bank= %s\n" % (fName))
    f.close()
    if verbose:
        sys.stderr.write("AMMOS_Energy: Generated %s.eammos\n" % (uName))
        sys.stderr.write("AMMOS_Energy: Generating %s_energy.mol2\n" % (uName))
    cmd = "%s %s.eammos >& /dev/null" % (AMMOSENE, uName)
    if verbose:
        sys.stderr.write("AMMOS cmd:  %s\n" % (cmd))
    status = os.system(cmd)
    # sys.exit(0)

    # Convertir au format final, sans oter les H
    # cmd = "%s -imol2 %s_minimized.mol2 -o%s %s 2> %s.log" % (BABEL_PATH, uName, o3dFileFormat, uName, uName)
    # os.system(cmd)


