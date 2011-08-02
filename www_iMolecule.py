#!/usr/bin/python
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

import sys, os, random, time, tempfile
from iMolecule.iMolecule import *

from iMolecule.Config import *

"""

1er et 2eme argument:
    -ismi smile     (le smile doit être entre guillemets)
                    (2 arguments suivants optionnels: 
                        -axEq stringAE  
                     stringAE etant une chaine de caracteres de type sortie de toNconf())
    -isdf file.sdf  (file.sdf : chemin absolu
                     le fichier file.sdf sera détruit avant la fin du script)

arguments suivants:
    -toNconf ou -toNconfWAxEq: renvoie la liste des smiles non-ambigus, avec ou sans les chaines ax/eq

"""
def checkSmiles(fname):
    from frowns import Smiles
    try:
        fSmiles = open(fname, "r")
        lines = fSmiles.readlines()
        fSmiles.close()
        asmiles = Smiles.smilin(lines[0].split()[0])
        return True
    except:
        return False

def checkSDF(fname):
    import tempfile
    iFile3DSdf = tempfile.mktemp(".smi","tmpFrog","./")
    # print iFile3DSdf
    babel_convert(fname, "sdf", iFile3DSdf, "smi", verbose = 0)
    rs = checkSmiles(iFile3DSdf)
    os.remove(iFile3DSdf)
    return rs

def checkMol2(fname):
    import tempfile
    iFile3DSdf = tempfile.mktemp(".smi","tmpFrog","./")
    # print iFile3DSdf
    babel_convert(fname, "mol2", iFile3DSdf, "smi", verbose = 0)
    rs = checkSmiles(iFile3DSdf)
    os.remove(iFile3DSdf)
    return rs


def genPmlScript(cName, allStates = False):
    f = open("pymol.pml","w")
    f.write("load %s\n" % cName)
    f.write("bg_color white\n" )
    f.write("center\n" )
    f.write("hide lines\n" )
    f.write("show sticks\n" )
    f.write("show spheres\n" )
    f.write("set stick_radius=0.1\n" )
    f.write("set sphere_scale=0.15\n" )
    f.write("set antialias, 2\n" )
    f.write("set ray_trace_mode, 1\n" )
    f.write("ray 200,200\n")
    f.write("png %s_1.png\n" % cName)
    f.write("set all_states, on\n")
    f.write("ray 200,200\n")
    f.write("png %s_all.png\n" % cName)
    f.write("quit\n" )
    f.close()

def doSnap():
    cmd = "%s pymol.pml -c > /dev/null" % PYMOL
    # sys.stderr.write("%s\n" % cmd)
    os.system(cmd)



def usage():
    sys.stderr.write("www_iMolecule.py <options> <-ismi,-isdf,-imol2> inputfile\n")
    sys.stderr.write("Options:\n")
    sys.stderr.write(" -ismi file: SMILES input file\n")
    sys.stderr.write(" -isdf file: SDF input file\n")
    sys.stderr.write(" -imol2 file: mol2 input file\n")
    sys.stderr.write(" -i3Dsdf file: 3D SDF input file (no generation from scratch, just multiconf)\n")
    sys.stderr.write(" -i3Dmol2 file: 3D mol2 input file (no generation from scratch, just multiconf)\n")
    sys.stderr.write(" -opdb file: PDB ouput file\n")
    sys.stderr.write(" -omol2 file: mol2 ouput file\n")
    sys.stderr.write(" -osdf file: SDF ouput file\n\n")
    sys.stderr.write(" -osmi file: SMILES ouput file\n")
    sys.stderr.write(" -logFile file: log file\n")
    sys.stderr.write(" -ounsolved file: file to list compounds not generated\n")
    sys.stderr.write(" -wrkPath file: working directory\n\n")


    sys.stderr.write(" -mono: generate only one conformer per stereo-isomer.\n")
    sys.stderr.write(" -multi x: generate x conformers per stereo-isomer.\n")
    sys.stderr.write(" -unambiguate: enable search for undefined stereo centers.\n")
    sys.stderr.write(" -vb: enable stage two Monte-Carlo.\n")
    sys.stderr.write(" -mcsteps x: # steps for stage two Monte-Carlo (limit to small values since costly. 100 is usually OK).\n")
    sys.stderr.write(" -emax x: energy slice relative to lowest energy compound. 50, 70 OK since approximate fast energy)\n")
    sys.stderr.write(" -eini x: energy value to declare compound of high energy (500 OK).\n")
    sys.stderr.write(" -automini: minimize compounds of high energy on exit (reasonably slow).\n")
    sys.stderr.write(" -mini: minimize compounds on exit (very slow).\n")
    sys.stderr.write(" -miniEach: minimize compounds during generation (very very slow!).\n")
    sys.stderr.write(" -gnb: output multi conf over all isomers (and not for each isomer).\n")
    sys.stderr.write(" -rmsd x: output only compounds different by more than x Angstroms RMSd.\n")
    sys.stderr.write(" -maxMols x: limit gereation to maxMols first compounds.\n\n")
    sys.stderr.write(" -split: each compound in a different file.\n\n")

    sys.stderr.write(" -pymol: generate pymol script and produce compound image.\n")
    sys.stderr.write(" -v: (very) verbose mode.\n")

    
id = ""

if __name__ == "__main__":


    if len(sys.argv) < 2:
        usage()
        sys.exit(0)

    iFileSmiles   = None
    iFileSdf      = None
    iFile3DSdf    = None
    iFileMol2     = None
    iFile3DMol2   = None
    iRefFileName  = None
    multi         = None
    o3dFileName   = None
    oFileSmiles   = None
    o3dFileFormat = None
    logFile       = None
    wrkPathTmp    = None
    unsolvedFileName = None
    nbBestResults = 100

    wrkPathTmp  = "."
    sdfFile     = None
    sdfTmpFile  = None
    verbose     = False
    nStereos    = False
    sampleStereo = True
    maxMols     = 5000
    mcSteps     = 100
    eMax        = 100
    eIni        = 300
    minimize    = False
    autoMinimize = False
    miniEach    = False
    ammp_energy = False
    clusterize  = False
    clus_trhld  = 0.8
    split       = False
    vibrate     = False
    doPymol     = False
    globalNBest = False
    i = 1
    exceptions = ""
    try:
        while i < len(sys.argv):
	    if ( sys.argv[i] == "-h" ):
                usage()
                sys.exit(0)
	    if ( sys.argv[i] == "-ismi" ) and ( i+1 < len(sys.argv) ):
	        if not os.path.isfile(sys.argv[i+1]):
		    exceptions += "[ERROR] file for input smile doesn't exist\n"
		    #raise Exception("file for input smile doesn't exist")
	        iFileSmiles = sys.argv[i+1]
	        i += 2 
	    elif ( sys.argv[i] == "-isdf" ) and ( i+1 < len(sys.argv) ):
	        if not os.path.isfile(sys.argv[i+1]):
                    exceptions += "[ERROR] file for input sdf doesn't exist\n"
                    #raise Exception("file for input sdf doesn't exist")
	        iFileSdf = sys.argv[i+1]
	        i += 2
	    elif ( sys.argv[i] == "-i3Dsdf" ) and ( i+1 < len(sys.argv) ):
	        if not os.path.isfile(sys.argv[i+1]):
                    exceptions += "[ERROR] file for input 3D sdf doesn't exist\n"
                    #raise Exception("file for input sdf doesn't exist")
	        iFile3DSdf = sys.argv[i+1]
	        i += 2
	    elif ( sys.argv[i] == "-imol2" ) and ( i+1 < len(sys.argv) ):
	        if not os.path.isfile(sys.argv[i+1]):
                    exceptions += "[ERROR] file for input mol2 doesn't exist\n"
                    #raise Exception("file for input sdf doesn't exist")
	        iFileMol2 = sys.argv[i+1]
	        i += 2
	    elif ( sys.argv[i] == "-i3Dmol2" ) and ( i+1 < len(sys.argv) ):
	        if not os.path.isfile(sys.argv[i+1]):
                    exceptions += "[ERROR] file for input 3D mol2 doesn't exist\n"
                    #raise Exception("file for input sdf doesn't exist")
	        iFile3DMol2 = sys.argv[i+1]
	        i += 2
	    elif sys.argv[i] == "-unambiguate":
	        nStereos = True	
                sampleStereo = False
	        i += 1
	    elif ( sys.argv[i] == "-osmi" ) and ( i+1 < len(sys.argv) ):
	        oFileSmiles = sys.argv[i+1]
                i += 2
	    elif ( sys.argv[i] == "-logFile" ) and ( i+1 < len(sys.argv) ):
                logFile = sys.argv[i+1]
                i += 2
            elif ( sys.argv[i] == "-wrkPath" ) and ( i+1 < len(sys.argv) ):
                wrkPathTmp = sys.argv[i+1]
                i += 2
	    elif sys.argv[i] == "-mono":
                multi = False
                i += 1
	    elif sys.argv[i] == "-pymol":
                doPymol = True
                i += 1
	    elif sys.argv[i] == "-vb":
                vibrate = True
                i += 1
	    elif sys.argv[i] == "-v":
                verbose = True
                i += 1
	    elif sys.argv[i] == "-gnb":
                globalNBest = True
                i += 1
	    elif ( sys.argv[i] == "-multi" ) and ( i+1 < len(sys.argv) ):
                multi = True
	        try:
		    nbBestResults = int(sys.argv[i+1])
	        except:
                    exceptions += "[ERROR] a number of best results should be given after -multi\n"
                    i += 1
                    continue
                i += 2
            elif ( sys.argv[i] == "-maxMols" ) and ( i+1 < len(sys.argv) ):
	    	try:
                    maxMols = int(sys.argv[i+1])
                except:
                    exceptions += "[ERROR] a maximum number of molecules to compute should be given after -maxMols\n"
                    i += 1
                    continue
                i += 2
	    elif ( sys.argv[i] == "-iref" ) and ( i+1 < len(sys.argv) ):
	        iRefFileName = sys.argv[i+1]
	        i += 2
	    elif ( sys.argv[i] == "-opdb" ) and ( i+1 < len(sys.argv) ):
	        o3dFileName = sys.argv[i+1]
		o3dFileFormat = "pdb"
	        i += 2
            elif ( sys.argv[i] == "-omol2" ) and ( i+1 < len(sys.argv) ):
                o3dFileName = sys.argv[i+1]
                o3dFileFormat = "mol2"
                i += 2
	    elif ( sys.argv[i] == "-osdf" ) and ( i+1 < len(sys.argv) ):
                o3dFileName = sys.argv[i+1]
                o3dFileFormat = "sdf"
                i += 2
	    elif ( sys.argv[i] == "-ounsolved" ) and ( i+1 < len(sys.argv) ):
                unsolvedFileName = sys.argv[i+1]
                i += 2
            elif ( sys.argv[i] == "-emax" ) and ( i+1 < len(sys.argv) ):
                eMax = float(sys.argv[i+1])
                i += 2
            elif ( sys.argv[i] == "-eini" ) and ( i+1 < len(sys.argv) ):
                eIni = float(sys.argv[i+1])
                i += 2
            elif ( sys.argv[i] == "-mini" ):
                minimize = True
	        i += 1
            elif ( sys.argv[i] == "-automini" ):
                autoMinimize = True
	        i += 1
            elif ( sys.argv[i] == "-miniEach" ):
                miniEach = True
	        i += 1
            elif ( sys.argv[i] == "-ammpEne" ):
                ammp_energy = True
	        i += 1
            elif ( sys.argv[i] == "-rmsd" ):
                clusterize = True
                clus_trhld = float(sys.argv[i+1])
	        i += 2
            elif ( sys.argv[i] == "-split" ):
                split = True
	        i += 1
            elif ( sys.argv[i] == "-mcsteps" ) and ( i+1 < len(sys.argv) ):
                mcSteps = int(sys.argv[i+1])
                i += 2
	    else:
                exceptions += "[ERROR] unrecognized option for www_iMolecule\n"
                break
    except Exception, err:
	if logFile:
	    logF = open(logFile, "w")
	    logF.write("[ERROR] " + str(err) + "\n")
	    logF.close()
    

    if iFileSmiles:
        if not checkSmiles(iFileSmiles):
            exceptions += "Sorry: SOME PROBLEM OCCURRED. Is input data smiles ?\n"
    if iFile3DSdf:
        if not checkSDF(iFile3DSdf):
            exceptions += "Sorry: SOME PROBLEM OCCURRED. Is 3D input data SDF ?\n"
    if iFileSdf:
        if not checkSDF(iFileSdf):
            exceptions += "Sorry: SOME PROBLEM OCCURRED. Is input data SDF ?\n"

    if iFile3DMol2:
        if not checkMol2(iFile3DMol2):
            exceptions += "Sorry: SOME PROBLEM OCCURRED. Is 3D input data mol2 ?\n"

    if iFileMol2:
        if not checkMol2(iFileMol2):
            exceptions += "Sorry: SOME PROBLEM OCCURRED. Is input data mol2 ?\n"

    if exceptions != "":
        if logFile:
            logF = open(logFile, "w")
            logF.write(exceptions)
            logF.close()
        sys.exit(0)


#
# Mol2 -> sdf
#

    if iFile3DMol2:
        if iFile3DMol2.count(".mol2"):
            iFile3DSdf = iFile3DMol2.replace(".mol2",".sdf")
        else:
            iFile3DSdf = "%s.sdf" % iFile3DMol2
        babel_convert(iFile3DMol2, "mol2", iFile3DSdf, "sdf", verbose = 0)
    if iFileMol2:
        if iFileMol2.count(".mol2"):
            iFileSdf = iFileMol2.replace(".mol2",".sdf")
        else:
            iFileSdf = "%s.sdf" % iFileMol2
        babel_convert(iFileMol2, "mol2", iFileSdf, "sdf", verbose = 0)



#     """
#     Start of the big try
#     """	
    try:
        if wrkPathTmp and wrkPathTmp[-1] != "/":
            wrkPathTmp = wrkPathTmp + "/"
        if logFile:
	    logF = open(logFile, "w")
	if unsolvedFileName:
	    unsolvedFile = open(unsolvedFileName, "w")
	mols = []
	ids = []
	if iFileSmiles:
#             """
#             We have a SMILES INPUT
#             we maintain mols (compound SMILES) and ids (id from file or index if not)
#             """
            fSmiles = open(iFileSmiles, "r")
	    lines = fSmiles.readlines()
	    fSmiles.close()
            i = 1
	    for line in lines:
		lineSplit = line.split()
		if not lineSplit: # Empty line: go on
		    continue
                mols.append(lineSplit[0])
	 	if len(lineSplit) == 2:
		    ids.append(lineSplit[1])
		else:
		    ids.append("mol" + str(i))
                i += 1
	    if not mols:
                raise Exception("no input data detected")
            #sdfFile = open(wrkPathTmp+"tmpSmiForSdf_WWW.smi", "w")
            #for i in range(len(mols)):
            #    sdfFile.write(mols[i] + " " + ids[i] + "\n")
            #sdfFile.close()
            #os.system(BABEL_PATH + " -ismi " + wrkPathTmp + "tmpSmiForSdf_WWW.smi -osdf " +
            #          wrkPathTmp + "tmpSdfFromSmi_WWW.sdf")
            #os.system("echo '$$$$' >> " + wrkPathTmp + "tmpSdfFromSmi_WWW.sdf")
	elif iFileSdf or iFile3DSdf:
#             """
#             We have a SDF INPUT
#             we maintain mols (compound SMILES) and ids (id from file or index if not)
#             mols are the sdf line for each compound
#             we ouput a temporary SDF file, purged from void lines
#             """

	    tmpMol = []
            if iFile3DSdf:
                fSdf = open(iFile3DSdf, "r")
            else:
                fSdf = open(iFileSdf, "r")
	    lines = fSdf.readlines()
	    fSdf.close()
            i = len(lines) - 1
            while i > 0: # Skip void lines at end of file
		if lines[i] in ["", "\n"]:      
                    del lines[i]
		    i = i - 1
                else:
                    break	
	    idFound = -1
            for i, line in enumerate(lines):
#                 """
#                 $$$$ is the delimiter to finish compound
#                 """
		if line.find("$$$$") != -1:
                    tmpMol.append("$$$$")
		    mols.append(tmpMol)
		    if idFound != -1:
			ids.append(lines[idFound].split()[0])
		    else:
                        if tmpMol[0][-1] == "\n":
                            ids.append(tmpMol[0][:-1])
                        else:
                            ids.append(tmpMol[0])
                    tmpMol = []
		    idFound = -1
		elif line.find(">  <ID>") != -1 and i < len(lines)-1:
		    idFound = i+1
		    tmpMol.append(line)
		else:
		    tmpMol.append(line)
	    if lines[-1].find("M  END") != -1:
		tmpMol.append("$$$$")
                mols.append(tmpMol)
		if idFound != -1:
                    ids.append(lines[idFound].split()[0])
                else:
                    if tmpMol[0][-1] == "\n":
                        ids.append(tmpMol[0][:-1])
                    else:
                        ids.append(tmpMol[0])
	    if not mols:
		raise Exception("some problem occurred while readin the sdf (check its format)")
            # sdfFile = open(wrkPathTmp + "tmpSdf_WWW.sdf", "w")

            sdfTmpFile = tempfile.NamedTemporaryFile(suffix='.sdf', prefix="Frog2-", dir=wrkPathTmp)
            sdfTmpFileName = sdfTmpFile.name
            for i in range(len(mols)):
                sdfTmpFile.write("".join(mols[i]))
                sdfTmpFile.write("\n")
            # sdfFile.close()
            sdfTmpFile.flush()
	else:
	    raise Exception("no input file")
        if logFile:
            logF.write("disambiguation         : %5s" % str(nStereos))
            if multi:
                logF.write(",  nb max conformations   : %5s" % str(nbBestResults))
                logF.write(",  per isomer             : %s\n" % str(globalNBest))
            else:
                logF.write("\n")
            logF.write("energetic treshold     : %5s" % str(eMax) )
            logF.write(",  bad energy threshold   : %5s" % str(eIni))
            # logF.write("nb Monte Carlo steps   = " + str(mcSteps))
            logF.write(",  stage 2 Monte Carlo    : %s\n" % str(vibrate))
            if clusterize:
                logF.write("clustering threshold   : %5s,  " % str(clus_trhld))
            logF.write("minimize               : %5s\n" % str(minimize))
            logF.write("\n")
	    logF.flush()
        #if iFileSmiles:
        #    sdfFile = open(wrkPathTmp + "tmpSdfFromSmi_WWW.sdf", "r")
        #else:
        if iFileSdf or iFile3DSdf:
            # sdfFile = open(wrkPathTmp + "tmpSdf_WWW.sdf", "r")
            sdfFile = open(sdfTmpFileName, "r")
            if iFile3DSdf:
                reader = MDL.sdin(sdfFile, stripHydrogens=0)
            else:
                reader = MDL.sdin(sdfFile, stripHydrogens=1)
        
        if verbose:
            sys.stderr.write("Input of parameters and data done ... Will manage %d mol(s)\n" % len(mols))


	for i in range(len(mols)):
	    try:
                nbC = 0
            	if i == maxMols and logFile:
	   	    logF.write("maximum authorized number of molecules processed reached\n")
		    break
		id = ids[i]
                if not id:
                    id = "mol_" + str(i+1)    
	    	if logFile:
		    if id:
                    	logF.write("Id %d: "%(i+1) + id + "\n")
                    else:
                    	logF.write("Id %d: molecule "%(i+1) + str(i+1) + "\n")
			logF.flush()
                if iFileSmiles and nStereos:
                    if verbose:
                        sys.stderr.write("Will load mol %d as smiles\n" % i)
                    iMol = iMolecule(smiles = mols[i], id = id)
                    if verbose:
                        sys.stderr.write("Will graph2molecule mol %d\n" % i)
                    iMol.graph2molecule(wrkPathTmp)
                    if verbose:
                        sys.stderr.write("Did graph2molecule mol %d\n" % i)
                elif iFileSdf:
                    if verbose:
                        sys.stderr.write("Will load mol %d as sdf\n" % i)
                    theMol = reader.next()[0]
                    if verbose:
                        sys.stderr.write("Successfully loaded sdf (cansmiles %s)\n" % theMol.cansmiles())                    
                    iMol = iMolecule(smiles = theMol.cansmiles(), id = id, verbose = verbose)
                    if logFile and iMol.bigRingWarning:
                        logF.write("   [WARNING] Detected large rings that could be managed poorly by Frog (see limitations of Frog for bridged rings).\n")                        
                    if verbose:
                        sys.stderr.write("... success\n")
                elif iFile3DSdf:
                    if verbose:
                        sys.stderr.write("Will read 3D SDF (%s)\n" % id)
                    frowns_sdf = reader.next()[0]
                        # for at in frowns_sdf.atoms:
                        #     print at.x, at.y, at.z
                        # sys.exit(0)
                    iMol = iMolecule(sdfMolecule = frowns_sdf, id = id, extractCoordsFromSdf=True)
                    if verbose:
                        sys.stderr.write("Read 3D SDF ...\n")
                    smiles = [iMol.molecule.cansmiles()]
	            axsEqs = [None]
                    if oFileSmiles:
                        smilesF = open(oFileSmiles, "a")
                        smilesF.write(smiles[0] + "   " + id + "_1\n")
		        smilesF.flush()
                        smilesF.close()
                    if multi:
                        if verbose:
                            sys.stderr.write("Will generate multiconf from 3D SDF\n")
                        nbC += iMol.to1confWWW(fileName = o3dFileName, mode = "multi", minimize = minimize, autoMinimize = autoMinimize, ammp_energy = ammp_energy, clusterize = clusterize, clus_trhld = clus_trhld, format = o3dFileFormat, nbBestResults = nbBestResults, energeticBarrer = eMax, eIni = eIni, mcsteps=mcSteps, visuSdfFileName = wrkPathTmp + "/mol-%d.sdf"%(i+1), split = split, refFileName = iRefFileName, vibrate = vibrate, miniEach = miniEach, verbose = verbose)
                        if doPymol:
                            if mode == "mono":
                                genPmlScript("mol-%d.sdf"%(i+1))
                            else:
                                genPmlScript("mol-%d.sdf"%(i+1), allStates = True)
                            doSnap()
                        if logFile:
                            if iMol.eOK:
                                logF.write("   %s conformation(s) calculated.\n\n" % (str(nbC)) )
                                logF.flush()
                            else:
                                logF.write("   %s conformation(s) calculated but have bad energies.\n\n" % (str(nbC)) )
                                logF.flush()
                    else:
	                logF.write("   no new conformation calculated from 3D sdf\n\n")	
			logF.flush()
                    continue
#                 """
#                 Here, we are using SMILES in all cases but iFile3DSdf.
#                 iMol is a iMolecule instance from SMILES or 3D sdf.
#                 """
                #iMol.graph2molecule(wrkPathTmp)
                #iMol = iMolecule(smiles = reader.next()[0].cansmiles(), id = id)
                #iMol.graph2molecule(wrkPathTmp)
	    	#if iFileSmiles:
		#    iMol = iMolecule(smiles = reader.next()[0].cansmiles(), id = id)
	        #else:
                #    
		#    fileSdfTmp = open(wrkPathTmp+"tmpSdfWWW.sdf", "w")
		#    fileSdfTmp.write("".join(mols[i]))
                #    fileSdfTmp.close()
	        #    iMol = iMolecule(sdf = wrkPathTmp + "tmpSdfWWW.sdf", id = id)
                if nStereos:
#                     """
#                     Here, we disambiguate the compounds (SMILES LEVEL)
#                     and we also identify axial / equatorial sites
#                     """
                    if verbose:
                        sys.stderr.write("Will toNsmiles\n")
    	            smiles, axsEqs, removed = iMol.toNsmiles(True, verbose = verbose)
                    if verbose:
                        sys.stderr.write("Did toNsmiles\n")
                    n = len(smiles)*len(axsEqs)
                    if oFileSmiles:
                        smilesF = open(oFileSmiles, "a")
                        for j, smile in enumerate(smiles):
                            smilesF.write(smiles[j] + "   " + id + "_" + str(j) + "\n")
                            smilesF.flush()
                    if logFile:
                        if n > 8:
                            rndSmiles = 8
                            if rndSmiles > len(smiles):
                                rndSmiles = len(smiles)
                            logF.write("   " + str(n) + " isomer(s) found (see the smiles file)\n")
                            logF.write("   [WARNING] number of isomers considered restrainted to %d (random selection is performed):\n" % rndSmiles)
                            if removed:
                                logF.write("   [WARNING] the chirality of %d stereo centers specified in the input has been withdrawn (see limitations of Frog for bridged rings).\n" % removed)
                                
			    # logF.flush()
                            # for s in smiles:
                            #     logF.write("   %s\n" % s)
                            # for s in axsEqs:
                            #     logF.write("   %s\n" % s)
                            axsEqs = [None]
                            smiles = random.sample(smiles, rndSmiles)
                        else:
                            logF.write("   " + str(n) + " isomer(s) found:\n")
                            if removed:
                                logF.write("   [WARNING] the chirality of %d stereo center(s) specified in the input has been withdrawn (see limitations of Frog for bridged rings).\n" % removed)
			    logF.flush()
	        else:
                    """
                    We do not disambiguate. 
                    """
                    if iFileSmiles:
                        smiles = [mols[i]]
                    else:
                        smiles = [iMol.molecule.cansmiles()]
	            axsEqs = [None]
                    if oFileSmiles:
                        smilesF = open(oFileSmiles, "a")
                        smilesF.write(smiles[0] + "   " + id + "_1\n")
		        smilesF.flush()


                if verbose:
                    sys.stderr.write("Initial smiles generation done ...\n")

                #j = 1
                for smile in smiles:
                    for axEq in axsEqs:
       	                if axEq and (axEq.find("A") != -1 or axEq.find("E") != -1) and logFile:
		            logF.write("       " + smile + "   " + axEq + "\n")
			    logF.flush()
		        elif logFile:
		            logF.write("      " + smile + "\n")
			    logF.flush()
		    #if oFileSmiles:
                    #    smilesF.write(smile + "   " + id + "_" + str(j) + "\n")
		    #    smilesF.flush()
		    #j += 1
                if oFileSmiles:
                    smilesF.close()
	        if multi == None:
                    if logFile:
	                if o3dFileName:
	                    logF.write("\ncannot write some conformations without the option -multi or -mono\n")
	                    logF.close()
                            break
	            continue
                
                if verbose:
                    sys.stderr.write("Will perform 3D generation ... sampleStereo is %d\n" % int(sampleStereo))

                j = 1
                nbC = 0
                highEnbC = 0
#                 """
#                 Now we gererate 3D conformations for each isomer.
#                 """
                if (globalNBest):
                    if len(smiles) > 1:
                        lnbBestResults = (nbBestResults / len(smiles) ) + 1
                    else:
                        lnbBestResults = nbBestResults
                else:
                    lnbBestResults = nbBestResults
                for smile in smiles:
	            for axEq in axsEqs:
		        iMol = iMolecule(smiles = smile,  id = id+"_"+str(j), axialEquatorial = axEq, sampleStereo = sampleStereo, verbose = verbose)
                        
		        #iMol.graph2molecule(wrkPathTmp)
                        j += 1 

                        """
                        Frog2.0:
                        We always generate 1 conformation, and additional file
                        """
                        # if verbose:
                        #     sys.stderr.write("Will perform 3D conformer generation ...\n")
                        # iMol.to1confWWW(fileName = o3dFileName, format = o3dFileFormat, nbBestResults = nbBestResults, energeticBarrer = eMax, mcsteps=mcSteps, visuSdfFileName = wrkPathTmp + "/mol-%d.sdf"%(i+1), verbose = verbose)
		        # if multi:
                        #     if verbose:
                        #         sys.stderr.write("Will perform 3D multiconformer generation ...\n")
		        #     nbC += iMol.toNconfWWW(fileName = o3dFileName, format = o3dFileFormat, nbBestResults = nbBestResults,
                        #                            energeticBarrer = eMax, mcsteps=mcSteps, visuSdfFileName = wrkPathTmp + "/mol-%d.sdf"%(i+1), verbose = verbose)

                        if verbose:
                            sys.stderr.write("Will perform 3D conformer generation for %s...\n" % id)
                        mode = "mono"
                        if multi:
                            mode = "multi"
                        lnbC = iMol.to1confWWW(fileName = o3dFileName, mode = mode, minimize = minimize, autoMinimize = autoMinimize, ammp_energy = ammp_energy, clusterize = clusterize, clus_trhld = clus_trhld, format = o3dFileFormat, nbBestResults = lnbBestResults, energeticBarrer = eMax, eIni = eIni, mcsteps=mcSteps, visuSdfFileName = wrkPathTmp + "/mol-%d.sdf"%(i+1), split = split, refFileName = iRefFileName, vibrate = vibrate, miniEach = miniEach, verbose = verbose)
                        if doPymol:
                            if mode == "mono":
                                genPmlScript("mol-%d.sdf"%(i+1))
                            else:
                                genPmlScript("mol-%d.sdf"%(i+1), allStates = True)
                            doSnap()
                        if verbose:
                            sys.stderr.write("Got %d confs, OK status: %s\n" % (lnbC, iMol.eOK))
                        if iMol.eOK:
                            nbC += lnbC
                        else:
                            highEnbC += lnbC

                if verbose:
                    sys.stderr.write("3D generation done ...\n")

	        if multi and logFile:
	            if nStereos:
	     	        logF.write("   %s conformation(s) calculated. For some isomer(s) %d have bad energies.\n\n" % (str(nbC+highEnbC), highEnbC) )
			logF.flush()
	            else:
		        logF.write("   %s conformation(s) calculated. For some isomer(s) %d have bad energies.\n\n" % (str(nbC+highEnbC), highEnbC) )
			logF.flush()
	        else:
	            if nStereos:
	                logF.write("   1 conformation(s) calculated for each isomer. %d have bad energies.\n\n" % highEnbC)	
			logF.flush()
	            else:
                        logF.write("   1 conformation(s) calculated\n\n")
			logF.flush()
                
            except SystemExit:
                pass
	    except MissingParameter, err:
                if logFile:
                    logF.write("  [ERROR] some energy parameters may be missing (is the molecule ADME/tox compliant?)\n")
		    logF.flush()
                if unsolvedFileName:
                    unsolvedFile.write(""+id+"\n")
		    unsolvedFile.flush()
            except TypeError, err:
                logF.write("  [ERROR] an internal problem occured. Please try again this compound\n")
                logF.flush()
            except Exception, err:
                print err
                if logFile:
                    logF.write("  [ERROR] " + str(err) + "\n")
		    logF.flush()
		if unsolvedFileName:
		    unsolvedFile.write(""+id+"\n")
		    unsolvedFile.flush()
	    else:
		pass

#     """
#     End of the big try
#     """
    except SystemExit:
        pass
    except Exception, err:
        print err
	if logFile:
            logF.write("  [ERROR] " + str(err) + "\n")
	    logF.flush()
    if logFile:
        try:
            logF.close()
	except:
            pass
    if unsolvedFileName:
        try:
	    unsolvedFile.close()
        except:
            pass
    try:
        sdfFile.close()
    except:
        pass
    try:
        sdfTmpFile.close()
    except:
        pass




