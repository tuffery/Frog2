/* ============================================================================
 *  This software is part of Frog, a chemo informatics class able to build 
 *  3D coordinates for small compounds
 *   Copyright (C) 2006-2007 P. Tuffery, B.O. Villoutreix, Th. Bohme Leite, D. Gomes, M. Miteva, J. Chomilier
 *
 *   Frog2 (C) 2009-2010 by P. Tuffery, M. Miteva, F. Guyon
 *
 *   Using this software, please cite:
 *       Frog2: Efficient 3D conformation ensemble generator for small compounds.
 *       Miteva MA, Guyon F, Tuffery P.
 *       Nucleic Acids Res. 2010 Jul;38(Web Server issue):W622-7. Epub 2010 May 5.
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *   ==========================================================================
 */

/* ###############################################################################
   #                                                                             #
   #         mol2clusterize.c                                                    #
   #                                                                             #
   #         authors P. Tuffery, F Guyon 2010                                    #
   #                                                                             #
   #         - input a mol2 that is a series of conformations of ONE compound    #
   #         atoms in the same order                                             #
   #         - perform clustering based on RMSd: not one output conformation     #
   #           should be close than RMSD threshold from another one              #
   #         - output filtered mol2 file (only diverse conformations)            #
   #                                                                             #
   #                                                                             #
   ############################################################################### */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>

#include "Types.h"
#include "argstr.h"
#include "TxtFile.h"
#include "Mol2.h"
#include "Zuker.h"
#include "RMSd.h"

/* ===============================================================
 * ---------- Main ----------------------------------------------- 
 * =============================================================== */
int main(int argc, char *argv[]) {

  FILE *oFd;
  int i, j, jj, k, kk, l, c, flag, posInVal=0, nbValRead, extremitesRead, numDiedreToRot, nbPasMC;
  char ligne[1025], val[1025];
  char buff[BUFSIZ], at0n[6],at1n[6],at2n[6],at3n[6];
  char confName[BUFSIZ];
  // double ***coordsConfs;      // coordonnees et energies
  // DtPoint3 **coordsConfs;
  double nrjTreshold;
  double nrjInitTreshold = 300.;
  DtStructConf *coordsConfs;
  DtMatrix4x4 M;

  DtPoint3 *newCoords;
  DtPoint3 *c1, *c2;

  int verbose = 0;

  double **rmsdVals;
  double drmsd;
  double ormsd;
  DtMolArray *pMa;
  int *flushp;
  int *refs;

  int aMol, aMol2;
  int refMol;
  int aAtom;
  int keepOn;

  int NOK;
  int nOut;

  srand(time(NULL));

  parseargstr(argc, argv);
  verbose = gVerbose;
  if (verbose)
    argstrPrmtrsSummary(argc, argv);
  
  if (IMol2FichName[0] != '\0') {
    pMa = lectMol2(IMol2FichName,  INT_MAX, verbose);
  }
  if (verbose) {
    fprintf(stderr,"Mol2Input done ... found %d mols\n", pMa->nMol);
  }

  /* Init RMSdARRAY */
  flushp    = calloc(pMa->nMol, sizeof(int));
  refs      = calloc(pMa->nMol, sizeof(int));
  rmsdVals = calloc(pMa->nMol, sizeof(double *));
  for (aMol = 0; aMol < pMa->nMol; aMol++) {
    rmsdVals[aMol] = calloc(pMa->nMol, sizeof(double));
  }
  

  /* Init workspace */
  c1 = calloc(pMa->pM[0].nAtoms, sizeof(DtPoint3));
  c2 = calloc(pMa->pM[0].nAtoms, sizeof(DtPoint3));

   if (verbose) {
    fprintf(stderr,"Memory allocation done ...\n");
  }

#if 0 
  for (aMol = 1; aMol < pMa->nMol; aMol++) {
    for (aAtom = 0; aAtom < pMa->pM[aMol].nAtoms; aAtom++) {
      fprintf(stderr,"mol %d: %.3lf %.3lf %.3lf\n", aMol, 
	      pMa->pM[aMol].atomDetails[aAtom].x,
	      pMa->pM[aMol].atomDetails[aAtom].y,
	      pMa->pM[aMol].atomDetails[aAtom].z);
	      
    }
  }
#endif

  for (aMol = 0; aMol < pMa->nMol; aMol++) {
    /* fill c1 */
    for (aAtom = 0; aAtom < pMa->pM[aMol].nAtoms; aAtom++) {
#if 0    
      memcpy(c1[aAtom], (DtPoint3 * ) &pMa->pM[aMol].atomDetails[aAtom].x, sizeof(DtPoint3));
#else
      c1[aAtom][0] = pMa->pM[aMol].atomDetails[aAtom].x;
      c1[aAtom][1] = pMa->pM[aMol].atomDetails[aAtom].y;
      c1[aAtom][2] = pMa->pM[aMol].atomDetails[aAtom].z;
#endif
      // fprintf(stderr,"mol1 (%d): %.3lf %.3lf %.3lf\n", aMol, c1[aAtom][0], c1[aAtom][1], c1[aAtom][2]);
    }
    // fprintf(stderr,"mol1 (%d): done ... \n", aMol);
    
    for (aMol2 = aMol+1; aMol2 < pMa->nMol; aMol2++) {

      /* fill c2 */
      // fprintf(stderr,"mol2 (%d %d atoms): \n", aMol2, pMa->pM[aMol2].nAtoms);
      for (aAtom = 0; aAtom < pMa->pM[aMol2].nAtoms; aAtom++) {
#if 0
	memcpy(c2[aAtom], (DtPoint3 * ) &pMa->pM[aMol2].atomDetails[aAtom].x, sizeof(DtPoint3));
#else
	c2[aAtom][0] = pMa->pM[aMol2].atomDetails[aAtom].x;
	c2[aAtom][1] = pMa->pM[aMol2].atomDetails[aAtom].y;
	c2[aAtom][2] = pMa->pM[aMol2].atomDetails[aAtom].z;
#endif
	// fprintf(stderr,"mol2 (%d) %.3lf %.3lf %.3lf\n", aMol2, c2[aAtom][0], c2[aAtom][1], c2[aAtom][2]);
      }
      
      /* fast RMSd ZUKER RETURNS SQUARED VALUES */
      frmsd(c1,c2, pMa->pM[aMol].nAtoms, &ormsd);
      drmsd = zuker_superpose(c1, c2, pMa->pM[aMol].nAtoms,  M);
      
      rmsdVals[aMol][aMol2] = drmsd;
      rmsdVals[aMol2][aMol] = rmsdVals[aMol][aMol2];
      if (verbose)
	fprintf(stderr, "%.3lf (%.3lf)", rmsdVals[aMol][aMol2], ormsd*ormsd);
    }
    if (verbose)
      fprintf(stderr, "\n");
  }


  /* clustering */
  for (aMol = 0; aMol < pMa->nMol; aMol++) {
    refs[aMol] = -1;
  }
  gRMSdTreshold *= gRMSdTreshold;  // Since ZUKER RETURNS DRMSd, not RMSd

  refMol = 0;
  keepOn = 0;
  do {
    refs[refMol] = refMol;
    if (verbose)
      fprintf(stderr,"refMol %d\n", refMol);
    for (aMol = 0; aMol < pMa->nMol; aMol++) {
      if (aMol == refMol) continue;
      if (flushp[aMol]) continue;
      if (rmsdVals[refMol][aMol] < gRMSdTreshold) {
	if (verbose)
	  fprintf(stderr,"%d tooClose from %d (%.3lf)\n",aMol, refMol, rmsdVals[refMol][aMol]);
	flushp[aMol] = 1;
	refs[aMol] = refMol;
      }
    }  
    keepOn = 0;
    aMol2 = -1;
    for (aMol = 0; aMol < pMa->nMol; aMol++) {
      if (refs[aMol] != -1) continue;
      keepOn = 1;
      if (aMol2 < 0)
	aMol2 = aMol;
      else {
	if (rmsdVals[refMol][aMol2]  < rmsdVals[refMol][aMol])
	  aMol2 = aMol;
      }
    }    
    refMol = aMol2;
  } while (keepOn);
  
  /* output */
  if (OMol2FichName[0] == '\000')
    oFd = stdout;
  else {
    oFd = fopen(OMol2FichName,"a");
    if (oFd == NULL) {
      if (verbose)
	fprintf(stderr,"Could not open %s\nReverting to stdout\n",OMol2FichName);
      oFd = stdout;
    }
  }

  for (aMol = 0; aMol < pMa->nMol; aMol++) {
    if (refs[aMol] == aMol) {
      if (verbose)
	fprintf(stderr,"%d is OK\n", aMol);
      outMol2(oFd,&pMa->pM[aMol], NULL, 1000000000., 1000000000., 0);
    }
  }

#if 0
  confNum = 0;
  // fprintf(stderr,"Will output %s\n", gId);
  while (confTmp != NULL) {
    if (confTmp->energie < bestConf->energie + nrjTreshold) {
      sprintf(confName,"%s_%d", gId, ++confNum);
      setupCoords(pMa->pM, confTmp->coords);
      outMol2(oFd,pMa->pM, confName, confTmp->energie, 1);
      confTmp = confTmp->nextConf;
    } else { break; }
  }
  fprintf(stdout,"%d %lf", confNum, bestConf->energie);
  if (verbose) fprintf(stderr, "%s: %d conformations generated\n", gId, confNum);
  // fprintf(stderr,"Memory cleanup 1\n");
#endif

  /* =============================================================
   * Memory cleanup
   * =============================================================
   */

  // fprintf(stderr,"Memory cleanup 2\n");
  if (verbose)
    fprintf(stderr,"Memory cleanup\n");
  free(rmsdVals);
  free(c1);
  free(c2);
  free(pMa->pM);
  free(pMa);

  return 0;
}


