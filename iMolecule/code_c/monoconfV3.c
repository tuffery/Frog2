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
   #         monoconfV3.c                                                        #
   #                                                                             #
   #         authors P. Tuffery, T. Bohme Leite                                  #
   #                                                                             #
   #         from v3.0 (09/2009, P. Tuffery):                                    #
   #                    accepts command line arguments instead of just stdin     #
   #                                                                             #
   ############################################################################### */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "Types.h"
#include "argstr.h"
#include "TxtFile.h"
#include "Mol2.h"

#include "Zuker.h"
#include "RMSd.h"

/* Calcul d'un angle dihedral ----------------------------------------- 
 *
 * a,b,c,d: Pointeurs sur coordonnees 3D (DtFloats)
 */
DtFloat dihedral(DtFloat a[3],DtFloat b[3],DtFloat c[3],DtFloat d[3], int asDegrees)
{
  DtPoint3 ab, bc, cd;
  DtFloat d01,d012,d12,d123,d23,d0123,RS,RS1;
  double num,den, arccos;

  ab[0] = (b[0] - a[0]);
  ab[1] = (b[1] - a[1]);
  ab[2] = (b[2] - a[2]);
  bc[0] = (c[0] - b[0]);
  bc[1] = (c[1] - b[1]);
  bc[2] = (c[2] - b[2]);
  cd[0] = (d[0] - c[0]);
  cd[1] = (d[1] - c[1]);
  cd[2] = (d[2] - c[2]);
  
  d012 = ab[0]*bc[0] + ab[1]*bc[1] + ab[2]*bc[2];
  d123 = cd[0]*bc[0] + cd[1]*bc[1] + cd[2]*bc[2];
  d0123 = ab[0]*cd[0] + ab[1]*cd[1] + ab[2]*cd[2];
  
  d01  = ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2];
  d12  = bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2];
  d23  = cd[0]*cd[0] + cd[1]*cd[1] + cd[2]*cd[2];
  
  num = (double)(d012 * d123 - d12 * d0123);
  den = (double)(d01*d12 - d012*d012)*(d12*d23 - d123*d123);
  arccos = num / sqrt(den);
  
  if (arccos > (double) 1.) { /* To make it robust !!! */
    arccos = 1.;
  }
  if (arccos < (double) -1.) {
    arccos = -1.;
  } 
  RS = (DtFloat) acos(arccos);
  if (asDegrees) RS = RTOD(RS);
  
  RS1 =   cd[0] * (ab[1] * bc[2] - ab[2] * bc[1]) +
    cd[1] * (bc[0] * ab[2] - ab[0] * bc[2]) +
    cd[2] * (ab[0] * bc[1] - ab[1] * bc[0]);
  
  if (RS1 > 0.) return RS;
  else return - RS;
}


/* ===============================================================
 * Matrice de rotation autour d'un axe a b. 
 * PREMULTIPLIE -> Y = XM (X vecteur ligne) 
 * angle is RADIAN
 * =============================================================== */
void MkArbitraryAxisRotMat4x4(DtPoint3 A,DtPoint3 B,DtFloat angle,DtMatrix4x4 M) {
  extern double sqrt(), sin(), cos();

  DtFloat rx, ry, rz, lnorme, rx2, rxry, rxrz, ry2, ryrz, rz2, 
  rxry1mcosa, rxrz1mcosa, ryrz1mcosa;

  DtFloat cosa, sina;

  cosa = cos((double) angle);
  sina = sin((double) angle);

  rx = B[0] - A[0];
  ry = B[1] - A[1];
  rz = B[2] - A[2];
  lnorme = (DtFloat) sqrt(rx*rx + ry*ry + rz*rz);
  rx /= lnorme;
  ry /= lnorme;
  rz /= lnorme;

  rx2 = rx*rx;
  ry2 = ry*ry;
  rz2 = rz*rz;
  rxry = rx*ry;
  rxrz = rx*rz;
  ryrz = ry*rz;
  rxry1mcosa = rxry * (1. - cosa);
  rxrz1mcosa = rxrz * (1. - cosa);
  ryrz1mcosa = ryrz * (1. - cosa);

  M[0][0] = rx2 + (1. - rx2) * cosa;
  M[1][0] = rxry1mcosa - rz * sina;
  M[2][0] = rxrz1mcosa + ry * sina;
  
  M[0][1] = rxry1mcosa + rz * sina;
  M[1][1] = ry2 + (1. - ry2) * cosa;
  M[2][1] = ryrz1mcosa - rx * sina;

  M[0][2] = rxrz1mcosa - ry * sina;
  M[1][2] = ryrz1mcosa + rx * sina;
  M[2][2] = rz2 + (1. - rz2) * cosa;

  M[3][0] =  A[0] * (1 - M[0][0]) - A[1] * M[1][0] - A[2] * M[2][0];
  M[3][1] = -A[0] * M[0][1] + A[1] * (1 - M[1][1]) - A[2] * M[2][1];
  M[3][2] = -A[0] * M[0][2] - A[1] * M[1][2] + A[2] * (1 - M[2][2]);

  M[0][3] = M[1][3] = M[2][3] = 0.;
  M[3][3] = 1.;
}


/* ===============================================================
 * Rotation un atome par rmat. 
 * =============================================================== */
void singleRotate(DtPoint3 p,DtMatrix4x4 rmat) {
  double x,y,z;

  x = p[0] * rmat[0][0] + p[1] * rmat[1][0] + p[2] * rmat[2][0] + rmat[3][0]; 
  y = p[0] * rmat[0][1] + p[1] * rmat[1][1] + p[2] * rmat[2][1] + rmat[3][1]; 
  z = p[0] * rmat[0][2] + p[1] * rmat[1][2] + p[2] * rmat[2][2] + rmat[3][2]; 
  p[0] = x;
  p[1] = y;
  p[2] = z;
#if 0
  printf("%f\n", x);
  printf("%f\n", y);
  printf("%f\n", z);
#endif
}


/* ===============================================================
 * Rotation un atome par rmat. 
 * Resultat dans out 
 * =============================================================== */
void singleRotateOut(DtPoint3 p,DtMatrix4x4 rmat,DtPoint3 out) {
  out[0] = p[0] * rmat[0][0] + p[1] * rmat[1][0] + 
           p[2] * rmat[2][0] + rmat[3][0]; 
  out[1] = p[0] * rmat[0][1] + p[1] * rmat[1][1] + 
           p[2] * rmat[2][1] + rmat[3][1]; 
  out[2] = p[0] * rmat[0][2] + p[1] * rmat[1][2] + 
           p[2] * rmat[2][2] + rmat[3][2]; 
}

#if 0

/* ===============================================================
 * Alloue l'espace memoire pour une matrice de XxY doubles 
 * =============================================================== */
double ** mallocDoubleTabXxY(int x, int y) {
	double **tab;
    tab = (double**)malloc(x*sizeof(double*));
    while( x-- > 0 )
		tab[x] = (double*)malloc(y*sizeof(double));
	return tab;
}

void printConfs(DtStructConf *firstConf, int nbAtoms) {
	int i;
	if( firstConf != NULL ) {
		printf("#NEW_CONF\n%f\n", firstConf->energie);
                for( i = 0; i < nbAtoms; i++)
			printf("%f %f %f\n", firstConf->coords[i][0], firstConf->coords[i][1], 
			       firstConf->coords[i][2]);
		printConfs(firstConf->nextConf, nbAtoms);
	}
}

#endif


/* ==================================================================
 * Get dihedral values for one conformation
 * Install in the actual angular values in dValues.
 * Also identify the closest canonical value.
 * 
 * initp: if not 0, we just assign the dValues
 *        else, we also compare with (supposed) value
 *
 * NOTE: dValues are DEGREES.
 * ==================================================================
 */
int confInitAngles(int nbDiedres,  diedre *diedres, DtStructConf *aConf, int initp, int verbose)
{
  DtPoint3 *coords;
  double orival;
  double cv, md, cd;
  int i;
  int asDegrees = 1;
  int rs = 0;
  int ac, cf; 
  

  coords = aConf->coords;

  /* values is actual value                    */
  /* dValues is value propagated incrementally */
  for( i = 0; i < nbDiedres; i++ ) {
    orival = aConf->dValues[i]; 
    aConf->values[i] = dihedral( coords[diedres[i].extremites[0]], 
				  coords[diedres[i].extremites[1]], 
				  coords[diedres[i].extremites[2]], 
				  coords[diedres[i].extremites[3]], asDegrees);
    if (initp) {
      aConf->dValues[i] = aConf->values[i];
    } else {
      aConf->dValues[i] = aConf->values[i];
      /* Here we can check adequation */
      if ((aConf->values[i] - orival > 0.1) && (fabs(aConf->values[i] + orival) - 360. > 0.1)) {
	if (verbose) {
	  fprintf(stderr,"confInitAngles: diedre %d (%d %d %d %d): %lf (%lf)\n", i, 
		  diedres[i].extremites[0], 
		  diedres[i].extremites[1], 
		  diedres[i].extremites[2], 
		  diedres[i].extremites[3],
		  aConf->values[i], orival);
	}
	rs = 1;
      }
    }
    // fprintf(stderr,"confInitAngles: Attempting for canonical value\n");
    /* identify the closest canoncical value from actual value */
    cf = 0;
    md = 1000.;
    for (ac = 0; ac < diedres[i].nbVals; ac++) {
      /* canonical angular value */
      cv = diedres[i].valeurs[ac];
      // fprintf(stderr,"%d: %lf / %lf\n", ac, cv, aConf->dValues[i]);
      /* angular difference to canonical value */
      cd = fabs(aConf->dValues[i] - cv);
      if (cd > 180.) cd = 360. - cd;
      
      /* If best ? */
      if (cd < md) {
	md = cd;
	cf = ac;
      }
    }
    // fprintf(stderr,"found %d\n", cf);
    aConf->cvals[i] = cf;
    if (verbose) {
      fprintf(stderr,"confInitAngles: diedre %d #cvals %d found %d for %lf (%lf)\n", i, 
	      diedres[i].nbVals, aConf->cvals[i], aConf->dValues[i], diedres[i].valeurs[aConf->cvals[i]]);
    }
  } /* end nbDiedres */
  return rs;
}

/* ==================================================================
 * Get dihedral values for one set of coordinates.
 * Put the values into the diedres record.
 * ==================================================================
 */
int initAngles(int nbDiedres,  diedre *diedres, DtPoint3 *coords, int verbose)
{
  double orival;
  int i;
  int asDegrees = 1;
  int rs = 0;

  for( i = 0; i < nbDiedres; i++ ) {
    orival = diedres[i].curval; 
    diedres[i].value = dihedral(coords[diedres[i].extremites[0]], 
				coords[diedres[i].extremites[1]], 
				coords[diedres[i].extremites[2]], 
				coords[diedres[i].extremites[3]], asDegrees);
    diedres[i].curval = diedres[i].value;
    if ((diedres[i].value - orival > 0.1) && (fabs(diedres[i].value + orival) - 360. > 0.1)) {
      if (verbose) {
	fprintf(stderr,"diedre %d (%d %d %d %d): %lf (%lf)\n", i, 
		diedres[i].extremites[0], 
		diedres[i].extremites[1], 
		diedres[i].extremites[2], 
		diedres[i].extremites[3],
		diedres[i].value, orival);
      }
      rs =  1;
    }
  }
  return rs;
}

void outAngles(FILE *f, int nbDiedres,  diedre *diedres)
{
  int i;
  
  fprintf(f,"outAngles:\n");
  for( i = 0; i < nbDiedres; i++ ) {
    fprintf(f,"diedre %d (%d %d %d %d): %lf (%lf)\n", i, 
	    diedres[i].extremites[0], 
	    diedres[i].extremites[1], 
	    diedres[i].extremites[2], 
	    diedres[i].extremites[3],
	    diedres[i].value);
  }
}

/* ===========================================================
 * bestConf         : conf initiale
 * lastConfAccepted : deniere conf acceptee
 * newConf          : espace de travail
 * We perform nbPasMC tries, and get as result the bestConf.
 *
 * entry:
 *  bestConf: current conformation. Must be valid crds and dvalues.
 *  lastConfAccepted, newConf: working space
 * 
 * inside: push bestConf as lastconfAccepted
 *         select dihedral
 *         draw angular perturbation according to inidegreRot
 *         perform rotation and generate newConf
 *         compute energy
 *         accept or not
 *           if accepted:
 *              push newConf as lastconfAccepted
 *              if lowest energy:
 *                 push newConf as bestConf
 *         Note: there is no pointer swap but copy.
 * return:
 *         bestConf is up2date.
 * ===========================================================
 */
void randomRotate2(DtStructConf *bestConf, diedre *diedres, 
		   DtEnergyFactors *pE, int **nonLies, 
		   int nbPasMC, int nbAtoms, int nbDiedres, 
		   double nrjTreshold, 
		   DtStructConf *lastConfAccepted, DtStructConf *newConf, 
		   int inidegreRot, int verbose) {
  
    int numDiedreToRot, i, j, jj, k, kk, accept;
    double Mrot[4][4], dist, vdwInter;
    double rotValue, drotValue;
    int aValue;
    int degreRot = inidegreRot;
    int NOK;

    if (verbose) {
      fprintf(stderr,"RandomRotate2: Initial Energy: %lf E Threshold is %lf. Nb MC steps: %d degreRot %d\n", 
	      bestConf->energie, nrjTreshold, nbPasMC, degreRot);
    }

    if ( bestConf->energie <  nrjTreshold) goto Lexit;

    /* Push coordinates in lastConf */
    for( i = 0; i < nbAtoms; i++ ) {
#if 0
      memcpy(lastConfAccepted->coords[i], bestConf->coords[i], sizeof(DtPoint3));
#else
      lastConfAccepted->coords[i][0] = bestConf->coords[i][0];
      lastConfAccepted->coords[i][1] = bestConf->coords[i][1];
      lastConfAccepted->coords[i][2] = bestConf->coords[i][2];
#endif      
    }
    lastConfAccepted->energie = bestConf->energie;
   
    /* Now we try for a move 
     *
     * Search for a dihedral to move and angular perturbation 
     */
    for( i = 0; i < nbPasMC; i++) {
      j=0;
      // Select a dihedral, but not a peptidic bond. (?? Maybe since managed by disambiguation layer ??)
      do {
	numDiedreToRot = (int)((double)nbDiedres*(rand() / (RAND_MAX + 1.0)));
	if( j++ > 100 )
	  return;
#if 0
      } while( diedres[numDiedreToRot].isPepLike );
#else
      } while( diedres[numDiedreToRot].nbVals < 3 );
#endif
      
      // Essai d'aller plus vite (rotameres)
      drotValue = (int)( ((double)2*degreRot+1.0) *(rand() / (RAND_MAX + 1.0)) - degreRot) ;
      rotValue = DTOR(drotValue);
      
      /* 
       * Perform rotation
       *
       */
      // initialisation des coordonnees des atomes qui ne vont pas tourner
      for( j = 0; j < diedres[numDiedreToRot].nbAtomsNot2move; j++ ) {
	jj = diedres[numDiedreToRot].atomsNot2move[j];
#if 0
	memccpy(newConf->coords[jj], lastConfAccepted->coords[jj], sizeof(DtPoint3));
#else
	newConf->coords[jj][0] = lastConfAccepted->coords[jj][0];
	newConf->coords[jj][1] = lastConfAccepted->coords[jj][1];
	newConf->coords[jj][2] = lastConfAccepted->coords[jj][2];
#endif
      }
      
      //creation de la matrice de rotation. on tire dans +- degreRot. RAND_MAX is big enough: 2147483647
      MkArbitraryAxisRotMat4x4(lastConfAccepted->coords[diedres[numDiedreToRot].extremites[1]], 
			       lastConfAccepted->coords[diedres[numDiedreToRot].extremites[2]],
			       // DTOR( (int)( ((double)2*degreRot+1.0) *(rand() / (RAND_MAX + 1.0)) - degreRot) ), Mrot);	
			       rotValue, Mrot);	
      
      //mise a jour des coordonnees des atomes apres rotation
      for( j = 0; j < diedres[numDiedreToRot].nbAtoms2move; j++ ) {
	jj = diedres[numDiedreToRot].atoms2move[j];
	singleRotateOut(lastConfAccepted->coords[jj], Mrot, newConf->coords[jj]);
      }
      
      /* 
       * Now we go for energy calculation
       *
       */
      newConf->energie = 0.;

      //redefinition des interactions entre atomes qui ont et qui n'ont pas tourne
      // Grossier: VdW pour les non lies.
      for( jj = 0; jj < nbAtoms-1; jj++ ) {
	for( kk = jj+1; kk < nbAtoms; kk++ ) {
	  if( !nonLies[jj][kk-jj-1] )
	    continue;
	  dist = sqrt( (newConf->coords[jj][0] - newConf->coords[kk][0])*(newConf->coords[jj][0] - newConf->coords[kk][0]) + 
		       (newConf->coords[jj][1] - newConf->coords[kk][1])*(newConf->coords[jj][1] - newConf->coords[kk][1]) + 
		       (newConf->coords[jj][2] - newConf->coords[kk][2])*(newConf->coords[jj][2] - newConf->coords[kk][2]) );
	  vdwInter = pE->epsilonInter[jj][kk-jj-1] * 
	    pow(pE->rEtoileInterUnZeroSept[jj][kk-jj-1]/(dist+pE->rEtoileInterZeroZeroSept[jj][kk-jj-1] ), 7.) *
	    ((pE->rEtoileInterPow7UnDouze[jj][kk-jj-1]/(pow(dist, 7.)+pE->rEtoileInterPow7ZeroDouze[jj][kk-jj-1]))-2.);
	  newConf->energie += vdwInter + pE->eStaticCtInter[jj][kk-jj-1]/dist;
	}
      }

      /* 
       * Do we accept MC step ?
       *
       */
      accept = 0;
      // fprintf(stderr,"newConf generated %lf\n", newConf->energie);
      if( newConf->energie < lastConfAccepted->energie )
	accept = 1;
      else if( exp( - newConf->energie / 0.59 ) > (rand() / (RAND_MAX + 1.0)) )
	accept = 1;
      
      /* 
       * If yes install new conformation as lastAccepted and possible bestConf
       *
       */
      if( accept ) {
	// fprintf(stderr,"substep %d accepted. Diedre %d curval %f (drot %f)\n", i, numDiedreToRot, diedres[numDiedreToRot].curval + drotValue, drotValue);
	diedres[numDiedreToRot].curval += drotValue;  // TO COMMENT
	newConf->dValues[numDiedreToRot] += drotValue;
	
#if 1
	NOK = initAngles(nbDiedres, diedres, newConf->coords, 0);
	NOK = confInitAngles(nbDiedres, diedres, newConf, 0, 0);
	if (NOK) {
	  fprintf(stderr,"substep %d diedre %d deltaRot %.2lf\n",i,numDiedreToRot, RTOD(rotValue)); 
	}
#endif	  
	if (verbose) {
	  fprintf(stderr,"newConf accepted %lf\n", newConf->energie);
	}
	
	for( j = 0; j < nbAtoms; j++ ) {
#if 0
	  memccpy(lastConfAccepted->coords[j], newConf->coords[j], sizeof(DtPoint3));
#else
	  lastConfAccepted->coords[j][0] = newConf->coords[j][0];
	  lastConfAccepted->coords[j][1] = newConf->coords[j][1];
	  lastConfAccepted->coords[j][2] = newConf->coords[j][2];
#endif
	}
	lastConfAccepted->energie = newConf->energie;
	for( j = 0; j < nbDiedres; j++ ) {
	  lastConfAccepted->dValues[j] = newConf->dValues[j];
	}
	if (verbose & 0) {
	  fprintf(stderr,"newConf accepted coucou0\n");
	}
	
	if(newConf->energie < bestConf->energie) {     
	  for( j = 0; j < nbAtoms; j++ ) {
#if 0
	    memccpy( bestConf->coords[j], newConf->coords[j], sizeof(DtPoint3));
#else
	    bestConf->coords[j][0] = newConf->coords[j][0];
	    bestConf->coords[j][1] = newConf->coords[j][1];
	    bestConf->coords[j][2] = newConf->coords[j][2];
#endif
	  }
	  bestConf->energie = newConf->energie;
	  
	  if (verbose & 0) {
	    fprintf(stderr,"newConf accepted coucou1\n");
	  }
	  
	  for( j = 0; j < nbDiedres; j++ ) {
	    bestConf->dValues[j] = newConf->dValues[j];
	  }
	  if (verbose & 0) {
	    fprintf(stderr,"newConf accepted coucou2\n");
	  }
	  
	  
	}
	if ( bestConf->energie <  nrjTreshold) break;
	if (verbose & 0) {
	  fprintf(stderr,"newConf installed\n");
	}
	
      } /* end accept */
    } /* end nBPasMC */

 Lexit:     
    if (verbose) {
      fprintf(stderr,"randomRotate2 ene: %lf (Wanted less than %lf)\n", bestConf->energie, nrjTreshold);
    }
    // return bestConf;      
}


/* ====================================================
 * P. Tuffery 2009
 * This parses the lines (either from file or stdin)
 * This setups the information for generation
 * ====================================================
 */
void parseLines(char    **lines, 
		int       nLines, 
		double   *nrjTreshold, // Parametres generaux
		int      *nbPasMC,
		int      *nbConfsMax,
		DtInfos **opI,
		int verbose
		) 
{
  DtPoint3 *coords;
  int    **nonLies;
  diedre  *diedres;
  int     *actuelsConfs;
  char    *pC, *pC2;
  char    *saveptr;
  char    *token;
  char     buff[BUFSIZ];
  int      i, j, k;
  int      curLine, breakp, nbValRead;
  int      flag, posInVal;
  DtInfos *pI;

  *opI = calloc(1, sizeof(DtInfos));
  pI = *opI;

  pI->nbTotConfs  = 1;
  // pI->nbConfs = atoi(lines[0]); /* to Add into iMolecule in all cases ... */
  pI->nbAtoms = atoi(lines[0]);
  *nrjTreshold = atof(lines[1]);
  *nbPasMC = atoi(lines[2]);
  
  /* Allocation memoire */
#if 0
  pI->coords   = mallocDoubleTabXxY(pI->nbAtoms, 3);
#else
  pI->coords   = calloc(pI->nbAtoms, sizeof(DtPoint3));
#endif
  //   pI->valeurs  = mallocDoubleTabXxY(pI->nbAtoms, 5); // Obsolete
  pI->vdwparms = calloc(pI->nbAtoms, sizeof(DtPoint4));
  pI->chrg     = calloc(pI->nbAtoms, sizeof(double));
  
  coords = pI->coords;
  //  valeurs = pI->valeurs;
  
  //lecture infos atomes + valeurs
  for( i = 0; i < pI->nbAtoms; i++ ) {
    sscanf(lines[3+i],"%lf%lf%lf%lf%lf%lf%lf%lf",
	   &coords[i][0],
	   &coords[i][1],
	   &coords[i][2],
	   &pI->chrg[i],
	   &pI->vdwparms[i][0],
	   &pI->vdwparms[i][1],
	   &pI->vdwparms[i][2],
	   &pI->vdwparms[i][3]
);
  }

  //lecture infos connectivite non liee
  pI->nonLies = (int**)malloc((pI->nbAtoms-1)*sizeof(int*));
  nonLies = pI->nonLies;
  for( i = 0; i < (pI->nbAtoms)-1; i++ ) {
    nonLies[i] = (int*)malloc(((pI->nbAtoms)-i-1)*sizeof(int));
    pC2 = &lines[3+(pI->nbAtoms)+i][0];
    // fprintf(stderr,"%s\n", pC2);
    pC = pC2;
    nbValRead = 0;
    breakp = 0;
    while (breakp == 0) {
      while (((*pC2) != ' ') && ((*pC2) != '\n') && ((*pC2) != '\0')) {
	++pC2;
      }
      if ((*pC2) == '\0') breakp = 1;
      *pC2 = '\000';
      sscanf(pC, "%d", &nonLies[i][nbValRead]);
      // fprintf(stderr,"%2d",nonLies[i][nbValRead]);
      nbValRead++;
      if (nbValRead == (pI->nbAtoms)-i-1) breakp = 1;

      pC2++;
      pC = pC2;
    }
    // fprintf(stderr,"\n");
    // fprintf(stderr,"Read %d values\n", nbValRead);
  }

  //lecture infos angles diedres
  pI->nbDiedres = atoi(lines[3+(pI->nbAtoms)+(pI->nbAtoms)-1]);
#if 0
  pI->diedres =  (diedre*)malloc((pI->nbDiedres)*sizeof(diedre));
  pI->actuelsConfs = (int*)malloc((pI->nbDiedres)*sizeof(int));
#else
  pI->diedres =  (diedre*)calloc((pI->nbDiedres),sizeof(diedre));
  pI->probDiedres =  (double *)calloc((pI->nbDiedres), sizeof(double));
  pI->actuelsConfs = (int*)calloc((pI->nbDiedres),sizeof(int));
#endif
  diedres = pI->diedres;
  actuelsConfs = pI->actuelsConfs;
  for( i = 0; i < pI->nbDiedres; i++ ) actuelsConfs[i] = 0;

  for (i = 0; i < (pI->nbDiedres); i++ ) {
    curLine = 3+(pI->nbAtoms)+(pI->nbAtoms)-1+1+(i*3);

    // 1e ligne: atom defining dihedral, # canonical values, <canonical values energies>.
    strcpy(buff, lines[curLine]);
    pC = &buff[0];
    token = strtok(pC, " ");
    diedres[i].extremites[0] = atoi(token);
    token = strtok(NULL, " ");
    diedres[i].extremites[1] = atoi(token);
    token = strtok(NULL, " ");
    diedres[i].extremites[2] = atoi(token);
    token = strtok(NULL, " ");
    diedres[i].extremites[3] = atoi(token);
    token = strtok(NULL, " ");
    // fprintf(stderr,"diedre %d atom1 %d atom2 %d atom3 %d atom4 %d\n",i,diedres[i].extremites[0],diedres[i].extremites[1],diedres[i].extremites[2],diedres[i].extremites[3]);
    diedres[i].nbVals  = atoi(token);
    pI->nbTotConfs = (pI->nbTotConfs > 1000000) ?  pI->nbTotConfs : pI->nbTotConfs * diedres[i].nbVals;
    diedres[i].valeurs = malloc(diedres[i].nbVals*sizeof(double));
    diedres[i].nrjs    = malloc(diedres[i].nbVals*sizeof(double));
    for (j=0; j<diedres[i].nbVals; j++) {
      token = strtok(NULL, " ");
      diedres[i].valeurs[j] = atof(token);

      token = strtok(NULL, " ");
      diedres[i].nrjs[j] = atof(token);
    }

    // 2e ligne: # atoms that undergo rotation <atom indexes >
    strcpy(buff, lines[curLine+1]);
    pC = &buff[0];
    token = strtok(pC, " ");
    diedres[i].nbAtoms2move = atoi(token);
    // fprintf(stderr,"diedre %d #atoms: %d\n",i,diedres[i].nbAtoms2move);
    diedres[i].nbAtomsNot2move = pI->nbAtoms - diedres[i].nbAtoms2move;
    diedres[i].atoms2move      = malloc(diedres[i].nbAtoms2move    * sizeof(int));
    if (verbose) fprintf(stderr,"parseLines: diedre %d atomsToMove %d nbAtomsNot2move %d\n", i, diedres[i].nbAtoms2move, diedres[i].nbAtomsNot2move);

    diedres[i].atomsNot2move   = malloc(diedres[i].nbAtomsNot2move * sizeof(int));
    for (j=0; j<diedres[i].nbAtoms2move; j++) {
      token = strtok(NULL, " ");
      diedres[i].atoms2move[j] = atoi(token);
      // fprintf(stderr,"%3d",diedres[i].atoms2move[j]);
    }
    // fprintf(stderr,"\n");

    // 3e ligne: is it peptidic like ? (or ... ?)
    diedres[i].isPepLike = atoi(lines[curLine+2]);

  }

  *nbConfsMax = atoi(lines[3+(pI->nbAtoms)+(pI->nbAtoms)-1+((pI->nbDiedres)*3)+1]);

  pI->nbRings  = atoi(lines[3+(pI->nbAtoms)+(pI->nbAtoms)-1+((pI->nbDiedres)*3)+1+1]);
  pI->rings  = calloc(pI->nbRings, sizeof(DtRing));
  for (i = 0; i < (pI->nbRings); i++ ) {
    curLine = 3+(pI->nbAtoms)+(pI->nbAtoms)-1+((pI->nbDiedres)*3)+1+1+i+1;

    // 1e ligne: atom defining dihedral, # canonical values, <canonical values energies>.
    strcpy(buff, lines[curLine]);
    pC = &buff[0];
    token = strtok(pC, " ");
    pI->rings[i].nbAtoms = atoi(token);
    pI->rings[i].atoms = calloc(pI->rings[i].nbAtoms, sizeof(int));
    for (j=0; j<pI->rings[i].nbAtoms; j++) {
      token = strtok(NULL, " ");
      pI->rings[i].atoms[j] = atof(token);
    }
  }

  
  // indexation des atomes qui ne bougent pas pour chaque diedre
  for( i = 0; i < pI->nbDiedres; i++ ) {
    posInVal = 0;
    for( j = 0; j < pI->nbAtoms; j++ ) {
      flag = 1;
      for( k = 0; k < diedres[i].nbAtoms2move; k++ ) {
	if( j == diedres[i].atoms2move[k] ) {
	  flag = 0;
	  break;
	}
      }
      if( flag )
	diedres[i].atomsNot2move[posInVal++] = j;
    }
  }

  // fprintf(stderr,"nbTotConf: %d\n", pI->nbTotConfs);
  // exit(0);
}

/* =======================================================
   Initialize probabilities for dihedrals
   depending on balance of atoms to rotate
   =======================================================
*/
void initProbDiedres(DtInfos *pI, int verbose)
{
  int aDiedre;
  double tot = 0.;

  if (verbose)
    fprintf(stderr,"initProbdiedres (%d)\n", pI->nbDiedres);
  for (aDiedre = 0; aDiedre < pI->nbDiedres; aDiedre++) {
    pI->probDiedres[aDiedre] = (float) (pI->diedres[aDiedre].nbAtoms2move * (pI->nbAtoms - pI->diedres[aDiedre].nbAtoms2move)) / (float) (pI->nbAtoms);
    /* bias to favour sp3 moves */
    if (pI->diedres[aDiedre].nbVals == 2) pI->probDiedres[aDiedre] /= 2.;
    if (verbose)
      fprintf(stderr,"diedre %d: %d %d %.3lf\n", aDiedre, pI->diedres[aDiedre].nbAtoms2move, 
	      pI->nbAtoms - pI->diedres[aDiedre].nbAtoms2move, pI->probDiedres[aDiedre]);
    tot += pI->probDiedres[aDiedre];
  }
  for (aDiedre = 0; aDiedre < pI->nbDiedres; aDiedre++) {
    pI->probDiedres[aDiedre] /= tot;
    if (aDiedre) pI->probDiedres[aDiedre] += pI->probDiedres[aDiedre-1];
  }

  if (verbose) {
    for (aDiedre = 0; aDiedre < pI->nbDiedres; aDiedre++) {
      fprintf(stderr,"diedre %d (%d - %d): %.3lf\n", aDiedre, pI->diedres[aDiedre].nbAtoms2move, pI->diedres[aDiedre].nbAtoms2move, pI->probDiedres[aDiedre]);
    }
    fprintf(stderr,"initProbdiedres: done ...\n");
  }
}

DtEnergyFactors *initEne(DtEnergyFactors *pE, int nbAtoms, int nbDiedres, DtPoint3 *coords, double *chrg, DtPoint4 *vdwparms, int **nonLies , int verbose)
{
  void *p;
  double dist;
  double vdwInter;
  int i,j, k;

  if (verbose) fprintf(stderr,"initEne: \n");

  if (pE == NULL) {
    pE = (DtEnergyFactors *) calloc(1, sizeof(DtEnergyFactors));
    if (verbose) fprintf(stderr,"initEne: Allocated pE\n");
  }

  if (verbose) fprintf(stderr,"initEne: nbAtoms %d\n", nbAtoms);

  if (verbose) fprintf(stderr,"initEne: will alloc rEtoile (%d)\n", nbAtoms*sizeof(double));
  p = malloc(sizeof( double )*(nbAtoms+1));
  if (verbose) fprintf(stderr,"initEne: did alloc rEtoile (%d)\n", nbAtoms*sizeof(double));
  pE->rEtoile = (double *) p;
  if (verbose) fprintf(stderr,"initEne: alloc step 0\n");
  pE->rEtoileInter = (double**)malloc((nbAtoms-1)*sizeof(double*));
  if (verbose) fprintf(stderr,"initEne: alloc step 0\n");
  pE->rEtoileInterUnZeroSept = (double**)malloc((nbAtoms-1)*sizeof(double*));
  pE->rEtoileInterZeroZeroSept = (double**)malloc((nbAtoms-1)*sizeof(double*));
  pE->rEtoileInterPow7 = (double**)malloc((nbAtoms-1)*sizeof(double*));
  if (verbose) fprintf(stderr,"initEne: alloc step 1\n");
  pE->rEtoileInterPow7UnDouze = (double**)malloc((nbAtoms-1)*sizeof(double*));
  pE->rEtoileInterPow7ZeroDouze = (double**)malloc((nbAtoms-1)*sizeof(double*));
  pE->epsilonInter = (double**)malloc((nbAtoms-1)*sizeof(double*));
  pE->eStaticCtInter = (double**)malloc((nbAtoms-1)*sizeof(double*));
  pE->energiesInter = (double***)malloc((nbAtoms-1)*sizeof(double**));
  if (verbose) fprintf(stderr,"initEne: nbDiedres %d\n", nbDiedres);
  for( i = 0; i < nbDiedres; i++ )
    pE->energiesInter[i] = (double**)malloc((nbAtoms-1)*sizeof(double*));

  if (verbose) fprintf(stderr,"initEne: 1\n");
  for( i = 0; i < nbAtoms-1; i++ ) {
    pE->rEtoile[i] = vdwparms[i][2]*pow(vdwparms[i][0], 0.25);
    pE->rEtoileInter[i] = (double*)malloc((nbAtoms-i-1)*sizeof(double));
    pE->rEtoileInterUnZeroSept[i] = (double*)malloc((nbAtoms-i-1)*sizeof(double));
    pE->rEtoileInterZeroZeroSept[i] = (double*)malloc((nbAtoms-i-1)*sizeof(double));
    pE->rEtoileInterPow7[i] = (double*)malloc((nbAtoms-i-1)*sizeof(double));
    pE->rEtoileInterPow7UnDouze[i] = (double*)malloc((nbAtoms-i-1)*sizeof(double));
    pE->rEtoileInterPow7ZeroDouze[i] = (double*)malloc((nbAtoms-i-1)*sizeof(double));
    pE->epsilonInter[i] = (double*)malloc((nbAtoms-i-1)*sizeof(double));
    pE->eStaticCtInter[i] = (double*)malloc((nbAtoms-i-1)*sizeof(double));
    for( j = 0; j < nbDiedres; j++ )
      pE->energiesInter[j][i] = (double*)malloc((nbAtoms-i-1)*sizeof(double));
  }
  pE->rEtoile[nbAtoms-1] = vdwparms[nbAtoms-1][2]*pow(vdwparms[nbAtoms-1][0], 0.25);
  if (verbose) fprintf(stderr,"initEne: 2\n");

  for( i = 0; i < nbAtoms-1; i++ ) {
    for( j = i+1; j < nbAtoms; j++ ) {
      if( !nonLies[i][j-i-1] )
	continue;
      pE->rEtoileInter[i][j-i-1] = 0.5*(pE->rEtoile[i]+pE->rEtoile[j])*
	                       (1.+0.2*(1.-exp(-12.*((pE->rEtoile[i]-pE->rEtoile[j]) / (pE->rEtoile[i]+pE->rEtoile[j]))
			       *((pE->rEtoile[i]-pE->rEtoile[j]) / (pE->rEtoile[i]+pE->rEtoile[j])))));
      pE->rEtoileInterUnZeroSept[i][j-i-1] = 1.07 * pE->rEtoileInter[i][j-i-1];
      pE->rEtoileInterZeroZeroSept[i][j-i-1] = 0.07 * pE->rEtoileInter[i][j-i-1];
      pE->rEtoileInterPow7[i][j-i-1] = pow(pE->rEtoileInter[i][j-i-1], 7.);
      pE->rEtoileInterPow7UnDouze[i][j-i-1] = 1.12 * pE->rEtoileInterPow7[i][j-i-1];
      pE->rEtoileInterPow7ZeroDouze[i][j-i-1] = 0.12 * pE->rEtoileInterPow7[i][j-i-1];
      pE->epsilonInter[i][j-i-1] = ((181.16*vdwparms[i][3]*vdwparms[j][3]*vdwparms[i][0]*vdwparms[j][0]) / 
				(pow(vdwparms[i][0]/vdwparms[i][1], 1./2.) + pow(vdwparms[j][0]/vdwparms[j][1], 1./2.))) * 
	                       (1./pow(pE->rEtoileInter[i][j-i-1], 6.));
      dist = sqrt( (coords[i][0]-coords[j][0])*(coords[i][0]-coords[j][0]) + 
		   (coords[i][1]-coords[j][1])*(coords[i][1]-coords[j][1]) +
		   (coords[i][2]-coords[j][2])*(coords[i][2]-coords[j][2]) );
      if( dist < 0.0000001 )
	dist = 0.0000001;
      vdwInter = pE->epsilonInter[i][j-i-1] * 
	         pow(pE->rEtoileInterUnZeroSept[i][j-i-1]/(dist+pE->rEtoileInterZeroZeroSept[i][j-i-1] ), 7.) *
	         ((pE->rEtoileInterPow7UnDouze[i][j-i-1]/(pow(dist, 7.)+pE->rEtoileInterPow7ZeroDouze[i][j-i-1]))-2.);	
      pE->eStaticCtInter[i][j-i-1] = 332.0716*chrg[i]*chrg[j];
      pE->energiesInter[0][i][j-i-1] = vdwInter + pE->eStaticCtInter[i][j-i-1]/dist;
      for( k = 1; k < nbDiedres; k++ )
	pE->energiesInter[k][i][j-i-1] = pE->energiesInter[0][i][j-i-1];
      pE->energieTot += pE->energiesInter[0][i][j-i-1];
    }	
    //printf("\n");
  }
  free(pE->rEtoile);

  if (verbose) fprintf(stderr,"initEne: done (%lf)\n", pE->energieTot);

  // fprintf(stderr, "Energy: %lf\n",  pE->energieTot);

  return pE;
}

/* ================================================================
 * Allocate an array of sze DtStructConf.
 * Init each for nCrds and nDiedres
 * initCrds: if 0 do not allocate coords.
 * return: pointer to array
 * ================================================================
 */
DtStructConf *newStructConfArray(int sze, int nCrds, int nDiedres, int initCrds)
{
  DtStructConf *p = NULL;
  int           i;

  p = calloc(sze, sizeof(DtStructConf));
  for (i=0;i<sze;i++) {
    if (initCrds)
      p[i].coords   = calloc(nCrds, sizeof(DtPoint3));
    p[i].dValues  = calloc(nDiedres, sizeof(double));
    p[i].values   = calloc(nDiedres, sizeof(double));
    p[i].cvals    = calloc(nDiedres, sizeof(int));
    p[i].nextConf = NULL;
  }
  return p;
}

void freeStructConf(DtStructConf *p, int freeCrds)
{
  int           i;
  if (freeCrds)
    free(p->coords);
  free(p->dValues);
  free(p->values);
  free(p->cvals);
  free(p);
}

void freeStructConfArray(DtStructConf *p, int sze, int freeCrds)
{
  int           i;
  for (i=0;i<sze;i++) {
    if (freeCrds)
      free(p[i].coords);
    free(p[i].dValues);
    free(p[i].values);
    free(p[i].cvals);
  }
  free(p);
}

// calcul de l'energie
double ene(DtPoint3 *coords, int nbAtoms, DtEnergyFactors *pE, int **nonLies)
{
  double eTot = 0.;
  double vdwInter;
  double dist;
  int jj, kk;
  
  for( jj = 0; jj < nbAtoms-1; jj++ ) {
    for( kk = jj+1; kk < nbAtoms; kk++ ) {
      if( !nonLies[jj][kk-jj-1] )
	continue;
      dist = sqrt( (coords[jj][0] - coords[kk][0])*(coords[jj][0] - coords[kk][0]) + 
		   (coords[jj][1] - coords[kk][1])*(coords[jj][1] - coords[kk][1]) + 
		   (coords[jj][2] - coords[kk][2])*(coords[jj][2] - coords[kk][2]) );
      if (dist < 0.01) fprintf(stderr, "dist %d %d %lf\n", jj,kk,dist);
      //      fprintf(stderr, "epsilonInter %d %d %lf\n", jj,kk,pE->epsilonInter[jj][kk-jj-1]);
      //      fprintf(stderr, "eStaticCtInter %d %d %lf\n", jj,kk,pE->eStaticCtInter[jj][kk-jj-1]);
      vdwInter = pE->epsilonInter[jj][kk-jj-1] * 
	pow(pE->rEtoileInterUnZeroSept[jj][kk-jj-1]/(dist+pE->rEtoileInterZeroZeroSept[jj][kk-jj-1] ), 7.) *
	((pE->rEtoileInterPow7UnDouze[jj][kk-jj-1]/(pow(dist, 7.)+pE->rEtoileInterPow7ZeroDouze[jj][kk-jj-1]))-2.);
      eTot += vdwInter + pE->eStaticCtInter[jj][kk-jj-1]/dist;
    }
  }
  // fprintf(stderr,"eTot: %lf\n", eTot);
  return eTot;
}


void outInfos(FILE *f, DtInfos *p, double nrjTreshold, int nbPasMC, int nbConfsMax)
{
  int i, j;

  fprintf(f, "%d\n", p->nbAtoms);
  fprintf(f, "%.2lf\n", nrjTreshold);
  fprintf(f, "%d\n", nbPasMC);
  
  for (i=0; i< p->nbAtoms; i++) {
    fprintf(stderr,"%.4lf %.4lf %.4lf %lf %lf %lf %lf %lf\n",
	    p->coords[i][0], p->coords[i][1], p->coords[i][2], 
	    p->chrg[i],
	    p->vdwparms[i][0], p->vdwparms[i][1], p->vdwparms[i][2], p->vdwparms[i][3]); 
  }
  outElim1213(f, p);
  fprintf(f,"%d\n", p->nbDiedres);
  for (i=0; i< p->nbDiedres; i++) {
    fprintf(f,"%3d %3d %3d %3d %d",
	    p->diedres[i].extremites[0], p->diedres[i].extremites[1], 
	    p->diedres[i].extremites[2], p->diedres[i].extremites[3],
	    p->diedres[i].nbVals);
    for (j=0; j< p->diedres[i].nbVals; j++) {
      fprintf(f,"%lf %lf ", p->diedres[i].valeurs[j], p->diedres[i].nrjs[j]);
    }
    fprintf(f,"\n%d ", p->diedres[i].nbAtoms2move);
    for (j=0; j< p->diedres[i].nbAtoms2move; j++) {
      fprintf(stderr,"%3d",p->diedres[i].atoms2move[j]);
    }
    fprintf(f,"\n%d\n", p->diedres[i].isPepLike);
  }
  fprintf(f,"%d\n", nbConfsMax);
  fprintf(f,"%d\n", p->nbRings);
  for (i=0; i< p->nbRings; i++) {
    fprintf(f,"%d ", p->rings[i].nbAtoms);
    for (j=0; j< p->rings[i].nbAtoms; j++) {
      fprintf(f,"%d ", p->rings[i].atoms[j]);
    }
    fprintf(f,"\n");
  }


}


/* ================================================================
 * Call AMMOS to get energy
 * ================================================================
 */
void ammp_ene(DtStructConf *aConf, 
	      DtMol2       *pM, 
	      char *crdFNamePrefix, /* mol2 file name prefix */
	      char *ammosFName,     /* Ad-hoc tmp file name for ammos input */
	      char *ammosHome,      /* AMMOS home for minimization          */
	      char *ammosCmd,       /* AMMOS command for minimization       */
	      // int   which,          /* One  of DcFULLENE, DcBONDENE, DcNONBONDENE */
	      int   verbose)
{
  FILE * f;
  int  status;
  char fName[BUFSIZ];
  char cmd[BUFSIZ];

  if (verbose)
    fprintf(stderr,"ammp_ene: ammos file name %s\n", ammosFName);
  /* Output mol2 file for ammos */
  sprintf(fName,"%s.mol2", crdFNamePrefix);
  f = fopen(fName,"w");
  setupCoords(pM, aConf->coords);
  outMol2(f,pM, "dummy", aConf->energie, 0., 1);
  fclose(f);

  /* Generate/write .ammos file */
  f =  fopen(ammosFName, "w");
  fprintf(f, "path_of_AMMOS_SmallMol= %s\n", ammosHome);
  fprintf(f, "bank= %s.mol2\n", crdFNamePrefix);
  fclose(f);

  /* Launch energy calculation */
  sprintf(cmd, "%s/%s %s 1> /dev/null 2>&1", ammosHome, ammosCmd, ammosFName);
  if (verbose)
    fprintf(stderr, "%s\n", cmd);

  status = system(cmd);
  if (status == -1) fprintf(stderr,"Could not compute ammp energy\n");

  /* Get coordinates and energy back */
  sprintf(fName,"%s_energy.mol2", crdFNamePrefix);
  lectMol2Crds(fName, aConf, pM, verbose);
}

/* ================================================================
 * Call AMMOS minimizer
 * ================================================================
 */
void minimize(DtStructConf *aConf, 
	      DtMol2       *pM, 
	      char *crdFNamePrefix, /* mol2 file name prefix */
	      char *ammosFName, /* Ad-hoc tmp file name for ammos input */
	      char *ammosHome,  /* AMMOS home for minimization          */
	      char *ammosCmd,   /* AMMOS command for minimization       */
	      int   verbose)
{
  FILE * f;
  int  status;
  char fName[BUFSIZ];
  char cmd[BUFSIZ];

  /* Output mol2 file for ammos */
  sprintf(fName,"%s.mol2", crdFNamePrefix);
  f = fopen(fName,"w");
  setupCoords(pM, aConf->coords);
  outMol2(f,pM, "dummy", aConf->energie, 0., 1);
  fclose(f);

  /* Generate/write .ammos file */
  f =  fopen(ammosFName, "w");
  fprintf(f, "path_of_AMMOS_SmallMol= %s\n", ammosHome);
  fprintf(f, "bank= %s.mol2\n", crdFNamePrefix);
  fclose(f);

  /* Launch minimization */
  sprintf(cmd, "%s/%s %s >& /dev/null", ammosHome, ammosCmd, ammosFName);
  status = system(cmd);
  if (status == -1) fprintf(stderr,"Could not minimize\n");

  /* Get coordinates and energy back */
  sprintf(fName,"%s_minimized.mol2", crdFNamePrefix);
  lectMol2Crds(fName, aConf, pM, verbose);
}

/* ================================================================
 * This will just calculate for one conformation
 * entry:
 *   bestConf: current conformation. Must be valid crds and dvalues.
 *
 * inside:
 *   The maximal number of steps is a function off the number of dihedrals.
 *   
 * ================================================================
 */
DtStructConf * monoConf( DtStructConf    *bestConf,
			 DtStructConf    *newConf,
			 DtPoint3        *newCoords, 
			 DtStructConf    *coordsConfs,
			 DtInfos         *pI, 
			 DtEnergyFactors *pE,
			 double           nrjTreshold,
			 int              nbConfsMax, 
			 int              nbPasMC,
			 DtStructConf    *aConf, 
			 DtStructConf    *aConf2,
			 int              vibrate,
			 int              verbose
			 )
{
  DtStructConf *confTmp;
  DtMatrix4x4 Mrot;
  int i, j, jj, k;
  int NOK;
  int numDiedreToRot;
  double deltaRot;

  if (verbose) {
    fprintf(stderr,"monoConf:\n");
  }
  if( bestConf->energie < (double)nrjTreshold ) {
    // fprintf(stderr,"OK for energy: will exit\n", bestConf->energie);
    return bestConf; /* 1 conformation */
  }

  /* ============================================================== *
   * Iterate to locate one conf acceptable 
   * ============================================================== */
  
  // Actually, we do at max: nbConfsMax iterations (could be optimized)
  nbConfsMax = pI->nbDiedres * 10;
  if (verbose) 
    fprintf(stderr,"Will perform supplementary %d tries (limited to %d)\n", pI->nbTotConfs, nbConfsMax);

  /* Initialize counter */
  for (i = 0; i < pI->nbDiedres; i++) {
    pI->actuelsConfs[i] = 0;
  }

  for( i = 1; i < pI->nbTotConfs; i++ ) {

    // Kind of a systematic attempt
    numDiedreToRot = pI->nbDiedres - 1; // last dihedral
    while( ++pI->actuelsConfs[numDiedreToRot] == pI->diedres[numDiedreToRot].nbVals ) // find the first dihedral not scanned so far
      pI->actuelsConfs[numDiedreToRot--] = 0;

    // What a complicated way to perform rotations !!
    for( j = 0; j < pI->diedres[numDiedreToRot].nbAtomsNot2move; j++ ) {
      jj = pI->diedres[numDiedreToRot].atomsNot2move[j];
      /* new coords is newConf coords, different from bestConf that is initialized with pI */
      newCoords[jj][0] = coordsConfs[numDiedreToRot].coords[jj][0];
      newCoords[jj][1] = coordsConfs[numDiedreToRot].coords[jj][1];
      newCoords[jj][2] = coordsConfs[numDiedreToRot].coords[jj][2];
    }

    // THIS IS UNEFFICIENT PRESENTLY. COULD BE LARGELY OPTIMIZED !!
    // Our strategy is to generate a big move, then to perform fast optimization using monteCarlo.
    // initAngles(pI->nbDiedres, pI->diedres, coordsConfs[numDiedreToRot]);
    deltaRot = pI->diedres[numDiedreToRot].valeurs[pI->actuelsConfs[numDiedreToRot]] -
      pI->diedres[numDiedreToRot].valeurs[pI->actuelsConfs[numDiedreToRot]-1];
    pI->diedres[numDiedreToRot].curval = pI->diedres[numDiedreToRot].value + deltaRot;
    // fprintf(stderr,"diedre %d deltaRot %.2lf\n", numDiedreToRot, deltaRot);
    coordsConfs[numDiedreToRot].dValues[numDiedreToRot] = coordsConfs[numDiedreToRot].dValues[numDiedreToRot] + deltaRot;

    MkArbitraryAxisRotMat4x4(coordsConfs[numDiedreToRot].coords[pI->diedres[numDiedreToRot].extremites[1]], 
			     coordsConfs[numDiedreToRot].coords[pI->diedres[numDiedreToRot].extremites[2]],
			     DTOR(deltaRot), 
			     Mrot);	 
    
    for( j = 0; j < pI->diedres[numDiedreToRot].nbAtoms2move; j++ ) {
      jj = pI->diedres[numDiedreToRot].atoms2move[j];
      singleRotateOut( coordsConfs[numDiedreToRot].coords[jj], Mrot, newCoords[jj] );
      // mise a jour des conformations intermediaires	
      if( numDiedreToRot != pI->nbDiedres - 1 ) {
	for( k = numDiedreToRot; k < pI->nbDiedres; k++) {
	  coordsConfs[k].coords[jj][0] = newCoords[jj][0];
	  coordsConfs[k].coords[jj][1] = newCoords[jj][1];
	  coordsConfs[k].coords[jj][2] = newCoords[jj][2];
	}
      }
    }
    if( numDiedreToRot != pI->nbDiedres - 1 ) {
      for( k = numDiedreToRot; k < pI->nbDiedres; k++) {
	coordsConfs[k].dValues[numDiedreToRot] = coordsConfs[numDiedreToRot].dValues[numDiedreToRot];
      }
    }

    if (verbose) {
      fprintf(stderr,"monoConf: Will confInitAngles ...\n");
    }

    NOK = initAngles(pI->nbDiedres, pI->diedres, newConf->coords, 0);
    NOK = confInitAngles(pI->nbDiedres, pI->diedres, &coordsConfs[numDiedreToRot], 0, 0);
#if 1
    if (NOK) {
      fprintf(stderr,"Angle inconsistency on entry of randomRotate2()\n");
    }
#endif
    // exit(0);
    
    // calcul de l'energie
    newConf->energie = ene(newConf->coords, pI->nbAtoms, pE, pI->nonLies);

    /* Vibrate the conformation */
    if (vibrate)
      randomRotate2(newConf, pI->diedres, pE, pI->nonLies,
		    nbPasMC, pI->nbAtoms, pI->nbDiedres, 
		    nrjTreshold, aConf, aConf2, 
		    10, verbose);

    if (verbose) {
      fprintf(stderr,"monoConf: Will confInitAngles after randomRotate2 ...\n");
    }


    NOK = confInitAngles(pI->nbDiedres, pI->diedres, newConf, 0, 0);
    if (NOK) {
      fprintf(stderr,"Angle inconsistency on return of randomRotate2()\n");
    }

#if 0
    // NOK = initAngles(pI->nbDiedres, pI->diedres, newConf->coords, 1);
    if (NOK) {
      fprintf(stderr,"Angle inconsistency on return of randomRotate2()\n");
    }
#endif

    if( newConf->energie < bestConf->energie ) {
      fprintf(stderr, "New best energy: %lf\n",  newConf->energie); 
      confTmp = bestConf;
      bestConf = newConf;
      newConf = confTmp;
      newCoords = confTmp->coords;
    }
    // nrjTreshold is absolute value. but depends on size of compound ...
    if( (bestConf->energie < nrjTreshold) || (i > nbConfsMax) )
      break;    
  }
  
  fprintf(stderr, "Final best energy: %lf\n",  bestConf->energie); 
  return bestConf;
}

/* ===============================================================
 * Insere la conformation confToInsert dans la file classee
 * par ordre d'energie.
 * nbConfs : taille maximale de la file
 * firstConf->first: pointer to the list off currently available 
 *                   conformations for final output.
 * firstConf->outted: a conformation available as further working space.
 *
 * So: we insert the new conformation in proper place in the list (rankd by increasing energy)
 *     we return if possible a conformation for further work.
 * NOTE: We only play with pointers. There is no copy here.
 * =============================================================== */
int inserer( DtStructConf *firstConf, 
	     DtStructConf *confToInsert, 
	     int nbConfs, 
	     double eTreshold,
	     DtStructFirstAndLastFilo *firstAndOutted) 
{
  DtStructConf *pConf, *nextConf, *curConf;
  double eMax;
  int lSze; 

  /* Empty list: confToInsert devient first */
  if( firstConf == NULL )  {
    firstAndOutted->first = confToInsert;
    return 1;
  }
  if (confToInsert->energie < firstConf->energie) {
    /* ConfToInsert is best energy */
    confToInsert->nextConf = firstConf;
    firstConf = confToInsert;
  } else {
    /* ConfToInsert is not best energy: ranked it */
    nextConf = firstConf->nextConf;
    pConf    = firstConf;
    while ((nextConf != NULL) && (nextConf->energie < confToInsert->energie)) {
      pConf = nextConf;
      nextConf = nextConf->nextConf;
    }
    pConf->nextConf = confToInsert;
    confToInsert->nextConf = nextConf;
  }

  /* NOW WE LOOK FOR OUTTED CONF */
  eMax = firstConf->energie + eTreshold;
#if 1
  pConf = firstConf;
  while (nbConfs && (pConf != NULL)) {
    pConf = pConf->nextConf;
    nbConfs--;
  }
  pConf = firstConf;
  lSze = 0;
  while (pConf != NULL) {
    if (pConf->energie < eMax) {
      lSze++;
    }
    pConf = pConf->nextConf;
  }

  if ((pConf != NULL) && (!nbConfs)) {
    pConf->nextConf = NULL;  /* Avoid infinite output for instance ... */
    firstAndOutted->outted = pConf;
  }
  firstAndOutted->first = firstConf;
  return lSze;
#else
  pConf = firstConf;
  lSze = 0;
  while (pConf != NULL) {
    if (pConf->energie < eMax) {
      lSze++;
    }
    if (lSze == nbConfs) {
      firstAndOutted->outted = pConf->nextConf;
      if (firstAndOutted->outted->nextConf != NULL) {
	pConf = firstAndOutted->outted->nextConf;
	firstAndOutted->outted->nextConf = NULL;
	while (pConf != NULL) {
	  curConf = pConf;
	  pConf = pConf->nextConf;
	  freeStructConf(curConf, 1);
	}
      }
      break;
    } 
    pConf = pConf->nextConf;
  }
  firstAndOutted->first = firstConf;
  return lSze;
#endif
}

/* ================================================================
 * Select one conformation in the pending list of conformations
 * ================================================================
 */
DtStructConf *randomConf(DtStructConf *firstConf, int *index)
{
  DtStructConf *pConf, *nextConf;
  int lSze = 0;
  int rndConf;

  // fprintf(stderr,"randomConf:\n");
  pConf = firstConf;
  while (pConf != NULL) {
    lSze++;
    pConf = pConf->nextConf;
  }
  do {
    rndConf = (int)( (double)lSze*(rand() / (RAND_MAX + 1.0)) );
  } while(rndConf > lSze);

  *index = rndConf;
  // fprintf(stderr,"randomConf: %d\n", rndConf);

  pConf = firstConf;
  while (rndConf > 0) {
    rndConf--;
    pConf = pConf->nextConf;
  }
  if (pConf == NULL) 
    fprintf(stderr,"Something went wrong in randomConf\n");
  return pConf;
}

/* =========================================
   Delete from the list of conformations
   those that have been flagged for removal by clusterize
   return: a list.
   =========================================
*/
DtStructConf *cleanupList(DtStructConf* cList,  int *refs)
{
  DtStructConf *curConf = cList;
  DtStructConf *parentConf;
  DtStructConf *nextConf;
  DtStructConf *oConf = cList;
  DtStructConf *flushList = NULL;
  int aMol = 0;

  curConf = cList;
  parentConf = NULL;
#if 0
  if (refs[0] == 0) {
    oConf = cList->nextConf;
    flushList = cList;
    parentConf = NULL;
  }
#endif
  aMol = 0;
  while (curConf != NULL) {
    nextConf = curConf->nextConf;
    if (refs[aMol] != aMol) {
      if (parentConf == NULL) {
	oConf = curConf->nextConf;
      } else {
	parentConf->nextConf = curConf->nextConf;
      }
      /* Free memory */
      freeStructConf(curConf, 1);
    } else {
      parentConf = curConf;
    }
    
    curConf = nextConf;
    aMol++;
  }
  return oConf;
}

/* ====================================================================
   Clusterize conformations according to RMSd diversity
   confList: a list of conformers
   threshold: RMSd threshold to aggregate conformers (Warning, it is RMSd not squaredRMSd

   return: a purged confList 
   ====================================================================
*/
DtStructConf *clusterize(DtStructConf *confList, int nAtoms, double threshold, int verbose)
{
  DtStructConf *pConf;
  DtStructConf *pConf2;
  DtPoint3 *c1, *c2;
  double **rmsdVals;
  int *flushp;
  int *refs;

  DtMatrix4x4 M;
  double drmsd;
  double ormsd;
  int nMol;
  int aMol;
  int aMol2;
  int aAtom;
  int refMol;
  int keepOn;

  nMol = 0;
  pConf = confList;
  while (pConf != NULL) {
    nMol++;
    pConf = pConf->nextConf;
  }

  /* Init RMSdARRAY */
  flushp    = calloc(nMol, sizeof(int));
  refs      = calloc(nMol, sizeof(int));
  rmsdVals = calloc(nMol, sizeof(double *));
  for (aMol = 0; aMol < nMol; aMol++) {
    rmsdVals[aMol] = calloc(nMol, sizeof(double));
  }
  
  /* Init workspace */
  c1 = calloc(nAtoms, sizeof(DtPoint3));
  c2 = calloc(nAtoms, sizeof(DtPoint3));
  
  if (verbose) {
    fprintf(stderr,"clusterize: Memory allocation done ...\n");
  }
  
  pConf = confList;
  for (aMol = 0; aMol < nMol; aMol++) {
    /* fill c1 */
    for (aAtom = 0; aAtom < nAtoms; aAtom++) {
      memcpy(c1[aAtom], pConf->coords[aAtom], sizeof(DtPoint3));
    }
    
    pConf2 = confList;
    aMol2 = 0;
    while ( (pConf2 != NULL) && (aMol2!=aMol+1) ) {
      pConf2 = pConf2->nextConf;
      aMol2++;
    }
      
    for (aMol2 = aMol+1; aMol2 < nMol; aMol2++) {

      /* fill c2 */
      for (aAtom = 0; aAtom < nAtoms; aAtom++) {
	memcpy(c2[aAtom], pConf2->coords[aAtom], sizeof(DtPoint3));
      }
      
      /* fast RMSd ZUKER RETURNS SQUARED VALUES */
#if 0     
      drmsd = zuker_superpose(c1, c2, nAtoms,  M);
      
      rmsdVals[aMol][aMol2] = drmsd;
      rmsdVals[aMol2][aMol] = rmsdVals[aMol][aMol2];
      fprintf(stderr, "%.3lf (%.3lf)", rmsdVals[aMol][aMol2], ormsd*ormsd);
#else
      /* frmsdd reports RMSd not squared RMSd */
      frmsd(c1,c2, nAtoms, &ormsd);
      rmsdVals[aMol][aMol2] = ormsd;
      rmsdVals[aMol2][aMol] = rmsdVals[aMol][aMol2];
      
#endif

      pConf2 = pConf2->nextConf;

    }
    // fprintf(stderr, "\n");
    pConf = pConf->nextConf;
  }


  /* clustering */
  for (aMol = 0; aMol < nMol; aMol++) {
    refs[aMol] = -1;
  }

  refMol = 0;
  keepOn = 0;
  do {
    refs[refMol] = refMol;
    //    fprintf(stderr,"refMol %d\n", refMol);
    for (aMol = 0; aMol < nMol; aMol++) {
      if (aMol == refMol) continue;
      if (flushp[aMol]) continue;
      if (rmsdVals[refMol][aMol] < threshold) {
	if (verbose)
	  fprintf(stderr,"%d tooClose from %d (%.3lf)\n",aMol, refMol, rmsdVals[refMol][aMol]);
	flushp[aMol] = 1;
	refs[aMol] = refMol;
      }
    }  
    keepOn = 0;
    aMol2 = -1;
    for (aMol = 0; aMol < nMol; aMol++) {
      if (refs[aMol] != -1) continue;

      /* Look for most distant conf from current reference */
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

  /* Now we cleanup the list */
  confList = cleanupList(confList, refs);

  return confList;
}

/* ================================================================
 * This will just calculate a series of conformers.
 *
 * entry:
 *   bestConf: current conformation. Must be valid crds, dvalues, cvals.
 *
 * inside:
 *   bestConf so far always used as reference for dValues.
 *   We use a buffer of nbDiedres conformations: coordsConfs
 *   On entry, all coordsConf elements should be identical to bestConf.
 *   coordsConfs propagated up to nbDiedres.
 *   At each step, we move the coordsConf of the rotating dihedral up to the nbDiedres coordsConf buffer
 *   in order to save rotation for further steps.
 *
 * exit:
 *   firstAndOutted->first pointer (to list) is returned.
 *   The list is sorted by increasing energies.
 *   Only nbOutConfs at max.
 * ================================================================
 */
DtStructConf *multiConf( DtStructConf             *bestConf,
			 DtStructFirstAndLastFilo *firstAndOutted,
			 DtStructConf             *coordsConfs,
			 DtInfos                  *pI, 
			 DtEnergyFactors          *pE,
			 double                    nrjTreshold,
			 int                       nbConfsMax,
			 int                       nbOutConfs,
			 int                       nbPasMC,
			 DtStructConf             *aConf, 
			 DtStructConf             *aConf2,
			 int                       vibrate,
			 int                       mini,
			 char                     *ammosPath,
			 char                     *tmpName,
			 DtMol2                   *pM,
			 int                       verbose
			 )
{
  DtStructConf  *confTmp;
  DtStructConf  *newConf;
  DtPoint3      *newCoords;
  DtMatrix4x4    Mrot;

  double deltaRot;
  double lowestE;

  int i, j, jj, k;
  int NOK;
  int nConfs;
  int numDiedreToRot;
  int a;

  /*
   * Setup entry bestConf into list.
   * That one should be very correct.
   *
   */

  fprintf(stderr,"multiConf: energy threshold %lf starting conf: \n", nrjTreshold);
  NOK = confInitAngles(pI->nbDiedres, pI->diedres, bestConf, 0, 0);
  firstAndOutted->outted = NULL; 

  // nConfs = inserer(firstAndOutted->first, bestConf, nbOutConfs, firstAndOutted);

  for( i = 0; i < pI->nbDiedres; i++ ) {
    for( j = 0; j < pI->nbAtoms; j++ ) {
      memcpy(coordsConfs[i].coords[j], bestConf->coords[j], sizeof(DtPoint3));
    }
    memcpy(coordsConfs[i].values, bestConf->values, sizeof(double) * pI->nbDiedres);
    memcpy(coordsConfs[i].dValues, bestConf->dValues, sizeof(double) * pI->nbDiedres);
    memcpy(coordsConfs[i].cvals, bestConf->cvals, sizeof(int) * pI->nbDiedres);
  }

  /* ============================================================== *
   * Iterate to locate one conf acceptable 
   * ============================================================== */
  
  // Actually, we do at max: nbConfsMax iterations (could be optimized)
  fprintf(stderr,"multiConf: Will perform  %d tries (limited to %d)\n", pI->nbTotConfs, nbConfsMax);
  if (nbConfsMax > pI->nbTotConfs) nbConfsMax = pI->nbTotConfs;

  for (i = 0; i < pI->nbDiedres; i++) {
    pI->actuelsConfs[i] = 0;
  }
  pI->actuelsConfs[pI->nbDiedres-1] = -1;

  lowestE = bestConf->energie;

  for( i = 0; i < nbConfsMax; i++ ) {

    if( firstAndOutted->outted == NULL ) {
      newConf           = (DtStructConf*)malloc(sizeof(DtStructConf));
      newConf->fullene  = 0.;
      newConf->bondene  = 0.;
      newConf->nonbondene  = 0.;
      newConf->energie  = 0.;
      newConf->coords   = calloc(pI->nbAtoms, sizeof(DtPoint3));
      newConf->dValues  = calloc(pI->nbDiedres, sizeof(double));
      newConf->values   = calloc(pI->nbDiedres, sizeof(double));
      newConf->cvals    = calloc(pI->nbDiedres, sizeof(int));
      newConf->nextConf = NULL;
    } else {
      newConf = firstAndOutted->outted;
    }
    newCoords = newConf->coords;
    
    // fprintf(stderr,"Attempt %d\n",i); 

    /*
     * Counter for systematic attempt
     * We start from last dihedral, then screen all possibilities back to first dihedral 
     */
    numDiedreToRot = pI->nbDiedres - 1; // last dihedral
    while(  (++pI->actuelsConfs[numDiedreToRot] == (pI->diedres[numDiedreToRot].nbVals) ) )// find the first dihedral not scanned so far
      {
	pI->actuelsConfs[numDiedreToRot] = 0;

	// We must reset the conformations to 0 as well ...
	if (numDiedreToRot && (coordsConfs[numDiedreToRot].cvals[numDiedreToRot])) {
	  deltaRot = pI->diedres[numDiedreToRot].valeurs[0] - coordsConfs[numDiedreToRot].dValues[numDiedreToRot];
	  // fprintf(stderr,"Should reset diedre %d to cval 0. Presently %d. Rotation %lf\n", numDiedreToRot, coordsConfs[numDiedreToRot].cvals[numDiedreToRot], deltaRot);
	  coordsConfs[numDiedreToRot].dValues[numDiedreToRot] = pI->diedres[numDiedreToRot].valeurs[0];
	  coordsConfs[numDiedreToRot].cvals[numDiedreToRot] = 0;

	  MkArbitraryAxisRotMat4x4(coordsConfs[numDiedreToRot].coords[pI->diedres[numDiedreToRot].extremites[1]], 
				   coordsConfs[numDiedreToRot].coords[pI->diedres[numDiedreToRot].extremites[2]],
				   DTOR(deltaRot), 
				   Mrot);	 
	  for( j = 0; j < pI->diedres[numDiedreToRot].nbAtoms2move; j++ ) {
	    jj = pI->diedres[numDiedreToRot].atoms2move[j];
	    singleRotate( coordsConfs[numDiedreToRot].coords[jj], Mrot);
	    // mise a jour des conformations intermediaires	
	    if( numDiedreToRot != pI->nbDiedres - 1 ) {
	      for( k = numDiedreToRot; k < pI->nbDiedres; k++) {
		coordsConfs[k].coords[jj][0] = coordsConfs[numDiedreToRot].coords[jj][0];
		coordsConfs[k].coords[jj][1] = coordsConfs[numDiedreToRot].coords[jj][1];
		coordsConfs[k].coords[jj][2] = coordsConfs[numDiedreToRot].coords[jj][2];
	      }
	    }
	  }
	  if( numDiedreToRot != pI->nbDiedres - 1 ) {
	    for( k = numDiedreToRot; k < pI->nbDiedres; k++) {
	      coordsConfs[k].dValues[numDiedreToRot] = coordsConfs[numDiedreToRot].dValues[numDiedreToRot];
	      coordsConfs[k].cvals[numDiedreToRot] =  coordsConfs[numDiedreToRot].cvals[numDiedreToRot];
	    }
	  }

	}

	// Back to previous wheel.
	numDiedreToRot--;
      }

#if 1
    fprintf(stderr,"numDiedreToRot %d, attempt %d over %d\ncvals: ", numDiedreToRot, i, nbConfsMax);
    for (a = 0; a < pI->nbDiedres; a++) {
      fprintf(stderr,"%d ",pI->actuelsConfs[a]);
    }
    fprintf(stderr,"\n");
#endif

    // What a complicated way to perform rotations !!
    for( j = 0; j < pI->diedres[numDiedreToRot].nbAtomsNot2move; j++ ) {
      jj = pI->diedres[numDiedreToRot].atomsNot2move[j];
      /* new coords is newConf coords, different from bestConf that is initialized with pI */
      newCoords[jj][0] = coordsConfs[numDiedreToRot].coords[jj][0];
      newCoords[jj][1] = coordsConfs[numDiedreToRot].coords[jj][1];
      newCoords[jj][2] = coordsConfs[numDiedreToRot].coords[jj][2];
    }
    
    // THIS IS UNEFFICIENT PRESENTLY. COULD BE LARGELY OPTIMIZED !!
    // Our strategy is to generate a big move, then to perform fast optimization using monteCarlo.
    // initAngles(pI->nbDiedres, pI->diedres, coordsConfs[numDiedreToRot]);
    deltaRot = pI->diedres[numDiedreToRot].valeurs[pI->actuelsConfs[numDiedreToRot]] -
      coordsConfs[numDiedreToRot].dValues[numDiedreToRot];
    // pI->diedres[numDiedreToRot].curval = pI->diedres[numDiedreToRot].value + deltaRot;
    // fprintf(stderr,"diedre %d deltaRot %.2lf\n", numDiedreToRot, deltaRot);
    coordsConfs[numDiedreToRot].dValues[numDiedreToRot] = coordsConfs[numDiedreToRot].dValues[numDiedreToRot] + deltaRot;
    coordsConfs[numDiedreToRot].cvals[numDiedreToRot] =  pI->actuelsConfs[numDiedreToRot];
    memcpy(newConf->cvals, coordsConfs[numDiedreToRot].cvals, sizeof(int)*pI->nbDiedres);
    
    MkArbitraryAxisRotMat4x4(coordsConfs[numDiedreToRot].coords[pI->diedres[numDiedreToRot].extremites[1]], 
			     coordsConfs[numDiedreToRot].coords[pI->diedres[numDiedreToRot].extremites[2]],
			     DTOR(deltaRot), 
			     Mrot);	 
    
    for( j = 0; j < pI->diedres[numDiedreToRot].nbAtoms2move; j++ ) {
      jj = pI->diedres[numDiedreToRot].atoms2move[j];
      singleRotateOut( coordsConfs[numDiedreToRot].coords[jj], Mrot, newCoords[jj] );
      // mise a jour des conformations intermediaires	
      if( numDiedreToRot != pI->nbDiedres - 1 ) {
	for( k = numDiedreToRot; k < pI->nbDiedres; k++) {
	  coordsConfs[k].coords[jj][0] = newCoords[jj][0];
	  coordsConfs[k].coords[jj][1] = newCoords[jj][1];
	  coordsConfs[k].coords[jj][2] = newCoords[jj][2];
	}
      }
    }
    if( numDiedreToRot != pI->nbDiedres - 1 ) {
      for( k = numDiedreToRot; k < pI->nbDiedres; k++) {
	coordsConfs[k].dValues[numDiedreToRot] = coordsConfs[numDiedreToRot].dValues[numDiedreToRot];
	coordsConfs[k].cvals[numDiedreToRot] =  coordsConfs[numDiedreToRot].cvals[numDiedreToRot];
      }
    }

    // NOK = initAngles(pI->nbDiedres, pI->diedres, newConf->coords, 0);
    NOK = confInitAngles(pI->nbDiedres, pI->diedres, &coordsConfs[numDiedreToRot], 0, 0);
#if 1
    if (NOK) {
      fprintf(stderr,"Angle inconsistency on entry of randomRotate2()\n");
    }
#endif
    // exit(0);
    
    // calcul de l'energie
    newConf->energie = ene(newConf->coords, pI->nbAtoms, pE, pI->nonLies);
    // fprintf(stderr,"Energy: %lf\n", newConf->energie);
    if (newConf->energie > lowestE + nrjTreshold) {
      if (vibrate) {
	randomRotate2(newConf, pI->diedres, pE, pI->nonLies,
		      nbPasMC, pI->nbAtoms, pI->nbDiedres, 
		      nrjTreshold, aConf, aConf2, 
		      10, verbose);
      }
    }

    if (mini) {
      minimize(newConf, pM, "toto", "ammos", ammosPath, DcAMMOSMINICMD, 1);
    }


    NOK = confInitAngles(pI->nbDiedres, pI->diedres, newConf, 0, 0);
    if (NOK) {
      fprintf(stderr,"Angle inconsistency on return of randomRotate2()\n");
    }
    fprintf(stderr,"newConf energy %.3lf max is %.3lf\n", newConf->energie, lowestE + nrjTreshold);
    // if (newConf->energie > lowestE + nrjTreshold) continue;

    firstAndOutted->outted = NULL;
    nConfs = inserer(firstAndOutted->first, newConf, nbOutConfs, nrjTreshold, firstAndOutted);
    fprintf(stderr,"Stored %d conformers\n",nConfs);
    lowestE = firstAndOutted->first->energie;

  }
  fprintf(stderr,"Exit loop\n");
  
  fprintf(stderr, "Final best energy: %lf\n",  firstAndOutted->first->energie); 
  return firstAndOutted->first;
}

void printCVals(int *cvals, int nbDiedres, char *msg)
{
  int k;
  fprintf(stderr,"cvals %s: ");
  
  for( k = 0; k < nbDiedres; k++) {
    fprintf(stderr,"%d ", cvals[k]);
  }
  fprintf(stderr,"\n");
}

/* ================================================================
 * Is the conformation in cvals too close from a previously sampled conformation?
 * ================================================================
 */
int tooClose(DtStructConf *firstConf, int *cvals, int nbDiedres)
{
  DtStructConf *pConf, *nextConf;
  int nTooSims = 0;
  int k, ldiff, parentSeen=0;
  int rs = 0;
  int cNum;
  char msg[BUFSIZ];

#if 0
  printCVals(cvals, nbDiedres, "tooClose check for: ");
#endif

  // fprintf(stderr,"randomConf:\n");
  pConf = firstConf;
  cNum = 0;

  while (pConf != NULL) {
#if 0
    sprintf(msg,"curConf %d : ", cNum);
    printCVals(pConf->cvals, nbDiedres, msg);
#endif
    ldiff = 0;
    for( k = 0; k < nbDiedres; k++) {
      ldiff += (pConf->cvals[k] == cvals[k]) ? 0 : 1;
    }
    if (!ldiff) nTooSims += 1;
#if 0
    // if (parentSeen && (ldiff < 2)) nTooSims += 1;
    if (parentSeen && (ldiff < 1)) nTooSims += 1;
    if ((!ldiff) && (!parentSeen)) parentSeen = 1;
#endif
    pConf = pConf->nextConf;
    cNum ++;
  }
#if 0
  fprintf(stderr,"tooClose: %d\n", nTooSims);
#endif
  return nTooSims;
}


/* ================================================================
 * This will just calculate a series of conformers.
 *
 * entry:
 *   bestConf: current conformation. Must be valid crds, dvalues, cvals.
 *
 * inside:
 *   bestConf so far always used as reference for dValues.
 *   newConf is the pointer to current conformation.
 *   coordsConf not used any longer. Kept for compatibility.
 *
 * exit:
 *   firstAndOutted->first pointer (to list) is returned.
 *   The list is sorted by increasing energies.
 *   Only nbOutConfs at max.
 * ================================================================
 */
DtStructConf *multiConfReduced( DtStructConf             *bestConf,
				DtStructFirstAndLastFilo *firstAndOutted,
				DtStructConf             *coordsConfs,
				DtInfos                  *pI, 
				DtEnergyFactors          *pE,
				double                    nrjTreshold,
				int                       nbConfsMax,
				int                       nbOutConfs,
				int                       nbPasMC,
				DtStructConf             *aConf, 
				DtStructConf             *aConf2,
				int                       vibrate,
				int                       mini,
				char                     *ammosPath,
				char                     *tmpName,
				DtMol2                   *pM,
				int                       verbose
				)
{
  DtStructConf  *curConf;
  DtStructConf  *newConf;
  DtPoint3      *newCoords;
  DtMatrix4x4    Mrot;

  double deltaRot;
  double lowestE;
  double drawProb;

  int i, j, jj, k;
  int NOK;
  int nConfs;
  int numDiedreToRot;
  int angleIndex;
  int confIndex;
  int nTries;

  /* ============================================================== *
   * Iterate to locate one conf acceptable 
   * ============================================================== */
  
  // Actually, we do at max: nbConfsMax iterations (could be optimized)
  fprintf(stderr,"multiConfReduced: Will perform  %d tries (limited to %d)\n", pI->nbTotConfs, nbConfsMax);

  if (nbConfsMax > pI->nbTotConfs) nbConfsMax = pI->nbTotConfs;

  /*
   * Setup entry bestConf into list.
   * That one should be very correct.
   *
   */
  NOK = confInitAngles(pI->nbDiedres, pI->diedres, bestConf, 0, 0);
  firstAndOutted->outted = NULL; 

#if 1
    fprintf(stderr,"\nInitial canonical values:");
    for( k = 0; k < pI->nbDiedres; k++) {
      fprintf(stderr,"%d:%d ",k, bestConf->cvals[k]);
    }
    fprintf(stderr,"\n");
#endif

  nConfs = inserer(firstAndOutted->first, bestConf, nbOutConfs, nrjTreshold, firstAndOutted);
#if 0
  fprintf(stderr,"Initial energy %lf threshold is %lf max E %lf\n",
	  firstAndOutted->first->energie, nrjTreshold, 
	  firstAndOutted->first->energie + nrjTreshold);
#endif
  i = 0;
  while( i++ < nbConfsMax ) {

    if (!((i+1)%10)) fprintf(stderr,".");
    if (!((i+1)%100)) fprintf(stderr," ");
    if (!((i+1)%1000)) fprintf(stderr,"\n");

    lowestE = firstAndOutted->first->energie;

    /* 
     * Init working conformation
     *
     */
    if( firstAndOutted->outted == NULL ) { /* We do not have a buffer conformation already */
      newConf           = (DtStructConf*)malloc(sizeof(DtStructConf));
      newConf->fullene  = 0.;
      newConf->bondene  = 0.;
      newConf->nonbondene  = 0.;
      newConf->energie  = 0.;
      newConf->coords   = calloc(pI->nbAtoms, sizeof(DtPoint3));
      newConf->dValues  = calloc(pI->nbDiedres, sizeof(double));
      newConf->values   = calloc(pI->nbDiedres, sizeof(double));
      newConf->cvals    = calloc(pI->nbDiedres, sizeof(int));
      newConf->nextConf = NULL;
    } else {
      newConf = firstAndOutted->outted;
    }
    newCoords = newConf->coords;
    
    nTries = 0;
    do {
      /* 
       * Select a starting conformation into list
       *
       */
      curConf = randomConf(firstAndOutted->first, &confIndex);
      /* 
       * Select a transformation
       *
       */
      /* Select a dihedral */
      do {
#if 0
	drawProb = rand() / (RAND_MAX + 1.0);
	if (drawProb < pI->probDiedres[0]) {
	  numDiedreToRot = 0;
	} else {
	  for (numDiedreToRot = 1; numDiedreToRot < pI->nbDiedres; numDiedreToRot++) {
	    if (drawProb < pI->probDiedres[numDiedreToRot]) break;
	  }
	}
#else
	numDiedreToRot = (int)( (double)pI->nbDiedres*(rand() / (RAND_MAX + 1.0)) );
#endif
      } while( pI->diedres[numDiedreToRot].nbVals == 1 );
      
      /* Select a canonical conformation different from present */
      if( pI->diedres[numDiedreToRot].nbVals == 2 ) {
	angleIndex = (curConf->cvals[numDiedreToRot] == 0) ? 1 : 0;
      } else {
	do {
	  angleIndex = (int)((double)pI->diedres[numDiedreToRot].nbVals*(rand() / (RAND_MAX + 1.0)));
	} while( angleIndex ==  curConf->cvals[numDiedreToRot]);
      }

      /* Update cvals */
      memcpy(newConf->cvals, curConf->cvals, sizeof(int)*pI->nbDiedres);
      newConf->cvals[numDiedreToRot] = angleIndex;
      nTries ++;
      // fprintf(stderr, "Try %d over %d\n", nTries, nbConfsMax * 2);
    } while (tooClose(firstAndOutted->first, newConf->cvals, pI->nbDiedres) && (nTries < nbConfsMax * 2));
    if (nTries >= nbConfsMax * 2) break;

#if 1
    fprintf(stderr,"Selected conf: %d:%d\n", numDiedreToRot, newConf->cvals[numDiedreToRot]);
#endif

#if 0
    fprintf(stderr,"Selected conf %p.\nInitial angle values: ",curConf);
    for( k = 0; k < pI->nbDiedres; k++) {
      fprintf(stderr,"%d:%.2f ",k, curConf->dValues[k]);
    }
#endif
#if 0
    fprintf(stderr,"\nParent conf %d canonical values:", confIndex);
    for( k = 0; k < pI->nbDiedres; k++) {
      fprintf(stderr,"%d:%d ",k, curConf->cvals[k]);
    }
    fprintf(stderr,"\n");
#endif
#if 0
    fprintf(stderr,"Current conf canonical values:");
    for( k = 0; k < pI->nbDiedres; k++) {
      fprintf(stderr,"%d:%d ",k, newConf->cvals[k]);
    }
    fprintf(stderr,"\n");
#endif

    /* Init angular values from parent conformer */
    memcpy(newConf->dValues, curConf->dValues, sizeof(double)*pI->nbDiedres);
    memcpy(newConf->values, curConf->values, sizeof(double)*pI->nbDiedres);

#if 0
    fprintf(stderr,"New conf %p.\nInitial angle values: ",newConf);
    for( k = 0; k < pI->nbDiedres; k++) {
      fprintf(stderr,"%d:%.2f ",k, newConf->dValues[k]);
    }
    fprintf(stderr,"\nInitial canonical values:");
    for( k = 0; k < pI->nbDiedres; k++) {
      fprintf(stderr,"%d:%d ",k, newConf->cvals[k]);
    }
    fprintf(stderr,"\n");
#endif

    if (verbose) {
      fprintf(stderr,"Selected diedre %d:%d peplike %d setup value %lf\n",numDiedreToRot, newConf->cvals[numDiedreToRot], pI->diedres[numDiedreToRot].isPepLike, pI->diedres[numDiedreToRot].valeurs[angleIndex]);
      // fprintf(stderr,"Will rotate diedre %d peplike %d setup value %lf\n",numDiedreToRot, pI->diedres[numDiedreToRot].isPepLike, pI->diedres[numDiedreToRot].valeurs[angleIndex]);
    }
    /* Calculate the actual delta rotation */
    deltaRot = - curConf->dValues[numDiedreToRot] + pI->diedres[numDiedreToRot].valeurs[angleIndex];
#if 0
    fprintf(stderr,"Move will be diedre %d to canonical %d delta rot %.2lf (%.2lf - %.2lf)\n",
	    numDiedreToRot, angleIndex, deltaRot,
	    pI->diedres[numDiedreToRot].valeurs[angleIndex], curConf->dValues[numDiedreToRot]);
#endif

    /* Propagate the delta rotation for memory */
    newConf->dValues[numDiedreToRot] = fmod(newConf->dValues[numDiedreToRot] + deltaRot, 360.);
    
    // What a complicated way to perform rotations !!
    for( j = 0; j < pI->diedres[numDiedreToRot].nbAtomsNot2move; j++ ) {
      jj = pI->diedres[numDiedreToRot].atomsNot2move[j];
      /* new coords is newConf coords, different from bestConf that is initialized with pI */
      newCoords[jj][0] = curConf->coords[jj][0];
      newCoords[jj][1] = curConf->coords[jj][1];
      newCoords[jj][2] = curConf->coords[jj][2];
    }
    
    // THIS IS UNEFFICIENT PRESENTLY. COULD BE LARGELY OPTIMIZED !!
    // Our strategy is to generate a big move, then to perform fast optimization using monteCarlo.
    MkArbitraryAxisRotMat4x4(curConf->coords[pI->diedres[numDiedreToRot].extremites[1]], 
			     curConf->coords[pI->diedres[numDiedreToRot].extremites[2]],
			     DTOR(deltaRot), 
			     Mrot);	 
    
    for( j = 0; j < pI->diedres[numDiedreToRot].nbAtoms2move; j++ ) {
      jj = pI->diedres[numDiedreToRot].atoms2move[j];
      singleRotateOut( curConf->coords[jj], Mrot, newCoords[jj] );
    }

    //    NOK = initAngles(pI->nbDiedres, pI->diedres, newConf->coords, 0);
    NOK = confInitAngles(pI->nbDiedres, pI->diedres, newConf, 0, 0);
#if 1
    if (NOK) {
      fprintf(stderr,"Angle inconsistency on entry of randomRotate2()\n");
    }
#endif
    // exit(0);

    /* 
     * calcul de l'energie
     *
     */
    newConf->energie = ene(newConf->coords, pI->nbAtoms, pE, pI->nonLies);
    if (verbose)
      fprintf(stderr,"Energy: %lf (acceptable %d / %lf)\n", newConf->energie, (newConf->energie < lowestE + nrjTreshold), lowestE + nrjTreshold);

    /* 
     * Small refinement
     *
     */
    if ((!pI->diedres[numDiedreToRot].isPepLike) && (newConf->energie > lowestE + nrjTreshold)) {
      if (vibrate)
	randomRotate2(newConf, pI->diedres, pE, pI->nonLies,
		      nbPasMC, pI->nbAtoms, pI->nbDiedres, 
		      nrjTreshold, aConf, aConf2, 
		      10, verbose);
    }

    if (mini) {
      minimize(newConf, pM, "toto", "ammos", ammosPath, DcAMMOSMINICMD, 1);
#if 0
      FILE *oFd;
      setupCoords(pM, newConf->coords);
      oFd = fopen("out.mol2","w");
      outMol2(oFd, pM, "dummy", newConf->energie, 0., 1);
      fclose(oFd);
      exit(0);
#endif
    }
    

    // fprintf(stderr,"energy %.2lf / %lf (%lf)\n", newConf->energie, lowestE + nrjTreshold, nrjTreshold);
    if (newConf->energie > lowestE + nrjTreshold) continue;

#if 0
    fprintf(stderr,"Generated conf %p.\nInitial angle values: ",newConf);
    for( k = 0; k < pI->nbDiedres; k++) {
      fprintf(stderr,"%d:%.2f ",k, newConf->dValues[k]);
    }
    fprintf(stderr,"\nInitial canonical values:");
    for( k = 0; k < pI->nbDiedres; k++) {
      fprintf(stderr,"%d:%d ",k, newConf->cvals[k]);
    }
    fprintf(stderr,"\n\n");
#endif

    // fprintf(stderr,"Did randomRotate2\n");
    NOK = confInitAngles(pI->nbDiedres, pI->diedres, newConf, 0, 0);
    if (NOK) {
      fprintf(stderr,"Angle inconsistency on return of randomRotate2()\n");
    }
#if 0
    fprintf(stderr,"Checked conf %p.\nInitial angle values: ",newConf);
    for( k = 0; k < pI->nbDiedres; k++) {
      fprintf(stderr,"%d:%.2f ",k, newConf->dValues[k]);
    }
    fprintf(stderr,"\nInitial canonical values:");
    for( k = 0; k < pI->nbDiedres; k++) {
      fprintf(stderr,"%d:%d ",k, newConf->cvals[k]);
    }
    fprintf(stderr,"\n\n");
#endif

    /* 
     * Store into list of solutions if OK.
     *
     */
    firstAndOutted->outted = NULL; 
    nConfs = inserer(firstAndOutted->first, newConf, nbOutConfs, nrjTreshold, firstAndOutted);
    fprintf(stderr, "Accepted. Got %d confs over %d\n", nConfs, nbOutConfs);
    if (nConfs == nbOutConfs) break;
  }
  
  fprintf(stderr, "Final best energy: %lf\n",  firstAndOutted->first->energie); 
  return firstAndOutted->first;
}



/* ===============================================================
 * ---------- Main ----------------------------------------------- 
 * =============================================================== */
int main(int argc, char *argv[]) {

  FILE *oFd;
  int i, j, jj, k, kk, l, c, flag, posInVal=0, nbValRead, extremitesRead, numDiedreToRot, nbPasMC;
  char ligne[1025], val[1025];
  char buff[BUFSIZ], at0n[6],at1n[6],at2n[6],at3n[6];
  char confName[BUFSIZ];
  char *tmpName;
  // double ***coordsConfs;      // coordonnees et energies
  // DtPoint3 **coordsConfs;
  double nrjTreshold;
  double nrjInitTreshold = 300.;
  DtStructConf *coordsConfs;

  DtPoint3 *newCoords;

  double dist, vdwInter, Mrot[4][4];
  DtStructConf *newConf, *confTmp, *bestConf;

  int verbose = 0;

  double deltaRot; 
	
  DtStructConf *randomRotateSCTmp1, *randomRotateSCTmp2; //pointeurs alloues 
  DtStructConf *aConf, *aConf2;

  DtStructFirstAndLastFilo *firstAndOutted;
	
  char   **lines;
  int      nLines;

  int angleIndex;
  double radiansToRot;
  int nbConfsMax;             // nb de conformations max a explorer
  int nbOutConfs = 100;
  int nbGenConfs = 100;
  int nbConfs;                // nb de conformations a renvoyer
  int confNum;

  double frogOriEne, ammpOriEne;
  DtEnergyFactors *pE = NULL;  // energy calculation
                
  DtInfos *pI;

  DtMolArray *pMa;

  int NOK;
  int nOut;
  int vibrate = 0;

  srand(time(NULL));

  parseargstr(argc, argv);
  argstrPrmtrsSummary(argc, argv);
  verbose = gVerbose;
  vibrate = gVibrate;

  nrjInitTreshold = gNrjInitTreshold;

  tmpName = tempnam("./","FrogTmp");

  if (IMol2FichName[0] != '\0') {
    pMa = lectMol2(IMol2FichName,  1, verbose); // Only consider the first conformer so far
  }
  if (verbose) {
    fprintf(stderr,"Mol2Input done ...\n");
  }

  mapMol2Info(pMa->pM, NULL);

#if 0  // Uncomment to debug
  outMol2(stderr,pMa->pM, 0., 1);
  exit(0);
#endif

  if (IDataFichName[0] != '\0') {
    /* either command line */
    fprintf(stderr,"command line mode\n");
    lines=txtFileRead(IDataFichName,&nLines);
  } else {
    fprintf(stderr,"pipe in mode\n");
    lines=stdinRead(&nLines);
  }
  if (DEBUG) {
    fprintf(stderr,"input data :\n");
    for (i=0; i<nLines; i++) {
      fprintf(stderr,"%s\n",lines[i]);
    }
  }

  parseLines(lines, nLines, &nrjTreshold, &nbPasMC, &nbConfsMax, &pI, verbose);
  if (verbose) {
    fprintf(stderr,"data input (%s) done ...\n", IDataFichName);
  }

  nbOutConfs = gConfMax;
  if (gConfMax == 0) nbOutConfs = nbConfsMax;
  nbGenConfs = nbOutConfs;
  if (gClusterize) nbGenConfs *= 1.5;
  fprintf(stderr,"Want %d will generate %d gCLusterize %d\n", nbOutConfs, nbGenConfs, gClusterize);

  /* Energy threshold = 100 */
  /*   nrjTreshold  = 300; */
  // outInfos(stderr, pI, nrjTreshold, nbPasMC, nbConfsMax);
  
  pI = info2Info(pMa->pM, pI);

  initProbDiedres(pI, verbose);

  // uncomment to debug
  // outInfos(stderr, pI, nrjTreshold, nbPasMC, nbConfsMax);
  // Extra info for debug
#if 0
  for (i=0; i< pI->nbDiedres; i++) {
    sscanf(pMa->pM->atoms[pI->diedres[i].extremites[0]],"%d%s",buff, at0n);
    sscanf(pMa->pM->atoms[pI->diedres[i].extremites[1]],"%d%s",buff, at1n);
    sscanf(pMa->pM->atoms[pI->diedres[i].extremites[2]],"%d%s",buff, at2n);
    sscanf(pMa->pM->atoms[pI->diedres[i].extremites[3]],"%d%s",buff, at3n);
    fprintf(stderr,"%3d %s %3d %s %3d %s %3d %s %d\n",
	    pI->diedres[i].extremites[0], at0n, 
	    pI->diedres[i].extremites[1], at1n, 
	    pI->diedres[i].extremites[2], at2n, 
	    pI->diedres[i].extremites[3], at3n, 
	    pI->diedres[i].nbVals);
  }
#endif


  // outElim1213(stderr, pI); // To check the elim1213 is OK
  // exit(0);
  if (verbose)
    fprintf(stderr,"nbAtoms: %d nbDiedres %d eTrhld: %d MCsteps %d nbTotConfs %d nbConfsMax %d nbOutConfs %d\n", 
	  pI->nbAtoms, pI->nbDiedres, nrjTreshold, nbPasMC, pI->nbTotConfs, nbConfsMax, nbOutConfs);

  // Initialisation energy calculations
  pE = initEne(pE, pI->nbAtoms, pI->nbDiedres, pI->coords, pI->chrg, pI->vdwparms, pI->nonLies, verbose );
  if (verbose)
    fprintf(stderr,"Energie inits done ...\n");

  // initialisation de la premiere conformation
  bestConf           = newStructConfArray(1, pI->nbAtoms, pI->nbDiedres, 0);
  bestConf->coords   = pI->coords;      // Coordonnees initiales en l'occurence
  bestConf->energie  = pE->energieTot; // Energie associee

#if 0
  bestConf = (DtStructConf*)malloc(sizeof(DtStructConf));
  bestConf->dValues  = calloc(pI->nbDiedres, sizeof(double));
  bestConf->values   = calloc(pI->nbDiedres, sizeof(double));
  bestConf->coords   = pI->coords;      // Coordonnees initiales en l'occurence
  bestConf->energie  = pE->energieTot; // Energie associee
  bestConf->nextConf = NULL;          // Debut d'un liste ??
#endif

  // Init angle calculation
  NOK = initAngles(pI->nbDiedres, pI->diedres, pI->coords, 1);
  NOK = confInitAngles(pI->nbDiedres, pI->diedres, bestConf, 1, 0);

  // initialisation des conformations intermediaires
  // coordsConfs = (DtPoint3 **)malloc(pI->nbDiedres*sizeof(DtPoint3 *));
  // coordsConfs = (DtStructConf *)malloc(pI->nbDiedres*sizeof(DtStructConf));
  coordsConfs = newStructConfArray(pI->nbDiedres, pI->nbAtoms, pI->nbDiedres, 1);

  for( i = 0; i < pI->nbDiedres; i++ ) {

#if 0
    coordsConfs[i].coords    = calloc(pI->nbAtoms, sizeof(DtPoint3));
    coordsConfs[i].dValues   = calloc(pI->nbDiedres, sizeof(double));
    coordsConfs[i].values    = calloc(pI->nbDiedres, sizeof(double));
    coordsConfs[i].nextConf  = NULL;
#endif

    for( j = 0; j < pI->nbAtoms; j++ ) {
      coordsConfs[i].coords[j][0] = pI->coords[j][0];
      coordsConfs[i].coords[j][1] = pI->coords[j][1];
      coordsConfs[i].coords[j][2] = pI->coords[j][2];
    }
    for( j = 0; j < pI->nbDiedres; j++ ) {
      coordsConfs[i].values[j]  = bestConf->values[j];
      coordsConfs[i].dValues[j] = bestConf->dValues[j];
      coordsConfs[i].cvals[j]   = bestConf->cvals[j];
    }
  }

#if 0
  // Uncomment for debug
  setupCoords(pMa->pM, pI->coords);
  outMol2(stdout,pMa->pM, 1);
  // exit(0);
#endif

  //allocations pour randomRotate
  aConf  = newStructConfArray(1, pI->nbAtoms, pI->nbDiedres, 1);
  aConf2 = newStructConfArray(1, pI->nbAtoms, pI->nbDiedres, 1);

#if 1
  if (verbose)
    fprintf(stderr,"Conformation inits done ... Energy is %lf Threshold is %d. Nb MC steps: %d\n", 
	    bestConf->energie, nrjTreshold, nbPasMC);
#endif

  // Small steps of 10 degrees to see.
  fprintf(stderr,"Attempt for first conformer less than %lf\n", nrjInitTreshold);
  randomRotate2(bestConf, pI->diedres, pE, pI->nonLies, 
		nbPasMC*10, pI->nbAtoms, pI->nbDiedres, // Was nbPasMC at origin
		nrjInitTreshold, aConf, aConf2, 
		10, verbose );
  fprintf(stderr, "New best energy: %lf\n",  bestConf->energie); 

  // printConfs(bestConf, pI->nbAtoms);
  NOK = initAngles(pI->nbDiedres, pI->diedres, bestConf->coords, 1);
  NOK = confInitAngles(pI->nbDiedres, pI->diedres, bestConf, 0, 0);

#if 0
  // Uncomment for debug
  setupCoords(pMa->pM, bestConf->coords);
  outMol2(stdout,pMa->pM, 1);
  //exit(0);
#endif

  // Prepare for crds buffering, angle propagation (new DtStructConf)
  newConf  = newStructConfArray(1, pI->nbAtoms, pI->nbDiedres, 1);
#if 0
  newCoords         = calloc(pI->nbAtoms, sizeof(DtPoint3));
  newConf->coords   = newCoords;
#else
  newCoords         = newConf->coords;
#endif
#if 0
  newConf           = (DtStructConf*)malloc(sizeof(DtStructConf));
  newConf->dValues  = calloc(pI->nbDiedres, sizeof(double));
  newConf->values   = calloc(pI->nbDiedres, sizeof(double));
  newConf->nextConf = NULL;
#endif

  firstAndOutted = (DtStructFirstAndLastFilo*)malloc(sizeof(DtStructFirstAndLastFilo));
  firstAndOutted->outted = NULL;
  firstAndOutted->first = NULL;

  bestConf = monoConf( bestConf,
		       newConf,
		       newCoords, 
		       coordsConfs,
		       pI, 
		       pE,
		       nrjInitTreshold,
		       nbConfsMax, 
		       nbPasMC,
		       aConf, 
		       aConf2,
		       1, // vibrate
		       verbose);


  if (gMini) {
    minimize(bestConf, pMa->pM, "toto", "ammos", gAmmosPath, DcAMMOSMINICMD, 1);
    if (verbose)
      fprintf(stderr,"Initial energy after minimization: %lf \n", bestConf->energie);
    // exit(0);
  }

  //  nrjTreshold = bestConf->energie + 100.;
  // nrjTreshold = 100.;
  if (verbose)
    fprintf(stderr,"Reference energy: %lf. nbPasMC: %d nrjTreshold: %lf)\n", bestConf->energie, nbPasMC, nrjTreshold);

  if (gMode & DcFULLMULTICONF) {
    bestConf = multiConf( bestConf,
			  firstAndOutted,
			  coordsConfs,
			  pI, 
			  pE,
			  nrjInitTreshold,
			  nbConfsMax,
			  nbGenConfs,
			  nbPasMC,
			  aConf, 
			  aConf2,
			  vibrate,
			  gMini,
			  gAmmosPath,
			  tmpName, 
			  pMa->pM,
			  verbose
			  );

  } else if (gMode & DcREDUCEDMULTICONF) {
    bestConf = multiConfReduced( bestConf,
				 firstAndOutted,
				 coordsConfs,
				 pI, 
				 pE,
				 nrjTreshold,
				 nbConfsMax,
				 nbGenConfs,
				 nbPasMC,
				 aConf, 
				 aConf2,
				 vibrate,
				 gMini,
				 gAmmosPath,
				 tmpName, 
				 pMa->pM,
				 verbose
				 );
  }

  // frogOriEne = bestConf->energie;
  //  minimize(bestConf, pMa->pM, "toto", "ammos", gAmmosPath, DcAMMOSMINICMD, 1);
  // fprintf(stderr,"Will ammp for first conf (frog Ene %lf)\n", bestConf->energie);
  if (verbose)
    fprintf(stderr,"Generation complete\n");
  ammp_ene(bestConf, pMa->pM, "toto", "ammos", gAmmosPath, DcAMMP_ENE_CMD, 1);

  if (verbose) {
    fprintf(stderr,"AMMP ENERGY: %lf (%lf %lf) frog ene %lf\n", bestConf->fullene, bestConf->bondene, bestConf->nonbondene, bestConf->energie);
    // exit(0);
  }

  
  /* Here, we could try to minimize automatically */
#if 0
    ammpOriEne = bestConf->energie;
    minimize(bestConf, pMa->pM, "toto", "ammos", gAmmosPath, DcAMMOSMINICMD, 1);
    fprintf(stderr,"frog Ene %lf ammp %lf minimized %lf\n", frogOriEne, ammpOriEne, bestConf->energie);
  }
#endif
  
  if (gClusterize)
    clusterize(bestConf, pI->nbAtoms, gRMSdTreshold, verbose);

  // printConfs(bestConf, pI->nbAtoms);
  NOK = initAngles(pI->nbDiedres, pI->diedres, bestConf->coords, 1);
  confTmp = bestConf;
  if (OMol2FichName[0] == '\000')
    oFd = stdout;
  else {
    oFd = fopen(OMol2FichName,"a");
    if (oFd == NULL) {
      fprintf(stderr,"Could not open %s\nReverting to stdout\n",OMol2FichName);
      oFd = stdout;
    }
  }
  confNum = 0;
  // fprintf(stderr,"Will output %s\n", gId);
  while (confTmp != NULL) {
    if (confTmp->energie < bestConf->energie + nrjTreshold) {
      sprintf(confName,"%s_%d", gId, ++confNum);
      setupCoords(pMa->pM, confTmp->coords);
      outMol2(oFd,pMa->pM, confName, confTmp->energie, confTmp->fullene, 1);
      confTmp = confTmp->nextConf;
      if (confNum == nbOutConfs) break;
    } else { break; }
  }
  // fprintf(stdout,"%d %lf", confNum, bestConf->energie);
fprintf(stdout,"%d %lf %lf %lf", confNum, bestConf->fullene, bestConf->nonbondene, bestConf->fullene - bestConf->nonbondene);
  if (verbose) fprintf(stderr, "%s: %d conformations generated\n", gId, confNum);
  // fprintf(stderr,"Memory cleanup 1\n");
  for( i = 0; i < pI->nbAtoms-1; i++ ) {
    free(pE->rEtoileInter[i]);
    free(pE->rEtoileInterPow7[i]);
    free(pE->epsilonInter[i]);
    free(pE->eStaticCtInter[i]);
    //free(aConf->coords[i]);
    //free(aConf2->coords[i]);

    free(pI->nonLies[i]);

  }

  /* =============================================================
   * Memory cleanup
   * =============================================================
   */

  // fprintf(stderr,"Memory cleanup 2\n");
  if (verbose)
    fprintf(stderr,"Memory cleanup\n");
  for( j = 0; j < pI->nbDiedres; j++ ) {
    for( i = 0; i < pI->nbAtoms-2; i++ ) {
      free(pE->energiesInter[j][i]);
      // free(coordsConfs[j][i]);
    }
  }
  //  for( j = 0; j < pI->nbDiedres; j++ )
  //  free(coordsConfs[j][pI->nbAtoms-1]);
  for( i = 0; i < pI->nbDiedres; i++ ) {
    free(pE->energiesInter[i]);
  }
  
  if (verbose)
    fprintf(stderr,"Memory cleanup stage 2\n");
  // fprintf(stderr,"Memory cleanup 3\n");

  free(pE->energiesInter);
  fprintf(stderr,"Memory cleanup stage 2.1\n");
  freeStructConfArray(coordsConfs, pI->nbDiedres, 1);
  fprintf(stderr,"Memory cleanup stage 2.1.0\n");
  freeStructConfArray(bestConf, 1, 0);
  freeStructConfArray(aConf, 1, 0);
  free(pE->rEtoileInter);

  free(pE->rEtoileInterPow7);
  free(pE->epsilonInter);
  free(pE->eStaticCtInter);
      fprintf(stderr,"Memory cleanup stage 2.2\n");
free(pI->chrg);
  free(pI->vdwparms);
  free(pI->nonLies);

  if (verbose)
    fprintf(stderr,"Memory cleanup stage 3\n");
  
  free(pI->actuelsConfs);
  free(pI->diedres->valeurs);
  free(pI->diedres->nrjs);
  free(pI->diedres->atoms2move);
  fprintf(stderr,"Memory cleanup stage 3.0\n");
  if (pI->diedres->atomsNot2move != NULL)
    free(pI->diedres->atomsNot2move);
  fprintf(stderr,"Memory cleanup stage 3.1\n");
  free(pI->diedres);
  

  return 0;
}


