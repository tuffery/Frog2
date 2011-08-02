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

/* ========================================================================== */
/*                                                                            */
/*  LigCSRre (C) 2008 - Small compounds 3D molecular similarity screening by: */
/*                                                                            */
/* Pierre Tuffery(1), Flavien Quintus(1), J. Grynberg(1), O. Sperandio(1),    */
/* M. Petitjean(2)                                                            */
/*  (1): MTi, INSERM UMR-S 973, Université Paris Diderot - Paris 7,           */
/*     	 F75013, Paris, France.                                               */
/* 	 (http://www.mti.univ-paris-diderot.fr)                               */
/* (2): CEA/DSV/iBiTec-S/SB2SM (CNRS URA 2096),                               */
/*       F91191, Gif-sur-Yvette, France.                                      */
/*                                                                            */
/* Grants to use: See the LICENSE.TXT file that comes with the package.       */
/*                                                                            */
/* ========================================================================== */
#include <stdio.h>
#include <stdlib.h>

#include "TxtFile.h"
#include "Types.h"

/* récupère le debut de la ligne (jusqu'aux crds) */
void getMol2AtomHead(char *line, char*head)
{
  char  buff[8];
  char  x[30];
  int   i;
  char *pC;

  sscanf(line,"%d %s %s", &i, buff, x);
  sprintf(head,"%4d %-5s", i, buff);
  /* pC = strstr(line, x); */
  /* if (pC != NULL) */
  /*   strncpy(head, line, pC-line); */
}

/* récupère la fin de la ligne (apres crds) */
void getMol2AtomTail(char *line, char*tail)
{
  char buff[30];
  char anom[5];
  char snom[5];
  float f;
  int i;
  char *pC;

  sscanf(line,"%d %29s %f %f %f %5s",&i,buff,&f,&f,&f,anom);
  anom[5] = '\0';
  sprintf(snom," %s ",anom);
  pC = strstr(line+20, snom);  /* + 20 to avoid cases where atom name is atom type */
  if (pC != NULL)
    strcpy(tail, pC);
}

/* récupère les coordonnées d'un atome dans une ligne de la section ATOM */
void getMol2AtomCoords(char *line,DtFloat *x,DtFloat *y,DtFloat *z)
{
  char buff[8];
  int i;

  sscanf(line,"%d %s %lf %lf %lf", &i, buff, x, y, z);
}

/* récupère le type de l'atome dans une ligne de la section ATOM */
void getMol2AtomType(char *line,char *anom)
{
  char buff[30];
  float f;
  int i;

  sscanf(line,"%d %29s %f %f %f %5s",&i,buff,&f,&f,&f,anom);

  anom[5] = '\0';
}

/* récupère le nom de l'atome dans une ligne de la section ATOM */
void getMol2AtomName(char *line,char *anom)
{
  char buff[30];
  float f;
  int i;

  sscanf(line,"%d %s",&i,anom);
  anom[5] = '\0';
}

void getMol2AtomNum(char *line, int *anum)
{
  sscanf(line,"%d", anum);
}

int getMol2LineStatus(char *d)
{
  char Entete[25];
  int i;
  float f;
  char s[10];

  strncpy(Entete, d, 9);
  Entete[9] = 0;
  if (strcmp(Entete, "@<TRIPOS>") == 0) {
    strncpy(Entete, d+9, 24);
    Entete[24] = 0;
    if (strcmp(Entete,"ATOM") == 0) return 1;
    if (strcmp(Entete,"MOLECULE") == 0) return 10;
    if (strcmp(Entete,"COMMENT") == 0) return 20;
    if (strcmp(Entete,"BOND") == 0) return 30;
  }

  if (sscanf(d,"%d %s %f %f %f %s", &i, &s, &f, &f, &f, &s) == 6) return 2; // atom data

  if (sscanf(d," REMARK ATOM %d", &i) == 1) return 21;
  if (sscanf(d," REMARK MATCH. %c", s) == 1) return 22;
  if (sscanf(d," REMARK WEIGHT. %f", &f) == 1) return 23;

  return 0;
}

void getMol2ResChainIde(char *line,char *chn)
{
  int i, id = 0;
  char c[30];
  float f;

  sscanf(line, "%d %s %f %f %f %s %d", &i, &c, &f, &f, &f, &c, &id);
  *chn = '0' + (id % 10);
}


/* ============================================================================
 * Remplit les tableaux de données utiles à Csr après lecture des fichiers mol2.
 * Alloue la mémoire pour la structure DtNuage.
 * ============================================================================
 */
DtMolArray *lectMol2(char *fname, int verbose)
{
  DtMolArray *pMa;
  DtMol2 *pM;
  int    *onNuages;

  char   **lines;
  int     *limits;

  int  nLines;
  int  aLine;
  int  lineStart;
  int  lineStat;
  int  nMol;
  int  molOpened;
//  int  nNuages;
  int  aMol;
  int  next;
  int  aAtom;
  int  nHydrogen;

  pMa = calloc(1, sizeof(DtMolArray));

  lines=txtFileRead(fname,&nLines);

  if (!lines) {
    exit(0);
  }
  if (verbose)
    fprintf(stderr,"File %s: Read %d lines\n", fname, nLines);


  /* -- Nombre de molecules -- */
  nMol    = 0;
  limits  = NULL;

  for (aLine=0; aLine < nLines; aLine++) {
    lineStat=getMol2LineStatus(lines[aLine]);

    if (lineStat == 10) { /* @<TRIPOS>MOLECULE */
      if (nMol) {
        limits[nMol-1] = aLine;
      }
      nMol++;
      limits = realloc(limits, sizeof(int)*nMol);
    }
  }

  if (nMol) {
    limits[nMol-1] = aLine;
  }

  if (verbose)
    fprintf(stderr,"File %s: Found %d Mols\n", fname, nMol);

  /* -- Allocate structures -- */
  pMa->nMol = nMol;
  pMa->pM = calloc(nMol, sizeof(DtMol2));


  /* -- Input des nuages -- */
  molOpened = 0;
  for (aMol = 0; aMol < nMol; aMol++) {

    /* -- Taille du nuage -- */
    pM           = &(pMa->pM)[aMol];
    pM->nHead    = 0;
    pM->nTail    = 0;
    pM->nAtoms   = 0;
    pM->nBonds   = 0;

    lineStart = (aMol) ? limits[aMol-1]: 0;

    // Get information about atoms and bonds
    sscanf(lines[lineStart+2],"%d%d",&pM->nAtoms, &pM->nBonds);
    for (aLine=lineStart; aLine < limits[aMol]; aLine++) {
      lineStat = getMol2LineStatus(lines[aLine]);
      // fprintf(stderr, "%s\n", lines[aLine]);
      if (lineStat==1) { /* @<TRIPOS>ATOM */
	// fprintf(stderr,"Found HEAD TAG\n");
	pM->nHead  = aLine - lineStart + 1;
      }
      else if (lineStat==30){ /* @<TRIPOS>BOND */
	// fprintf(stderr,"Found BOND TAG\n");
	pM->nTail = limits[aMol] - aLine - pM->nBonds - 1;
      } 
    }

    fprintf(stderr,"Mol %d: head %d atoms %d bonds %d tail %d\n", aMol, pM->nHead, pM->nAtoms, pM->nBonds, pM->nTail);
    /* -- Allocation memoire -- */
    pM->atomDetails  = calloc(pM->nAtoms, sizeof(DtAtmLine));
    pM->atoms        = calloc(pM->nAtoms, sizeof(char*));
    pM->index        = calloc(pM->nAtoms, sizeof(int));
    pM->bonds        = calloc(pM->nBonds, sizeof(char*));
    pM->head         = calloc(pM->nHead,  sizeof(char*));
    pM->tail         = calloc(pM->nTail,  sizeof(char*));
    /* NCONNECT NOT USED YET */

    /* -- Lecture donnees -- */
    if (verbose)
      fprintf(stderr,"HEAD\n");
    for (aLine = 0; aLine < pM->nHead; aLine ++) {
      pM->head[aLine] = lines[aLine];
      if (verbose)
	fprintf(stderr,"%s\n", lines[aLine]);
    }

    if (verbose)
      fprintf(stderr,"ATOMS\n");
    for (aLine = pM->nHead; aLine < pM->nHead + pM->nAtoms; aLine ++) {
      pM->atoms[aLine - pM->nHead] = lines[aLine];
      if (verbose)
	fprintf(stderr,"%s\n", lines[aLine]);
      getMol2AtomHead(lines[aLine], pM->atomDetails[aLine - pM->nHead].head);
      // fprintf(stderr,"HEAD: \"%s\"\n", pM->atomDetails[aLine - pM->nHead].head);
      getMol2AtomCoords(lines[aLine], 
			&pM->atomDetails[aLine - pM->nHead].x, 
			&pM->atomDetails[aLine - pM->nHead].y, 
			&pM->atomDetails[aLine - pM->nHead].z);
      // fprintf(stderr,"XYZ: %lf %lf %lf\n", pM->atomDetails[aLine - pM->nHead].x, 
      // pM->atomDetails[aLine - pM->nHead].y, 
      // pM->atomDetails[aLine - pM->nHead].z);
      getMol2AtomTail(lines[aLine], pM->atomDetails[aLine - pM->nHead].tail);
      // fprintf(stderr,"TAIL: \"%s\"\n", pM->atomDetails[aLine - pM->nHead].tail);
      // exit(0);
    }

    if (verbose)
      fprintf(stderr,"BONDS\n");
    lineStart = pM->nHead+pM->nAtoms+1;
    for (aLine = lineStart; aLine < lineStart + pM->nBonds; aLine ++) {
      pM->bonds[aLine - lineStart] = lines[aLine];
      if (verbose)
	fprintf(stderr,"%s\n", lines[aLine]);
    }

    if (verbose)
      fprintf(stderr,"TAIL\n");
    lineStart = pM->nHead + pM->nAtoms + 1 + pM->nBonds;
    for (aLine = lineStart; aLine < limits[aMol]; aLine ++) {
      pM->tail[aLine - lineStart] = lines[aLine];
      if (verbose)
	fprintf(stderr,"%s\n", lines[aLine]);
    }

    if (verbose)
      fprintf(stderr,"END\n");

    break;
  }

  pMa = mol2mol2(pMa);

  /* -- PAS DE LIBERATION MEMOIRE DE TXTFILE CAR ON POINTE SUR LES LIGNES -- */
  //  txtFileFree(lines,nLines);
  if (limits)
    free(limits);

  if (verbose)
    fprintf(stderr,"Currently %d mols\n", nMol);

  return pMa;
}

/*
 * Output to a mol2 file 
 * 
 * useDetails: use atomDetails instead of atoms lines
 *
 */
void outMol2(FILE *f, DtMol2 *pM, int useDetails)
{
  int aLine;
  for (aLine=0; aLine< pM->nHead; aLine++) {
    fprintf(f,"%s\n",pM->head[aLine]);
  }
  for (aLine=0; aLine< pM->nAtoms; aLine++) {
    if (useDetails) {
      fprintf(f,"%s %10.4lf %10.4lf %10.4lf %s\n",
	      pM->atomDetails[aLine].head,
	      pM->atomDetails[aLine].x,
	      pM->atomDetails[aLine].y,
	      pM->atomDetails[aLine].z,
	      pM->atomDetails[aLine].tail
	      );
    } else {
      fprintf(f,"%s\n",pM->atoms[aLine]);
    }
  }
  fprintf(f,"%s\n","@<TRIPOS>BOND");
  for (aLine=0; aLine< pM->nBonds; aLine++) {
    fprintf(f,"%s\n",pM->bonds[aLine]);
  }
  for (aLine=0; aLine< pM->nTail; aLine++) {
    fprintf(f,"%s\n",pM->tail[aLine]);
  }
}

/*
 * Equivalences between a Mol2 and an information file
 * It is based on coordinates for heavy atoms,
 * on bonds for hydrogens
 */
void mapMol2Info(DtMol2 *pM, DtInfos *pI)
{
  int aAtm;
  /* TODO */
  /* BY DEFAULT MOL2 AND INFO ARE SAME ORDER */
  fprintf(stderr,"mapMol2Info\n");
  for (aAtm = 0; aAtm < pM->nAtoms; aAtm++) {
    pM->index[aAtm] = aAtm;
  }
}

/*
 * Install coordinates in atomDetails
 * 
 */
void setupCoords(DtMol2 *pM, DtPoint3 *crds)
{
  int aAtm;

  for (aAtm = 0; aAtm < pM->nAtoms; aAtm++) {
    pM->atomDetails[aAtm].x = crds[pM->index[aAtm]][0];
    pM->atomDetails[aAtm].y = crds[pM->index[aAtm]][1];
    pM->atomDetails[aAtm].z = crds[pM->index[aAtm]][2];
  }  
}


void bonds(DtMol2 *pM)
{
  int  *pnBonds;
  int **pbonds;
  int   i, j, ffrom, tto;

  /* memory allocation */
  pM->bondDetails.nBonds = calloc((pM->nAtoms+1), sizeof(int));
  pM->bondDetails.bonds  = calloc((pM->nAtoms+1), sizeof(int *));
  //  for (i=0; i<pM->nBonds; i++) {
  for (i=0; i<pM->nAtoms+1; i++) {
    pM->bondDetails.bonds[i] = calloc((pM->nAtoms+1), sizeof(int)); // Warning: number from 1 ??
  }  
  pnBonds  = pM->bondDetails.nBonds;
  pbonds  = pM->bondDetails.bonds;

  //  fprintf(stderr,"Will parse %d bond lines (expecting %d)\n",pM->nBonds, pnBonds);
  /* Parse mol2 bond lines */
  for (i=0; i<pM->nBonds; i++) {
    // fprintf(stderr,"%s\n",pM->bonds[i]);
    sscanf(pM->bonds[i],"%d%d%d",&j,&ffrom,&tto);
    ffrom --;
    tto --;
    // fprintf(stderr,"%d %d %d\n",j, ffrom, tto);
    pbonds[ffrom][pnBonds[ffrom]++] = tto;
    // fprintf(stderr,"pnBonds %d %d %p %p\n",tto, ffrom,  pbonds[tto], pnBonds);
    pbonds[tto][pnBonds[tto]++] = ffrom;
    // fprintf(stderr,"done\n");
  }

  // fprintf(stderr,"Done ...\n");

}

void setElim1213 (DtMol2 *pM, DtInfos *pI)
{
  int   i, j, k, ffrom, tto;
  int   aDiedre;
  int  *nBonds;
  int **bonds;
  char *seen;
  char *buff;

  /* memory allocation */
  nBonds = pM->bondDetails.nBonds;
  bonds  = pM->bondDetails.bonds;
#if 0
  nBonds = calloc(pM->nAtoms, sizeof(int));
  bonds  = calloc(pM->nAtoms, sizeof(int *));
  for (i=0; i<pM->nBonds; i++) {
    bonds[i] = calloc(pM->nAtoms, sizeof(int));
  }  

  fprintf(stderr,"Will parse %d bond lines\n",pM->nBonds);
  /* Parse mol2 bond lines */
  for (i=0; i<pM->nBonds; i++) {
    // fprintf(stderr,"%s\n",pM->bonds[i]);
    sscanf(pM->bonds[i],"%d%d%d",&j,&ffrom,&tto);
    ffrom --;
    tto --;
    // fprintf(stderr,"%d %d %d\n",j, ffrom, tto);
    bonds[ffrom][nBonds[ffrom]++] = tto;
    bonds[tto][nBonds[tto]++] = ffrom;
  }
#endif

#if 0
  for (i=0; i<pM->nAtoms; i++) {
    fprintf(stderr,"%d %d: ",i, nBonds[i]);
    for (j=0; j<nBonds[i]; j++) {
      fprintf(stderr,"%3d",bonds[i][j]);
    }
    fprintf(stderr,"\n");
  }
#endif
  
  fprintf(stderr,"Will setup 1213\n");

  /* 12 13 */
  for( i = 0; i < (pI->nbAtoms)-1; i++ ) {
    for( j = 0; j < (pI->nbAtoms)-i-1; j++ ) {
      pI->nonLies[i][j] = 1;
    }
  }

  for (i=0; i<pM->nAtoms; i++) {
    if (pI->atmRing[i] != -1) {
      for (j=i+1; j<pM->nAtoms; j++) {
	if (pI->atmRing[i] == pI->atmRing[j]) {
	  pI->nonLies[i][j-i-1] = 0;
	}
      }
    }

    // fprintf(stderr,"atom %d: ", i);
    for (j=0; j<nBonds[i]; j++) {
      if (bonds[i][j] > i) {
	pI->nonLies[i][bonds[i][j]-i-1] = 0; // 12
	// fprintf(stderr,"12 with %d, ", bonds[i][j]);
      }
      // fprintf(stderr,"checking 13 from %d ",bonds[i][j]); 
      for (k=0; k<nBonds[bonds[i][j]]; k++) {
	if (bonds[bonds[i][j]][k] > i) {
	  pI->nonLies[i][bonds[bonds[i][j]][k]-i-1] = 0; // 13
	  // fprintf(stderr,"13 with %d, ", bonds[bonds[i][j]][k]);
	}
      }
    }
    // fprintf(stderr,"\n");
  }
}

void outElim1213(FILE *f, DtInfos *oI)
{
  int i, j;

  for( i = 0; i < (oI->nbAtoms)-1; i++ ) {
    fprintf(f,"%d: ", i);
    for( j = 0; j < (oI->nbAtoms)-i-1; j++ ) {
      fprintf(f,"%2d",oI->nonLies[i][j]);
    }
    fprintf(f,"\n");
  }

}


/* ===============================================================
 * ------------- Atoms depending on a single bond ----------------
 * at1 and at2: the bond. Atoms searched from at2.
 * Returns 0 if cycle.
 * seen && buf are working areas of size nbAtoms. Must be allocated BEFORE.
 * seen return moving atoms (at2 not included)
 * =============================================================== */
int divide(int nbAtoms, int *nBonds, int **bonds, int at1, int at2, int *seen, int *buf)
{
  int i, jdex, k, kb, kdex;
  int rs = 0;

  memset(buf,-1,sizeof(int)*nbAtoms);
  memset(seen,0,sizeof(int)*nbAtoms);

  seen[at2] = 1;
  rs ++;
  kb = 0;
  for(i = 0; i<nBonds[at2]; i++) {
    jdex = bonds[at2][i];
    if (jdex != at1) {
      seen[jdex] = 1;
      rs ++;
      buf[kb] = jdex;
      ++kb;
    }
  }

  k = 0;
  while((jdex = buf[k]) != -1) {
    ++k;
    for(i = 0; i<nBonds[jdex]; i++) {
      kdex = bonds[jdex][i];
      if (kdex == at1) {
	return 0;
      }
      if (!seen[kdex]) {
	seen[kdex] = 1;
	rs ++;
	buf[kb] = kdex;
	++kb;
      }
    }
  }
  //  fprintf(stderr,"rs end: %d\n", rs);

  // seen[at2] = 0;
  return rs;
}

/* ===============================================================
 * ------------- Ring closing on a single bond ----------------
 * at1 and at2: the bond. Atoms searched from at2.
 * Returns ring size if cycle.
 * seen && buf are working areas of size nbAtoms. Must be allocated BEFORE.
 * seen return ring atoms (at2 not included)
 * =============================================================== */
int ringPath(int nbAtoms, int *nBonds, int **bonds, int at1, int at2, int *seen, int *buf)
{
  int i, jdex, k, kb, kdex;
  int rs = 0;

  memset(buf,-1,sizeof(int)*nbAtoms);
  memset(seen,0,sizeof(int)*nbAtoms);

  seen[at2] = 1;
  rs ++;
  kb = 0;
  for(i = 0; i<nBonds[at2]; i++) {
    jdex = bonds[at2][i];
    if (jdex != at1) {
      seen[jdex] = 1;
      rs ++;
      buf[kb] = jdex;
      ++kb;
    }
  }

  k = 0;
  while((jdex = buf[k]) != -1) {
    ++k;
    for(i = 0; i<nBonds[jdex]; i++) {
      kdex = bonds[jdex][i];
      if (kdex == at1) {
	seen[at1] = 1;
	rs ++;

	/* OK, BUT WE NEED TO REMOVE EXTRA ATOMS NOT IN RINGS ... */
	
	return rs;
      }
      if (!seen[kdex]) {
	seen[kdex] = 1;
	rs ++;
	buf[kb] = kdex;
	++kb;
      }
    }
  }

  return 0;
}

/* ====================================================
 * Find rings on the basis of each bond.
 * We return rings on a form similar to that of divide.
 * This is not correct !!! We do not have rings yet.
 * ====================================================
 */
void findRings(int nbAtoms, int *nBonds, int **bonds, int *nRings, int **rings, int *seen, int *buf)
{
  int at1, at2;
  int aBond, aCycle, cycleSze, newCycle;
  int i, isSame, aRing;
  int *ringSze;

  ringSze = calloc(nbAtoms, sizeof(int));

  *nRings = 0;
  fprintf(stderr,"findRings:\n");
  for (at1 = 0; at1 < nbAtoms; at1++) {
    for (aBond = 0; aBond < nBonds[at1]; aBond++) {
      at2 = bonds[at1][aBond];
      if (at2 < at1) continue;
      if (!divide(nbAtoms, nBonds, bonds,  at1, at2, seen, buf)) {
	rings[at1][at2] = rings[at2][at1] = 1;
      }
#if 0
      cycleSze = ringPath( nbAtoms, nBonds, bonds,  at1, at2, seen, buf);
      if (cycleSze) {
	fprintf(stderr,"Found ring from atoms %d %d\n", at1, at2);
	for (i= 0; i< nbAtoms; i++) {
	  fprintf(stderr,"%2d", seen[i]);
	}
	fprintf(stderr,"\n");
	newCycle = 1;
	for (aCycle = 0; aCycle < *nRings; aCycle++) {
	  if (ringSze[aCycle] != cycleSze) continue;
	  isSame  = 1;
	  for (i= 0; i< nbAtoms; i++) {
	    if (seen[i] != rings[aRing][i]) {
	      isSame = 0;
	      break;
	    }
	  }
	  if (isSame) {
	    newCycle = 0;
	    break;
	  }
	}
	if (newCycle) {
	  aRing = *nRings;
	  for (i= 0; i< nbAtoms; i++) {
	    rings[aRing][i] = seen[i];
	  }
	  (*nRings)++;
	}
      }
#endif 
    }
  }

  fprintf(stderr,"findRings: found %d rings\n", *nRings);
  for (at1 = 0; at1 < nbAtoms; at1++) {
    for (at2= 0; at2< nbAtoms; at2++) {
      fprintf(stderr,"%2d", rings[at1][at2]);
    }
    fprintf(stderr,"\n");
  }
  
}

#if 0
void  fixRingElim1213(DtMol2 *pM, DtInfos *oI, DtInfos *pI, int *revIndex, int *atmIndex)
{
  // oI[ atmIndex[i] ] = pI(i]
  // pI[ revIndex[j] ] = oI[j]

  int i, j;

  for( i = 0; i < (oI->nbAtoms)-1; i++ ) {
    if (revIndex[i] == -1) continue;
    for( j = 0; j < (oI->nbAtoms)-i-1; j++ ) {
      if (oI->nonLies[i][j]) {
	/* Check if was nonLie in old matrix: heavy atm correspondance + heavy atom bond for H */
	atDex = j+i+1;
	/* H case */
	/* bound to  atom pM is oI */
	nBondTo = pM->bondDetails[atDex].nBonds;
	if (nBondTo == 1) {
	  bondTo = pM->bondDetails[atDex].bonds[0];
	  /* Equivalent atom in old info: */
	  bondToEquiv  = revIndex[ bondTo ];
	}
	
      }
    }
    fprintf(f,"\n");
  }

}
#endif


/*
 * Build a new info oI from old info pI and new mol2
 * It is based on coordinates for heavy atoms,
 * on bonds for hydrogens
 */
DtInfos * info2Info(DtMol2 *pM, DtInfos *pI)
{
  DtInfos *oI;
  double   dx, dy, dz;
  int      aAtm;
  int    **nonLies;
  int      i, j, k;
  int     *atmIndex;
  int     *revIndex;
  int     *unmatched;
  int     *seen, *buff;
  int    **rings, nbRings;

  /* 1. Allocate oI */
  oI = calloc(1, sizeof(DtInfos));

  oI->nbAtoms    = pM->nAtoms;     // new number of atoms
  oI->nbDiedres  = pI->nbDiedres;  // same number of dihedrals
  oI->nbTotConfs = pI->nbTotConfs; // same combinatorial size

  /* Coordinates are from the mol2 */
  oI->coords    = calloc(oI->nbAtoms, sizeof(DtPoint3));
  oI->atmRing   = calloc(oI->nbAtoms, sizeof(int));
  for( i = 0; i < (oI->nbAtoms); i++ ) {
    oI->coords[i][0] = pM->atomDetails[i].x;
    oI->coords[i][1] = pM->atomDetails[i].y;
    oI->coords[i][2] = pM->atomDetails[i].z;
    /* no ring for that atom */
    oI->atmRing[i]   = -1; 
  }

  /* 2. Which old atoms are which new atoms : atmIndex pI -> oI */
  /*    Which new atoms are which old atoms : revIndex oI -> pI */
  // rings used for internal ring detection (not correct)
  // while oI->rings and atmRing rely on external ring detection (frowns)
  rings     = calloc(oI->nbAtoms, sizeof(int *));
  seen      = calloc(oI->nbAtoms, sizeof(int));
  buff      = calloc(oI->nbAtoms, sizeof(int));
  atmIndex  = calloc(pI->nbAtoms, sizeof(int));
  revIndex  = calloc(oI->nbAtoms, sizeof(int));
  unmatched = calloc(oI->nbAtoms, sizeof(int));
  nbRings   = 0;
  for( i = 0; i < (oI->nbAtoms); i++ ) {
    revIndex[i] = -1; // -1 for no equivalence
    unmatched[i] = 1;
    rings[i] = calloc(oI->nbAtoms, sizeof(int));
  }
  // oI[ atmIndex[i] ] = pI(i)
  // pI[ revIndex[j] ] = oI[j]
  for ( i = 0; i < (pI->nbAtoms); i++ ) {
    atmIndex[i] = -1; // -1 for no equivalence
    for( j = 0; j < (oI->nbAtoms); j++ ) {
      dx = pI->coords[i][0] - oI->coords[j][0];
      dy = pI->coords[i][1] - oI->coords[j][1];
      dz = pI->coords[i][2] - oI->coords[j][2];
      if ((dx*dx + dy*dy + dz*dz) < 0.05) {
	atmIndex[i] = j;
	revIndex[j] = i;
	unmatched[j] = 0;
	break;
      }
    }
  }

  /* Ring propagation from old info */
  oI->nbRings = pI->nbRings;
  oI->rings   = calloc(oI->nbRings, sizeof(DtRing));
  for( i = 0; i < (oI->nbRings); i++ ) {
    oI->rings[i].nbAtoms = pI->rings[i].nbAtoms;
    oI->rings[i].atoms   = calloc(oI->rings[i].nbAtoms,sizeof(int));
    k = 0;
    for ( j = 0; j < oI->rings[i].nbAtoms; j++ ) {
      // fprintf(stderr,"old ring atom %d %d -> %d\n", k, pI->rings[i].atoms[j], atmIndex[pI->rings[i].atoms[j]]);
      if (atmIndex[pI->rings[i].atoms[j]] >= 0) {
	oI->rings[i].atoms[k] = atmIndex[pI->rings[i].atoms[j]];
      }
      oI->atmRing[ oI->rings[i].atoms[k++] ] = i;
    }
    oI->rings[i].nbAtoms = k;
  }
  
#if 0
  // Uncommment to debug 
  for ( i = 0; i < (oI->nbAtoms); i++ ) {
    fprintf(stderr,"atom %d equiv old %d ring %d\n", i, revIndex[i], oI->atmRing[i]);
  }
#endif
#if 0
  /* Check for debug NOT FOR PRODUCTION since H valse */
  for ( i = 0; i < (pI->nbAtoms); i++ ) {
    if (atmIndex[i] == -1) {
      fprintf(stderr,"Could not find atom equivalence for %lf %lf %lf\n", 
	      pI->coords[i][0],
	      pI->coords[i][1],
	      pI->coords[i][2]
	      ); // -1 for no equivalence
    }
  }
  for ( i = 0; i < (oI->nbAtoms); i++ ) {
    if (revIndex[i] == -1) {
      fprintf(stderr,"Could not find reverse equivalence for %lf %lf %lf\n", 
	      oI->coords[i][0],
	      oI->coords[i][1],
	      oI->coords[i][2]
	      ); // -1 for no equivalence
    }
  }
#endif

  /* The atoms that are different (no equiv) are necessarily hydrogens */
  //  oI->valeurs  = mallocDoubleTabXxY(oI->nbAtoms, 5); // Obsolete
  oI->chrg      = calloc(oI->nbAtoms, sizeof(double));
  oI->vdwparms  = calloc(oI->nbAtoms, sizeof(DtPoint4));
  oI->nonLies   = (int**)calloc((oI->nbAtoms-1), sizeof(int*));
  nonLies = oI->nonLies;
  for( i = 0; i < (oI->nbAtoms)-1; i++ ) {
    nonLies[i]  = (int*)calloc( (oI->nbAtoms)-i-1, sizeof(int));
  }
  oI->diedres      =  (diedre*) malloc((oI->nbDiedres)*sizeof(diedre));
  oI->actuelsConfs = (int*) malloc((oI->nbDiedres)*sizeof(int));
  for (i = 0; i < (pI->nbDiedres); i++ ) {
    oI->actuelsConfs[i] = 0;

    oI->diedres[i].nbVals  = pI->diedres[i].nbVals;  // This is unaffected
    oI->diedres[i].valeurs = malloc(oI->diedres[i].nbVals*sizeof(double));
    oI->diedres[i].nrjs    = malloc(oI->diedres[i].nbVals*sizeof(double));
    for (j=0; j<oI->diedres[i].nbVals; j++) {
      oI->diedres[i].valeurs[j] = pI->diedres[i].valeurs[j];
      oI->diedres[i].nrjs[j] = pI->diedres[i].nrjs[j];
    }
  }

  /* Now we install atom based information */
#if 0
  for( i = 0; i < (oI->nbAtoms)-1; i++ ) {
#else
  for( i = 0; i < (oI->nbAtoms); i++ ) {
#endif
    if (revIndex[i] != -1) {
      oI->chrg[i] = pI->chrg[revIndex[i]];
      oI->vdwparms[i][0] = pI->vdwparms[revIndex[i]][0];
      oI->vdwparms[i][1] = pI->vdwparms[revIndex[i]][1];
      oI->vdwparms[i][2] = pI->vdwparms[revIndex[i]][2];
      oI->vdwparms[i][3] = pI->vdwparms[revIndex[i]][3];
    } else {
      oI->chrg[i] = 0.;  // We assume not charged ...
      oI->vdwparms[i][0] = 0.25; 
      oI->vdwparms[i][1] = 0.8; 
      oI->vdwparms[i][2] = 4.2; 
      oI->vdwparms[i][3] = 1.209; 
    }
  }  

  /* Now, we setup ring information */
  for( i = 0; i < (oI->nbRings); i++ ) {
#if 0
    fprintf(stderr,"Old ring: ");
    for ( j = 0; j < pI->rings[i].nbAtoms; j++ ) {
      fprintf(stderr,"%d ",pI->rings[i].atoms[j]);
    }
    fprintf(stderr,"\nNew ring: ");
    for ( j = 0; j < pI->rings[i].nbAtoms; j++ ) {
      fprintf(stderr,"%d ",atmIndex[ pI->rings[i].atoms[j] ]);
    }
    fprintf(stderr,"\n");
#endif
    for ( j = 0; j < oI->rings[i].nbAtoms; j++ ) {
      oI->rings[i].atoms[j] = atmIndex[ pI->rings[i].atoms[j] ];
    }
  }

  /* Now we setup the nonLie 1-2 1-3 matrix ( from scratch or derived ???) */
  /* From scratch is better/safer                                          */
  /* Pick up either from MacLib or from iMolecule                          */
  fprintf(stderr,"Will setElim1213 from Info2Info\n");
  bonds(pM);
  setElim1213 (pM, oI); /* This manage rings correctly, IF setup in oI. */

  //  findRings(oI->nbAtoms, pM->bondDetails.nBonds, pM->bondDetails.bonds, &nbRings, rings, seen, buff);
  // exit(0);

  /* So we check against old info to correct ring information */
  // fixRingElim1213(pM, oI, pI, revIndex, atmIndex);

  // fprintf(stderr,"Will outElim1213\n");
  // outElim1213(stderr, oI);
  // exit(0);

  /* Now we setup information about dihedrals: who is rotating who is not  */
  for (i = 0; i < (pI->nbDiedres); i++ ) {
    oI->diedres[i].extremites[0] = atmIndex[ pI->diedres[i].extremites[0] ];
    oI->diedres[i].extremites[1] = atmIndex[ pI->diedres[i].extremites[1] ];
    oI->diedres[i].extremites[2] = atmIndex[ pI->diedres[i].extremites[2] ];
    oI->diedres[i].extremites[3] = atmIndex[ pI->diedres[i].extremites[3] ];
#if 0
    fprintf(stderr,"diedre %d converted into:\n",i);
    fprintf(stderr,"%d -> %d\n",pI->diedres[i].extremites[0], oI->diedres[i].extremites[0]);
    fprintf(stderr,"%d -> %d\n",pI->diedres[i].extremites[1], oI->diedres[i].extremites[1]);
    fprintf(stderr,"%d -> %d\n",pI->diedres[i].extremites[2], oI->diedres[i].extremites[2]);
    fprintf(stderr,"%d -> %d\n",pI->diedres[i].extremites[3], oI->diedres[i].extremites[3]);
#endif
    oI->diedres[i].value     = pI->diedres[i].value;
    oI->diedres[i].curval    = pI->diedres[i].curval;
    oI->diedres[i].isPepLike = pI->diedres[i].isPepLike;
    oI->diedres[i].nbAtoms2move  = divide(oI->nbAtoms, 
					  pM->bondDetails.nBonds, 
					  pM->bondDetails.bonds,
					  oI->diedres[i].extremites[1], 
					  oI->diedres[i].extremites[2], 
					  seen, 
					  buff);
    oI->diedres[i].nbAtomsNot2move = oI->nbAtoms - oI->diedres[i].nbAtoms2move;
    oI->diedres[i].atoms2move      = malloc(oI->diedres[i].nbAtoms2move    * sizeof(int));
    oI->diedres[i].atomsNot2move   = malloc(oI->diedres[i].nbAtomsNot2move * sizeof(int));
    oI->diedres[i].nbAtoms2move    = 0;
    oI->diedres[i].nbAtomsNot2move = 0;
    for (j = 0; j<oI->nbAtoms; j++) {
      if (seen[j]) {
	oI->diedres[i].atoms2move[oI->diedres[i].nbAtoms2move ++ ] = j;
      } else {
	oI->diedres[i].atomsNot2move[oI->diedres[i].nbAtomsNot2move ++ ] = j;
      }
    }

#if 0
    fprintf(stderr, "diedre %d (%d %d): %d atoms to move\n", 
	    i, 
	    oI->diedres[i].extremites[1], 
	    oI->diedres[i].extremites[2], 
	    oI->diedres[i].nbAtoms2move);
    for (j=0; j<oI->nbAtoms; j++) {
      fprintf(stderr,"%2d",seen[j]);
    }
    fprintf(stderr,"\n");
    for (j=0; j<oI->diedres[i].nbAtoms2move; j++) {
      fprintf(stderr,"%3d",oI->diedres[i].atoms2move[j]);
    }
    fprintf(stderr,"\n");
#endif

  }


  /* We are done ... */
  free(atmIndex);
  free(revIndex);
  return oI;
}

