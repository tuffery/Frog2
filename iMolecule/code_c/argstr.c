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


/* ================================================================
 * 
 * Pierre Tuffery, November 1993
 * Version 1.5
 * ================================================================
 */
/* argstr.c ------------------ */

/* Parse les arguments de la ligne de commande
 * Suppose toutes les variables globales.
 */

#include <stdlib.h>
#include <stdio.h>  
#include <string.h>

#include "Types.h"

char IDataFichName[BUFSIZ] = "\0";
char IMol2FichName[BUFSIZ] = "\0";
char OMol2FichName[BUFSIZ] = "\0";
char OSdfFichName[BUFSIZ]  = "\0";
char gId[BUFSIZ]           = "\0";
char gAmmosPath[BUFSIZ]     = "\0";

double gNrjTreshold = 50;     /* seuil energie              */
double gNrjInitTreshold = 300;  /* seuil energie pour la premiere conf */

int    gClusterize = 0;
double gRMSdTreshold = 0.8;  /* seuil RMSd for clustering */

int gNbPasMC = 1000;       /* nombre de pas monteCarlo   */
int gNbAtoms = 0;          /* nombre d'atomes du composé */
int gMode    = DcMONOCONF; /* Mono conf is the default   */
int gConfMax = 1;          /* Maximum confs returned     */
int gVerbose = 0;          /* verbose mode               */   
int gVibrate = 0;          /* vibrate mode               */   
int gMini    = 0;          /* minimize mode              */
int DEBUG    = 0;          /* DEBUG mode                 */

void delimiter(FILE *f)
{
  fprintf(f,"----------------------------------------------------\n");
}

/* =====================================================
 * output a summary of program usage 
 * =====================================================
 */
void argstrPrmtrsUsage(int argc,char *argv[])
{
  char *pC;
  /* Get program name (this is very dirty !!)*/
  if ((pC = strrchr(argv[0],'/')) != NULL) pC++;
  else  pC = argv[0];

  delimiter(stderr);
  fprintf(stderr,"%s\n",pC);
  delimiter(stderr);
  fprintf(stderr,"Accepted command line arguments:\n");
  fprintf(stderr,"  -idata file   : the data file (as for stdin)\n");
  fprintf(stderr,"  -imol2 file   : the compound file (mol2 format)\n");
  fprintf(stderr,"  -omol2 file   : the generated compound file (mol2 format)\n");
  fprintf(stderr,"  -osdf  file   : the generated compound file (SDF format)\n");
  fprintf(stderr,"  -eTrhld float   : the energy threshold to accept conformation\n");
  fprintf(stderr,"  -eIni   float : the energy threshold to accept conformation\n");
  fprintf(stderr,"  -mcSteps int  : numbre of monteCarlo steps\n");
  fprintf(stderr,"  -mono         : reduced mutli conf search\n");
  fprintf(stderr,"  -rmulti int   : reduced mutli conf search\n");
  fprintf(stderr,"  -fmulti int   : full multi conf search (return x confs max)\n");
  fprintf(stderr,"  -mini PATH    : minimize using AMMOS; PATH is PATH to AMMOS distribution\n");
  fprintf(stderr,"  -id string    : name of compound to generate conformer identifiers\n");
  fprintf(stderr,"  -rmsd x       : rmsd thershold for clustering (mol2clusterize)\n");  
  fprintf(stderr,"  -vb           : vibrate mode (small perturbations around canonical values)\n");  
  fprintf(stderr,"  -v            : verbose\n");
  fprintf(stderr,"  -h            : display this message\n");
  fprintf(stderr,"\n");
}


/* =====================================================
 * output a summary of true parameter values
 * =====================================================
 */
void argstrPrmtrsSummary(int argc,char *argv[])
{
  char *pC;
  /* Get program name */
  if ((pC = strrchr(argv[0],'/')) != NULL) pC++;
  else  pC = argv[0];

  delimiter(stderr);
  fprintf(stderr,"%s\n",pC);
  delimiter(stderr);

  if (IDataFichName[0] != '\0')
    fprintf(stderr," data file          : %s\n", IDataFichName);

  if (IMol2FichName[0] != '\0')
    fprintf(stderr," input mol2 file    : %s\n", IMol2FichName);

  if (OMol2FichName[0] != '\0')
    fprintf(stderr," output mol2 file   : %s\n", OMol2FichName);

  if (OSdfFichName[0] != '\0')
    fprintf(stderr," output SDF file    : %s\n", OSdfFichName);

  fprintf(stderr," energy threhold    : %.1lf\n",gNrjTreshold);
  fprintf(stderr," initial energy     : %.1lf\n",gNrjInitTreshold);
  fprintf(stderr," monteCarlo steps   : %d\n",gNbPasMC);
  fprintf(stderr," generation mode    : %d\n",gMode);
  fprintf(stderr," confs returned     : %d\n",gConfMax);
  fprintf(stderr," id                 : %s\n",gId);
  fprintf(stderr," vibrate            : %d\n",gVibrate);
  fprintf(stderr," minimize           : %d\n",gMini);
  fprintf(stderr," verbose            : %d\n",gVerbose);
  delimiter(stderr);
}

/* ------------------------------------------------------------
 * Les arguments de la ligne de commande. --------------------- 
 * ------------------------------------------------------------
 */
void parseargstr(int argc,char *argv[])
{
  extern int open_file();
  extern double atof();
  extern long atol();
  
  int i,j;
  
  for (i=0;i<argc;i++) {
    if (!strcmp(argv[i],"-h")) {
      argstrPrmtrsUsage(argc, argv);
      exit(0);
    }
    if (!strcmp(argv[i],"--help")) {
      argstrPrmtrsUsage(argc, argv);
      exit(0);
    }

    if (!strcmp(argv[i],"-idata")) {    /* data file name (equivalent to data passed on stdin) */
      strcpy(IDataFichName, argv[i+1]);
      i+=1;
      continue;
    }

    if (!strcmp(argv[i],"-imol2")) {    /* data file name (equivalent to data passed on stdin) */
      strcpy(IMol2FichName, argv[i+1]);
      i+=1;
      continue;
    }

    if (!strcmp(argv[i],"-omol2")) {    /* ouptut file name (mol2 formatted) */
      strcpy(OMol2FichName, argv[i+1]);
      i+=1;
      continue;
    }

    if (!strcmp(argv[i],"-osdf")) {    /* ouptut file name (SDF formatted) */
      strcpy(OSdfFichName, argv[i+1]);
      i+=1;
      continue;
    }

    if (!strcmp(argv[i],"-id")) {    /* molecule identifier */
      strcpy(gId, argv[i+1]);
      i+=1;
      continue;
    }

    if (!strcmp(argv[i],"-eTrhld")) {   /* energy threshold */
      gNrjTreshold = atof(argv[i+1]);    /* WARNING: IT IS NO LOINGER AN INTEGER !! */
      i++;
      continue;
    }

    if (!strcmp(argv[i],"-eIni")) {   /* energy threshold for first conf */
      gNrjInitTreshold = atof(argv[i+1]);    /* WARNING: IT IS NO LOINGER AN INTEGER !! */
      i++;
      continue;
    }

    if (!strcmp(argv[i],"-mcSteps")) {    /* nombre de pas monteCarlo */
      gNbPasMC = atoi(argv[i+1]);
      i++;
      continue;
    }

    if (!strcmp(argv[i],"-mini")) {    /* minimize using ammmos */
      gMini = 1;
      strcpy(gAmmosPath, argv[i+1]);
      continue;
    }

    if (!strcmp(argv[i],"-ammos")) {    /* minimize using ammmos */
      strcpy(gAmmosPath, argv[i+1]);
      continue;
    }

    if (!strcmp(argv[i],"-mono")) {    /* mono conf generation */
      gMode = 0;
      continue;
    }

    if (!strcmp(argv[i],"-rmulti")) {    /* reduced multi conf generation */
      gMode = 1;
      gConfMax = atoi(argv[i+1]);    /* WARNING: IT IS AN INTEGER !! */
      i++;
      continue;
    }

    if (!strcmp(argv[i],"-rmsd")) {    /* clusterize threshold */
      gClusterize   = 1;
      gRMSdTreshold = atof(argv[i+1]);    /* WARNING: IT IS AN DOUBLE !! */
      i++;
      continue;
    }

    if (!strcmp(argv[i],"-fmulti")) {    /* full multi conf generation */
      gMode = 2;
      gConfMax = atoi(argv[i+1]);    /* WARNING: IT IS AN INTEGER !! */
      i++;
      continue;
    }

    if (!strcmp(argv[i],"-vb")) { /* vibrate on */
      gVibrate = 1;          
    }

    if (!strcmp(argv[i],"-v")) { /* verbose on */
      gVerbose = 1;          
    }

    if (!strcmp(argv[i],"-dbg")) { /* debug on */
      DEBUG = 1;          
    }
  }
}

