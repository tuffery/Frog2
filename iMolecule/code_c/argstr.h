#ifndef __argstr_h__
#define __argstr_h__

#include "Types.h"

extern char IDataFichName[BUFSIZ];
extern char IMol2FichName[BUFSIZ];
extern char OMol2FichName[BUFSIZ];
extern char OSdfFichName[BUFSIZ];
extern char gId[BUFSIZ];   /* Conformer identifier       */
extern char gAmmosPath[BUFSIZ];

extern double gNrjTreshold;   /* seuil energie              */
extern double gNrjInitTreshold; /* Energie premiere conformation */

extern int    gClusterize;
extern double gRMSdTreshold;  /* seuil RMSd for clustering */

extern int gNbPasMC;       /* nombre de pas monteCarlo   */
extern int gNbAtoms;       /* nombre d'atomes du composé */
extern int gMode;          /* Generation mode, one of 0 (mono), 1 (rmulti), 2 (fmulti) */
extern int gConfMax;       /* Maximum confs returned     */
extern int gVibrate;       /* vibrate mode               */   
extern int gMini;          /* minimize mode              */

extern int gVerbose;       /* verbose mode               */   
extern int DEBUG;          /* DEBUG mode                 */

extern void delimiter(FILE *f);
extern void argstrPrmtrsUsage(int argc,char *argv[]);
extern void argstrPrmtrsSummary(int argc,char *argv[]);
extern void parseargstr(int argc,char *argv[]);

#endif
