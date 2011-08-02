#ifndef __types_h__
#define __types_h__

#define DcMONOCONF         0
#define DcREDUCEDMULTICONF 1
#define DcFULLMULTICONF    2

#define DcAMMOSMINICMD "bin/AMMOS_SmallMol_sp4.py"
#define DcAMMP_ENE_CMD "bin/AMMOS_SmallMol_ene.py"
#define DcAMMP_BOND_ENE_CMD "bin/AMMOS_SmallMol_bond_ene.py"
#define DcAMMP_NONBOND_ENE_CMD "bin/AMMOS_SmallMol_nonbond_ene.py"

typedef char DtChar;
typedef DtChar 		DtString[BUFSIZ];

typedef double DtFloat;
typedef double DtPoint3[3];
typedef double DtPoint4[4];
typedef char   DtStr3[4];
typedef char   DtStr4[5];
typedef char   DtStr5[6];

typedef DtFloat            DtMatrix4x4[4][4];

typedef struct satmLine {
  char head[30];
  double x, y, z;
  char tail[50];
} DtAtmLine;

typedef struct sbonds {
  int  *nBonds;
  int **bonds;
} DtBonds;

typedef struct sMol2 {
  int    nHead;
  int    nAtoms;
  int    nBonds;
  int    nTail;
  int   *index;
  char **head;
  char **atoms;
  char **bonds;
  char **tail;
  DtAtmLine *atomDetails;
  DtBonds    bondDetails;
} DtMol2;

typedef struct sMolArray {
  int nMol;
  DtMol2 *pM;
} DtMolArray;

#define DTOR(angle) 	   ((angle) * 0.017453292519943295)
#define RTOD(angle) 	   ((angle) * 57.295779513082323)
#define PI 3.14159265358979323846
#define dPI (PI*2.)

#define DcBONDENE    1
#define DcNONBONDENE 2
#define DcFULLENE    3

typedef struct structConf {
  double    energie; // Energy (internal to frog)
  double    fullene; // external energy (e.g. AMMP)
  double    bondene;    // external energy (e.g. AMMP) valence terms
  double    nonbondene; // external energy (e.g. AMMP) non bonded terms
  double   *dValues; // Dihedral actual values
  double   *values;  // Dihedral expected values
  int      *cvals;   // canonical conformation associated (among conformers per dihedral)
  DtPoint3 *coords;  // Atom coordinates
  struct structConf *nextConf; 
} DtStructConf;

#if 0
typedef struct coordsMol {
  double energie;
  DtPoint3 *coords;
  struct coordsMol *nextConf;
} coordsMol;
#endif

typedef struct DtStructFirstAndLastFilo {
  DtStructConf *first;
  DtStructConf *outted;
} DtStructFirstAndLastFilo;

typedef struct diedre {
  int     extremites[4]; // 4 atoms defining the dihedral
  double  value;      // Initial dihedral angle value in conformation
  double  curval;     // Actual dihedral angle value in conformation
  int     nbVals;
  int     nbAtoms2move;
  int     nbAtomsNot2move;
  double *valeurs;  // Canonical values
  double *nrjs;     // Canonical energy values
  int    *atoms2move;
  int    *atomsNot2move;
  int     isPepLike;
} diedre;

typedef struct energyFactors{
  double   *rEtoile;
  double  **rEtoileInter;
  double  **rEtoileInterUnZeroSept;
  double  **rEtoileInterZeroZeroSept;
  double  **rEtoileInterPow7;
  double  **rEtoileInterPow7UnDouze;
  double  **rEtoileInterPow7ZeroDouze;
  double  **epsilonInter;
  double  **eStaticCtInter;
  double ***energiesInter;
  double    energieTot;
} DtEnergyFactors;

typedef struct sring{
  int       nbAtoms;
  int      *atoms;
} DtRing;

typedef struct sinfo { // Information from .data file
  int       nbAtoms;
  DtPoint3 *coords;
  //  double  **valeurs; // chrg + 4 vdwparams  OBSOLETE
  double   *chrg;
  DtPoint4 *vdwparms;
  int     **nonLies;
  int       nbDiedres;
  diedre   *diedres;
  double   *probDiedres;
  int      *actuelsConfs;
  int       nbTotConfs;
  int      *nbConfsMax;
  int       nbRings;
  DtRing   *rings;
  int      *atmRing;
} DtInfos;

#endif
