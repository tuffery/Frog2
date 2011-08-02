#ifndef __Mol2_h__
#define __Mol2_h__

// extern DtMolArray *lectMol2(char *fname, int verbose);
extern DtMolArray *lectMol2(char *fname, int maxRead, int verbose);
extern void lectMol2Crds(char *fname, DtStructConf *aConf, DtMol2 *pM, int verbose);
// extern void lectMol2Crds(char *fname, DtStructConf *aConf, DtMol2 *pM, int which, int verbose);
// extern void outMol2(FILE *f, DtMol2 *pM, char *cname, double energy, int useDetails);
extern void outMol2(FILE *f, DtMol2 *pM, char *cname, double energy, double fullEnergy, int useDetails);
extern void mapMol2Info(DtMol2 *pM, DtInfos *pI);
extern void setupCoords(DtMol2 *pM, DtPoint3 *crds);
extern DtInfos * info2Info(DtMol2 *pM, DtInfos *pI);

#endif
