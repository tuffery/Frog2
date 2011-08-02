/* include file for uff (modified ) preammp
*  
*  modifications form original include a,b terms rather than
* well depth and vdw terms
*  atom types effective masses and charges are read in with the atoms
*  in the residue dictionary,  charge and mass are
*   optional and a default value
*   will be used unless specified.
*
*  unified atoms   c2u c3u  included for unified atoms
*
* this file has the  data structures defined in it 
*
*/
/* copyright 1993 robert w. harrison
*  this notice may not be removed
*  the author(s) and his sponsoring institution (TJU)
*  reserve all commercial rights including the right 
*  to modify this notice.  This code may be
*  freely distributed for non-commercial scientific
*  use.
*/

#define NAME_LENGTH 8

typedef struct {
	char type[NAME_LENGTH+1];
	float r;  /* combining bond radius */
	float theta;  /* combining angle term */ 
	float a,b;  /* a,b terms for vdw */
	float Z;   /* effective charge */
	float angle_inc;
	float X;  /* electronegativity */
	float jaa; /* self coloumb */
	float mass, charge; /* guess */
	float V,U; /* torsion V and U terms for estimation */
	float hybrid; /* default hybrid force ignored if <= 0 */
	}   ATOMDATA;

/* data types for special bonds, and the like */
/* note the constants are treated symbolically */
#define SIG_DIGIT 10
typedef struct {
	char type1[NAME_LENGTH+1],type2[NAME_LENGTH+1];  /*atom types */
	int special;
	char r[SIG_DIGIT+1],k[SIG_DIGIT+1];              /* bond parameters */
	} BONDDATA;
typedef struct {
	char type1[NAME_LENGTH+1],type2[NAME_LENGTH+1],type3[NAME_LENGTH+1];  /*atom types */
	int special;
	char theta[SIG_DIGIT+1],k[SIG_DIGIT+1];              /* angle parameters */
	} ANGLEDATA;
typedef struct {
	char type1[NAME_LENGTH+1],type2[NAME_LENGTH+1],type3[NAME_LENGTH+1],type4[NAME_LENGTH+1];  /*atom types */
	int special;
	char hite[SIG_DIGIT+1],k[SIG_DIGIT+1];              /* hybrid parameters */
	} HYBRIDDATA;
typedef struct {
	char type1[NAME_LENGTH+1],type2[NAME_LENGTH+1],type3[NAME_LENGTH+1],type4[NAME_LENGTH+1];  /*atom types */
	int special;
	char k[SIG_DIGIT+1],off[SIG_DIGIT+1];              /* torsion parameters */
	char n[SIG_DIGIT+1] ;
	} TORSIONDATA;
