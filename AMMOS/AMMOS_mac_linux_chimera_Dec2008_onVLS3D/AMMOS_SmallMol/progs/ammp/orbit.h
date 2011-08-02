
/* 1center 1s */
#define Or1   0
/* 1center 1s arbitrary origin */
#define Or1o  1
/* 2 center (on a bond) 1s */
#define Or2   2
/* 2 center p orbit on a bond rather than defined by a plane  */
#define Or2p   7
/* 3 center (on an angle diagonal) 1s */
#define Or3   3
/* 4 center unsymmetric (one sided) on a plane normal (ie nitrogen) */
#define Or4s  4
/* 4 center symmetric (one sided) on a plane normal  */
#define Or4p  6
/* multicentered on atoms (i.e. pi orbital for benzene) */
#define Orm  5

#ifdef GRACELESS
/*#define exp(x) ((x) > -100.? exp(x): 0.)
*/
#define exp(x) ((x) > -500.? exp(x): 0.)
#endif


typedef struct {
        int osn ; /* orbital serial number */
        int type; /* how am i defined ? */
	int ncenter;
        ATOM *myatom;
        ATOM *a1,*a2,*a3,*a4,*a5; /* orbits to define geometry */
	float normal ; /* self-term in the normalization */
        float along,along2 ;  /* fraction along geometry direction */
        float x,y,z; /* origin offset or normal direction */
        float rx[6],ry[6],rz[6];  /* my origin in full space */
        float rx1,ry1,rz1;  /* second origin in full space (for Or4p) */
        int spin,phase; /* spin (any) and phase (of Or4p only ) */
	int ipair;
        int n;  /* how many terms */
	int active;
        float a[6],r[6]; /* my weights and radii for gaussians */
	float rl[6];
/* couple is the orbits which are in the same spin
* system.  Normally the orbital is only itself, but
* for sp2 systems this may be larger
* if ncouple is < 0  it is a paired orbital (a basis member)
*/
	void *couple[10];
	int  ncouple;
        void *next;
	void *gang;
        }  ORBITAL ;


