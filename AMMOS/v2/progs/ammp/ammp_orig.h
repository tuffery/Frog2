/* ammp.h
*  include file for ammp and routines associated with it
*
*  define the ATOM structure
*/
/*  excluded stuff is for a list of bonded and 1-3 atoms not to be
*   used in the non-bonded evaluation 
*   the purpose of this is to speed up the nonbonded calculations since
*   the bond list varies only once and a while
*/
/* shadow version use 4dimensions (6 may be best) */
/* dont use double precision */
/* do use it 
#define float double
*/

/*#define NEXCLUDE 32 */
/*#define NEXCLUDE 64  */
#define NEXCLUDE 128  

/*
#define QUARTIC
#define  QUINTIC 
#define CUBIC  
*/

#ifdef QUARTIC
#define CUBIC 
#endif

#ifdef QUINTIC
#define CUBIC 
#define QUARTIC
#endif

/*#define NCLOSE 700 */
/*#define NCLOSE 500  */
/*#define NCLOSE 250   */
#define NCLOSE 400  
/*#define NCLOSE 100 */
/* 100 is enough for a 6A sphere  (usually) */

/* on the SGI use a different malloc */
#ifdef SGI
#include <sys/types.h>
#include <malloc.h>
#endif
#ifdef DECALPHA
#include <malloc.h>
#endif
#define critical_precision double

typedef struct{
	float x,y,z,w;
	critical_precision fx,fy,fz,fw;
	int serial;
	float q,a,b,mass;
	float na; /* atomic number == 0 if no orbitals */
	float rdebye;
	void *next;
	char active;
	char name[9];
	float chi,jaa;
	float vx,vy,vz,vw,dx,dy,dz,dw;
	float gx,gy,gz,gw;
	float VP,px,py,pz,pw,dpx,dpy,dpz,dpw; 
/*      float dpxx,dpxy,dpxz,dpyy,dpyz,dpzz; */
/* place holders for interpolation on V */
	float qxx,qxy,qxz,qyy,qyz,qzz;
	float qxw,qyw,qzw,qww;
#ifdef CUBIC
	float qxxx,qxxy,qxxz,qxyy,qxyz,qxzz;
	float qyyy,qyyz,qyzz,qzzz;
#endif
#ifdef QUARTIC
	float qxxxx,qxxxy,qxxxz,qxxyy,qxxyz,qxxzz;
	float qxyyy,qxyyz,qxyzz,qxzzz;
	float qyyyy,qyyyz,qyyzz,qyzzz,qzzzz;
#endif
#ifdef QUINTIC
	float qxxxxx,qxxxxy,qxxxxz,qxxxyy,qxxxyz,qxxxzz;
	float qxxyyy,qxxyyz,qxxyzz,qxxzzz;
	float qxyyyy,qxyyyz,qxyyzz,qxyzzz,qxzzzz;
	float qyyyyy,qyyyyz,qyyyzz,qyyzzz,qyzzzz,qzzzzz;
#endif
/* interpolation on force */
	void *close[NCLOSE];
	void *excluded[NEXCLUDE];
	char exkind[NEXCLUDE];
/* bitmap way to do it
#if (NEXCLUDE <33 )
	unsigned exkind;
#else
	unsigned long exkind;
#endif
*/
	int  dontuse;
	} ATOM;



#include "numeric.h" 
