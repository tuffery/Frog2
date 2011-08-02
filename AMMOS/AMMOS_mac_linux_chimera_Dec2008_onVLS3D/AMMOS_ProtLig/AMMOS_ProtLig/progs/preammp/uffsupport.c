/* uffsupport.c
*
* functions which return uff force constants and the like
*
*
** the rules for geometry
*
**
**1) every pair of bonds with one common atom determine an angle
**2) every set of three bonds where 1 and 3 share half of 2 
**  i.e.  a-b, b-c, c-d
** defines a torsion.
**3) special chiral centers are hand coded for each residue
**4) Planar hybrids are hand coded for unified atom models, but
**  in general occur when three bonds exist with a common c_2 or o_2
**  atom.
**5) inversion hybrids are generally ignored, but are defined for a
**   a_n atom when n bonds exist with a as a common atom.
*
*/
#define ANSI 1
/* misc includes - ANSI and some are just to be safe */
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#ifdef ANSI
#include <stdlib.h>
#endif
#include "uff.h"
/* uff.h defines NAME_LENGTH == length of an atom id 
*              (Name_length+1 space reserved for \0)
*   and SIG_DIGIT == size of number or identifier in special
*    structure (SIG_DIGIT+1 space reserved for \0)
*/ 

/* copyright 1993 robert w. harrison
*  this notice may not be removed
*  the author(s) and his sponsoring institution (TJU)
*  reserve all commercial rights including the right
*  to modify this notice.  This code may be
*  freely distributed for non-commercial scientific
*  use.
*/

/* this below is not great but is the
* easiest way to do i */
/* atom types are defined globally */
#define MAX_ATOM 300
extern ATOMDATA  kinds[MAX_ATOM];
extern int inkinds;


float uff_bondlength( t1,t2,order )
int t1,t2;
float order;
{
float rad,rbo,ren;
rad = kinds[t1].r + kinds[t2].r;
/* 0.1332 is a magic number !! */
rbo = -0.1332* rad * log( (float)order);
ren =  sqrt(kinds[t1].X) -sqrt(kinds[t2].X) ;
ren = ren*ren* kinds[t1].r *kinds[t2].r
	/(kinds[t1].r*kinds[t1].X+kinds[t2].X*kinds[t2].r); 
return ( rad + rbo + ren);

}

float uff_bondforce( t1,t2,rij )
int t1,t2;
float rij; 
{
rij = rij*rij*rij;
return ( 664.12*kinds[t1].Z*kinds[t2].Z/rij/2);
/* the paper gives 664.12 shouldn't it be 2*332.17752 = 664.35504? */
}

/* angle force is centered on t2 */
float uff_angleforce( t1,t2,t3,rij,rjk )
int t1,t2,t3;
float rij,rjk;
{
float rik; /* j is t2 */
float cth;
float kijk;

cth = cos( kinds[t2].theta*3.14159265/180. );
rik = sqrt(rij*rij + rjk*rjk -2*rij*rjk*cth);
/* kijk = rij*rjk; */
kijk = 1.;
kijk = kijk*( 3*kijk*(1.-cth*cth)/rik/rik - cth)/rik/rik/rik;
kijk = kijk*664.12*kinds[t1].Z*kinds[t3].Z + 2*kinds[t2].angle_inc;
return(kijk*.5);
}
