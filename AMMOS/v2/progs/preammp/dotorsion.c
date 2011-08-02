/*
*  v 2  add ring awarness
*  3 rings have no torsions
*  4,5,... ? rings have increased torsions
*
* 3/5/03
*/
/* preammp.c
*
* prepare ammp files using UFF type potentials
*  dictionary entries for bonds
* 	
*          nbond  (general and special)
*          as  res.atom res.atom order
*          or  bond res.atom res.atom r k  for specials
*
*
**rules for geometry
**
**1) every pair of bonds with one common atom determine an angle
**2) every set of three bonds where 1 and 3 share half of 2 
**  i.e.  a-b, b-c, c-d
** defines a torsion.
**3) special chiral centers are hand coded for each residue
**4) Planar torsions are hand coded for unified atom models, but
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

/* atom kinds are defined globally */
#define MAX_ATOM 300
extern ATOMDATA  kinds[MAX_ATOM];
extern int inkinds;

/* maximum number of atoms is a residue */
#define MAX_RES 300

extern float x[MAX_RES], y[MAX_RES], z[MAX_RES];
extern float a[MAX_RES],b[MAX_RES],q[MAX_RES];
extern float mass[MAX_RES]; 
extern int serial[MAX_RES],defined[MAX_RES];
extern int mykind[MAX_RES];
extern char name[MAX_RES][20];
extern int inres;

/* space for special stuff */
#define MAX_SP 300


/* the array for the bonds */
extern  int bondlist[1000][2],inbondlist; /* more than enough */
extern  float rijlist[1000],orderlist[1000];
extern int torsionlist[1000][4],intorsionlist;


static int ringlist[1000];
static int anglelist[1000][3],inanglelist;
/* this anglelist is different from the extern one in that it is the
* list of all the angles that are generated from bonds */


do_torsion(dict, ntorsion, output )
FILE *dict,*output;
int ntorsion; /* the number of special angles can be zero */
{ /* start of routine */
int i,j,k,l,m;
int right, left ;
int i1,i2,i3,i4, match_label();
char line[256],token[7][80];
int bondedr[1000],inbondedr;
int bondedl[1000],inbondedl;
int n;
float offset,kijkl;
/* initialize and build  the ring lists */
	init_ringlist();
	init_anglelist();
	three_four_five_ringlist();

intorsionlist = 0;
/* do the special angles if any */
for( l=0; l< ntorsion; l++)
{
if( fgets(line,256,dict) == NULL)
{ /* premature end of file is a fatal error */
fprintf(stderr," FATAL error, end of file in a dictionary ");
exit(0);
}
/* now parse out the first term */
i = 0;
while((line[i] == ' ' || line[i]  == '\t') && line[i] != '\0' ) i++;

j = 0;
while( (line[i] != ' '  && line[i] != ',') && (line[i] != '\0'
 && line[i] != '\n')) token[0][j++] = line[i++];
token[0][j] = '\0';
while((line[i] == ' ' || line[i]  == '\t') && line[i] != '\0' ) i++;
j = 0;
while( (line[i] != ' '  && line[i] != ',') && (line[i] != '\0'
 && line[i] != '\n')) token[1][j++] = line[i++];
token[1][j] = '\0';
while((line[i] == ' ' || line[i]  == '\t') && line[i] != '\0' ) i++;
j = 0;
while( (line[i] != ' '  && line[i] != ',') && (line[i] != '\0'
 && line[i] != '\n')) token[2][j++] = line[i++];
token[2][j] = '\0';
while((line[i] == ' ' || line[i]  == '\t') && line[i] != '\0' ) i++;
j = 0;
while( (line[i] != ' '  && line[i] != ',') && (line[i] != '\0'
 && line[i] != '\n')) token[3][j++] = line[i++];
token[3][j] = '\0';
while((line[i] == ' ' || line[i]  == '\t') && line[i] != '\0' ) i++;
j = 0;
while( (line[i] != ' '  && line[i] != ',') && (line[i] != '\0'
 && line[i] != '\n')) token[4][j++] = line[i++];
token[4][j] = '\0';

while((line[i] == ' ' || line[i]  == '\t') && line[i] != '\0' ) i++;
j = 0;
while( (line[i] != ' '  && line[i] != ',') && (line[i] != '\0'
 && line[i] != '\n')) token[5][j++] = line[i++];
token[5][j] = '\0';
j = 0;
while( (line[i] != ' '  && line[i] != ',') && (line[i] != '\0'
 && line[i] != '\n')) token[6][j++] = line[i++];
token[6][j] = '\0';

i1 = match_label( &token[0][0]);
i2 = match_label( &token[1][0]);
i3 = match_label( &token[2][0]);
i4 = match_label( &token[3][0]);

torsionlist[intorsionlist][0] = i1;
torsionlist[intorsionlist][1] = i2;
torsionlist[intorsionlist][2] = i3;
torsionlist[intorsionlist][3] = i4;
intorsionlist += 1;

 fprintf( output,"torsion %d %d %d %d %s %s %s ;\n",
serial[i1],serial[i2],serial[i3],serial[i4],&token[4][0],&token[5][0],
&token[6][0]);
}/* end of special torsion loop */

for( i=0; i< inbondlist; i++)
{ /* for each bond */
inbondedl = 0;
for( k=0; k< inbondlist; k++)
{  if( k != i)
	{
	if( bondlist[k][0] == bondlist[i][1])
	{ bondedl[inbondedl] = bondlist[k][1]; 
	inbondedl += 1;}
  if( bondlist[k][1] == bondlist[i][1])
	{ bondedl[inbondedl] = bondlist[k][0]; 
	inbondedl += 1;}
	}
}
inbondedr = 0;
for( k=0; k< inbondlist; k++)
{  if( k != i)
	{
	if( bondlist[k][0] == bondlist[i][0])
	{ bondedr[inbondedr] = bondlist[k][1]; 
	inbondedr += 1;}
  if( bondlist[k][1] == bondlist[i][0])
	{ bondedr[inbondedr] = bondlist[k][0]; 
	inbondedr += 1;}
	}
}

/* skip if terminal or theta == 180 */
/*if( inbondedr > 0 && inbondedl > 0
 && (kinds[mykind[right]].theta <179. && kinds[mykind[left]].theta< 179.))
*/
if( inbondedr > 0 && inbondedl > 0 )
{ /* only do torsions when there are some */
/* sp2 atoms do not possess kinds.V > 0 */
/* sp1 atoms have 180 as theta */
right = bondlist[i][0];
left = bondlist[i][1];  /* note right and left are arb. */
if( kinds[mykind[right]].V > 0. &&
	kinds[mykind[left]].V > 0.)
	{ /* sp3 sp3 */
	n = 3; offset = 0.;
	kijkl = sqrt(kinds[mykind[left]].V*kinds[mykind[right]].V);
	if( ringlist[i] == 3) kijkl *= 6.;
	if( ringlist[i] == 4) kijkl *= 3.;
	if( ringlist[i] == 5) kijkl *= 1.5;
	} else 
if( (kinds[mykind[right]].V > 0. &&
	kinds[mykind[left]].V < 0.) ||
( kinds[mykind[right]].V < 0. &&
	kinds[mykind[left]].V > 0.)  )
	{ /* sp3 sp2 */
/*  n is complicated
*
*  sp3 sp3 sp2 x  n=3  ,180.
* sp2 sp3 sp2 x n = 6  ,180.
*/
	n = 3; offset = 180.;
	for( l=0; l< inbondedr ;l++)
	for( k=0; k< inbondedl ;k++)
	{
		i1 = bondedr[l];
		i4 = bondedl[k];
	if((kinds[mykind[right]].V > 0. && kinds[mykind[i1]].V < 0. )
			|| (kinds[mykind[left]].V> 0. &&kinds[mykind[i4]].V < 0.))
	n = 6; offset = 180.;
	}
	kijkl = 1.;
	kijkl = sqrt(kinds[mykind[left]].U*kinds[mykind[right]].U);
	kijkl = 5.*kijkl*(1.+ 4.18*log( orderlist[i]));
	} else 
if( kinds[mykind[right]].V < 0. &&
	kinds[mykind[left]].V < 0.)
	{ /* sp2 sp2 */
	n = 2; offset = 180.;
	kijkl = sqrt(kinds[mykind[left]].U*kinds[mykind[right]].U);
	kijkl = 5.*kijkl*(1.+ 4.18*log( orderlist[i]));
	} 

if( kinds[mykind[left]].theta > 170. ) kijkl = 0.;
if( kinds[mykind[right]].theta > 170. ) kijkl = 0.;

	kijkl = kijkl/inbondedr/inbondedl;
	i2 = right;
	i3 = left;
	for( l=0; l< inbondedr ;l++)
	for( k=0; k< inbondedl ;k++)
	{
		i1 = bondedr[l];
		i4 = bondedl[k];
fprintf( output," torsion %d %d %d %d %f %d %f ;\n",
	serial[i1],serial[i2],serial[i3],serial[i4],
	kijkl,n,offset);
	}



} /* end of inbondedr > 0 if */
}/* end of bond loop */

} /* end of do_torsion */


int init_ringlist()
{
	int i;
	for( i=0; i< inbondlist; i++)
		ringlist[i] = 0;
}

int init_anglelist()
{
	int i,j;
	inanglelist = 0;
	for( i=0; i< inbondlist; i++)
	{
		int i1,i2,j1,j2,j3;
		i1 = bondlist[i][0];
		i2 = bondlist[i][1];
		for( j=0; j< inbondlist; j++)
		{  if( j == i) continue; 
		j1 = bondlist[j][0];
		j2 = bondlist[j][1];

		if( j1 == i1)
		{ 
			anglelist[inanglelist][0] = j2;
			anglelist[inanglelist][1] = i1;
			anglelist[inanglelist++][2] = i2;
		}
		if( j2 == i1)
		{ 
			anglelist[inanglelist][0] = j1;
			anglelist[inanglelist][1] = i1;
			anglelist[inanglelist++][2] = i2;
		}
		if( j1 == i2)
		{ 
			anglelist[inanglelist][0] = j2;
			anglelist[inanglelist][1] = i2;
			anglelist[inanglelist++][2] = i1;
		}
		if( j2 == i2)
		{ 
			anglelist[inanglelist][0] = j1;
			anglelist[inanglelist][1] = i2;
			anglelist[inanglelist++][2] = i1;
		}
		}/* j */
	}/* i */
}

int three_four_five_ringlist( )
{
	int i,j,k,l, i1,i2,j1,j2;
	int rightbond[100],leftbond[100];
	int inright,inleft;
	for( i=0; i< inbondlist; i++)
	{
/* build the right/left bond list */
	inright = 0;
	inleft = 0;
		i1 = bondlist[i][0];
		i2 = bondlist[i][1];
		for( j=0; j< inbondlist; j++)
		{
		if( j == i) continue;
		j1 = bondlist[j][0];
		j2 = bondlist[j][1];
		
		if( j1 == i1) rightbond[inright++] = j2;
		if( j2 == i1) rightbond[inright++] = j1;
		if( j1 == i2) leftbond[inleft++] = j2;
		if( j2 == i2) leftbond[inleft++] = j1;
		}/* j */
/* no rings if no bonds */
		if( inleft == 0 || inright == 0) continue;
/* now to check for a threebond there is an identical atom in right and left */
		for( k=0; k< inright; k++)
		for( l = 0; l < inleft; l++)
		{  if( rightbond[k] == leftbond[l] ) 
			{ ringlist[i] = 3; break; }
		}/* k,l */

		if( ringlist[i] > 0 ) continue; /* cant be 3 and 4 and ... */
		for( k=0; k< inright; k++)
		for( l = 0; l < inleft; l++)
		{
			int m;
			for(m=0; m < inbondlist; m++)
			{
			if( m == i) continue;
			j1 = bondlist[m][0];
			j2 = bondlist[m][1];
			if( (rightbond[k] == j1 && leftbond[l] == j2) ||
			    (rightbond[k] == j2 && leftbond[l] == j1) )
			{
				ringlist[i] = 4; break;
			}
			}/* m */
			if( ringlist[i] == 4) break;
		}/* k,l */
		if( ringlist[i] > 0 ) continue; /* cant be 3 and 4 and ... */
		for( k=0; k< inright; k++)
		for( l = 0; l < inleft; l++)
		{
			int m;
			for(m=0; m < inanglelist; m++)
			{
			j1 = anglelist[m][0];
			j2 = anglelist[m][2];
			if( (rightbond[k] == j1 && leftbond[l] == j2) ||
			    (rightbond[k] == j2 && leftbond[l] == j1) )
			{
				ringlist[i] = 5; break;
			}
			}/* m */	
		}/* k,l */
	}/* i */
}
