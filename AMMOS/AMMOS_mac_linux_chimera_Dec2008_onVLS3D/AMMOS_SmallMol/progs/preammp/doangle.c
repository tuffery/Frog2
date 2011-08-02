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

/* atom types are defined globally */
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

extern BONDDATA sbond[MAX_SP];

/* the array for the bonds */
extern  int bondlist[1000][2],inbondlist; /* more than enough */
extern  float rijlist[1000];
extern int anglelist[1000][3], inanglelist;
extern float aijklist[1000];


do_angle( dict, nangle, output )
FILE *dict,*output;
int nangle; /* the number of special angles can be zero */
{ /* start of routine */
int i,j,k,l;
int i1,i2,i3, match_label();
float uff_angleforce();
char line[256],token[5][BUFSIZ];
int bonded[1000],inbonded;
float rbond[1000];
inanglelist = 0;

/*
for( i=0; i< inbondlist ; i++)
{ printf( "%d %d %f\n", bondlist[i][0],bondlist[i][1],rijlist[i]); }
*/


/* do the special angles if any */
for( l=0; l< nangle; l++)
{
if( fgets(line,256,dict) == NULL)
{ /* premature end of file is a fatal error */
fprintf(stderr," FATAL error, end of file in a dictionary ");
exit(0);
}
/*printf("%s",line); */
token[0][0] = '\0';
token[1][0] = '\0';
token[2][0] = '\0';
token[3][0] = '\0';
token[4][0] = '\0';
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

i1 = match_label( &token[0][0]);
i2 = match_label( &token[1][0]);
i3 = match_label( &token[2][0]);

anglelist[inanglelist][0] = i1;
anglelist[inanglelist][1] = i2;
anglelist[inanglelist][2] = i3;
aijklist[inanglelist] = -1.;
inanglelist += 1;

 fprintf( output,"angle %d %d %d %s %s ;\n",
serial[i1],serial[i2],serial[i3],&token[3][0],&token[4][0]);
}/* end of special angle loop */
/* now for every atom find those atoms bonded to it */
for( i=0;i< inres; i++)
{/* for each atom in the residue */
inbonded = 0;
for( k=0; k< inbondlist; k++)
{  if( bondlist[k][0] == i)
	{ bonded[inbonded] = bondlist[k][1];
	  rbond[inbonded] = rijlist[k];
	inbonded += 1;}
  if( bondlist[k][1] == i)
	{ bonded[inbonded] = bondlist[k][0];
	  rbond[inbonded] = rijlist[k];
	inbonded += 1;}
/* duplicate checking -should never find any */

/*
	for( j=0; j< inbonded-1; j++)
	 if( bonded[j] == bonded[inbonded] ) {inbonded -= 1; break;}
*/

}
if( inbonded > 1 )
{/* if inbonded == 1 there are no angle centered on me */

/*
for( j=0; j< inbonded; j++)
{
printf(" %d %d %f\n",j, bonded[j],rbond[j]);
}	
*/



for( j=0; j< inbonded; j++)
{/*j*/
i2 = i; i1 = bonded[j];

/*for( k=j+1; k< inbonded; k++)*/
for( k=0; k< inbonded; k++)
{ /*k*/
i3 = bonded[k];
if(k == j ) i3 = -1;
/* even though this process will never generate an overlaping angle
*  we still must check because of specials */
for( l=0; l < inanglelist; l ++)
{/*l*/
if( i2 == anglelist[l][1] )
{/* if not centered cannot be equal !!! */
if( ( i1 == anglelist[l][0] && i3 == anglelist[l][2]) || 
    ( i1 == anglelist[l][2] && i3 == anglelist[l][0]) )
{
	i3 = -1; break;
}
}
}/*l*/
if( i3 > 0 ) /* its oK */
{
anglelist[inanglelist][0] = i1;
anglelist[inanglelist][1] = i2;
anglelist[inanglelist][2] = i3;
aijklist[inanglelist] = 
	uff_angleforce(mykind[i1],mykind[i2],mykind[i3],rbond[j],rbond[k]);
fprintf( output," angle %d %d %d %f %f ;\n",
serial[i1],serial[i2],serial[i3], 
aijklist[inanglelist],
kinds[mykind[i2]].theta);
inanglelist += 1;

}

}/*k*/
}/*j*/
}/* end of if inbonded > 1 */

} /* end of the for each atom in the residue (i) loop */


} /* end of do_angle */

