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
extern int hybridlist[1000][4],inhybridlist;


do_hybrid(dict, nhybrid, output )
FILE *dict,*output;
int nhybrid; /* the number of special angles can be zero */
{ /* start of routine */
int i,j,k,l,m;
int i1,i2,i3,i4, match_label();
char line[256],token[6][80];
int bonded[1000],inbonded;

inhybridlist = 0;
/* do the special angles if any */
for( m=0; m< nhybrid; m++)
{
if( fgets(line,256,dict) == NULL)
{ /* premature end of file is a fatal error */
fprintf(stderr," FATAL error, end of file in a dictionary ");
exit(0);
}
/* printf("%s",line); */
token[0][0] = '\0';
token[1][0] = '\0';
token[2][0] = '\0';
token[3][0] = '\0';
token[4][0] = '\0';
token[5][0] = '\0';
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


i1 = match_label( &token[0][0]);
i2 = match_label( &token[1][0]);
i3 = match_label( &token[2][0]);
i4 = match_label( &token[3][0]);

hybridlist[inhybridlist][0] = i1;
hybridlist[inhybridlist][1] = i2;
hybridlist[inhybridlist][2] = i3;
hybridlist[inhybridlist][3] = i4;
inhybridlist += 1;

 fprintf( output,"hybrid %d %d %d %d %s %s ;\n",
serial[i1],serial[i2],serial[i3],serial[i4],&token[4][0],&token[5][0]);
}/* end of special hybrid loop */

/* now for every atom find those atoms bonded to it */
for( i=0;i< inres; i++)
{/* for each atom in the residue */
i1 = i;
/*if( kinds[mykind[i1]].theta > 118. &&
    kinds[mykind[i1]].theta < 122.)
*/
if( kinds[mykind[i1]].hybrid > 0. )
{/* it is possible to do me */
inbonded = 0;
for( k=0; k< inbondlist; k++)
{  if( bondlist[k][0] == i)
	{ bonded[inbonded] = bondlist[k][1]; 
	inbonded += 1;}
  if( bondlist[k][1] == i)
	{ bonded[inbonded] = bondlist[k][0]; 
	inbonded += 1;}
/* duplicate checking -should never find any */
/*
	for( j=0; j< inbonded-2; j++)
	 if( bonded[j] == bonded[inbonded] ) {inbonded -= 1; break;}
*/
}
/*
for( j=0; j< inbonded; j++)
printf(" %d ",bonded[j]);
printf(" to %d %d\n",i1,inres);
*/
if( inbonded == 3)
{ /* only if fully bonded */
i2 = bonded[0];
i3 = bonded[1];
i4 = bonded[2];


for( k=0;k<inhybridlist;k++)
{
if( i1 == hybridlist[k][0])
{
if( i2 == hybridlist[k][1] || i2 == hybridlist[k][2]
	|| i2 == hybridlist[k][3] )
{
if( i3 == hybridlist[k][1] || i3 == hybridlist[k][2]
	|| i3 == hybridlist[k][3] )
{
if( i4 == hybridlist[k][1] || i4 == hybridlist[k][2]
	|| i4 == hybridlist[k][3] ) 
	{ i2 = -1; break; }
}
}
}
}


if( i2 >= 0)
{
fprintf( output," hybrid %d %d %d %d %f %f;\n",
serial[i2],serial[i3],serial[i4],serial[i1],kinds[mykind[i1]].hybrid, 0.);
}

}  /* only if fully bonded */
}/* end of 120 degree angle check */

} /* end of the for each atom in the residue (i) loop */


} /* end of do_hybrid */

