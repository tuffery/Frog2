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
extern  int bondlist[1000][2]; /* more than enough */
extern  float rijlist[1000];
extern  float orderlist[1000];


do_bonds( data,nbond, output)
FILE *data,*output;
int nbond;
{
int i,j,ibond;
char line[256];
char token[5][80];
int match_label(),i1,i2;
float iord;
float uff_bondlength(), uff_bondforce();

for( ibond = 0; ibond< nbond; ibond++)
{ /* start of ibond loop */
if( fgets( line,256,data) == NULL)
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

if( strcmp( &token[0][0] ,"bond" ) == 0)
{ /* its a special bond  with user definition */
	i1 = match_label( &token[1][0]);
	i2 = match_label( &token[2][0]);
	fprintf( output," bond %d %d %s %s ;\n",
	serial[i1],serial[i2],&token[3][0],&token[4][0]); 
	bondlist[ibond][0] = i1;
	bondlist[ibond][1] = i2;
	for( j=0; j< ibond; j++)
	{
	if(( bondlist[j][0] == i1 && bondlist[j][1] == i2)
	 || ( bondlist[j][1] == i1 && bondlist[j][0] == i2))
	{
	fprintf(stderr," warning - bond %s %s is duplicated\n",
	&token[1][0],&token[2][0]);
	}
	}
/* make up an appx bond length */
	rijlist[ibond] = uff_bondlength( mykind[i1],mykind[i2],1.);
	orderlist[ibond] = 1.;
	
} else { /* otherwise its a normal bond */
	i1 = match_label( &token[0][0]);
	i2 = match_label( &token[1][0]);
	sscanf( &token[2][0],"%f",&iord);
	if( iord < 1) iord = 1;
	orderlist[ibond] = iord;
	bondlist[ibond][0] = i1;
	bondlist[ibond][1] = i2;
	for( j=0; j< ibond; j++)
	{
	if(( bondlist[j][0] == i1 && bondlist[j][1] == i2)
	 || ( bondlist[j][1] == i1 && bondlist[j][0] == i2))
	{
	fprintf(stderr," warning - bond %s %s is duplicated\n",
	&token[0][0],&token[1][0]);
	}
	}
	rijlist[ibond] = uff_bondlength( mykind[i1],mykind[i2],iord);
	fprintf( output," bond %d %d %f %f ;\n",
	serial[i1],serial[i2],rijlist[ibond],
	uff_bondforce( mykind[i1],mykind[i2],rijlist[ibond]));
}  /* end of bond kind if else then */


}/* end of ibond loop */
/*
*for(i=0; i< ibond; i++)
*{  printf(" %d %d %f\n",bondlist[i][0],bondlist[i][1],rijlist[i]);
*}
*/

} /* end of routine do_bonds*/

int match_label( who )
char *who;
{
 int i;
 char whole[80];
 char *cp;

 cp = who;
i = 0;
  while( *cp != '\0')
	{ whole[i++] = tolower(*cp++);}
 whole[i] = '\0';
 for( i=0; i< inres ; i++) /* inres,name are defined extern */
 {
	if( strcmp( &(whole[0]), &name[i][0]) == 0) return ( i);
}
	fprintf(stderr,"warning atom >%s< in geometry but not defined",who);
	return (-1);
}
