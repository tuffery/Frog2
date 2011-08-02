/* preammp.c
*
* prepare ammp files using UFF type potentials
* 
*  required files
*
*  atomdictionary  -  defines atom dictionary 
*                  - default masses and charges
*                  - default a,b terms read as distance and minimum
*                  - data is used to calculate the bonds...
*   per line type (%s), r,theta,sigma,emin,Z,angle_inc,X,jaa,mass,charge,V,U,hybrid
*    r = combining radius, theta = angle in degrees
*   sigma = vdw radius, emin = vdw depth, Z = effect charge, X = electronegati
*   jaa = self coloumb term
*    mass = mass, Charge = default charge, v,u are torsion terms
*       v > 0 , sp3, v< 0 sp2,sp1
*   hybrid is default hybrid  (c2,n2,ca default of 6 was used is not good)
*
*  residuedictionary - one per residue
*                    lists atoms
*   
*	natom then natom records (generalized )
*    as '  name res.atom  type type mass mass charge charge a a b b sigma sigma e distance'
*  where name type mass charge a b sigma e  are keywords
* either a b or sigma e or none should be used
*   name and type are required the others are optional 
*                     the lists the bonds
*
*           r k theta... constants will be parsed as strings so tha
*           symbolic constants may be used
* 	
*          nbond  (general and special)
*          as  res.atom res.atom order
*          or  bond res.atom res.atom r k  for specials
*
*          nangle  (special only)
*             res.atom res.atom res.atom  theta k
*         
*          nhybrid  (special only)
*         res.atom res.atom res.atom res.atom hite k
*
*          ntorsion  (special only)
*         res.atom res.atom res.atom res.atom k,n,offset
*
* in as much as possible the following rules will be used
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
ATOMDATA  kinds[MAX_ATOM];
int inkinds;

/* maximum number of atoms is a residue */
#define MAX_RES 300

float x[MAX_RES], y[MAX_RES], z[MAX_RES];
float a[MAX_RES],b[MAX_RES],q[MAX_RES];
float mass[MAX_RES]; 
int serial[MAX_RES],defined[MAX_RES];
int mykind[MAX_RES];
char name[MAX_RES][20];
int inres;

/* space for special stuff */
#define MAX_SP 300
/*
*BONDDATA sbond[MAX_SP];
*ANGLEDATA sangle[MAX_SP];
*HYBRIDDATA shybrid[MAX_SP];
*TORSIONDATA storsion[MAX_SP];
*/

/* the array for the bonds */
 int bondlist[1000][2],inbondlist; /* more than enough */
 float rijlist[1000],orderlist[1000];
 int anglelist[1000][3],inanglelist;
 float aijklist[1000];
int hybridlist[1000][4], inhybridlist;
 int torsionlist[1000][4],intorsionlist;
int natoms,nbond,nangle,nhybrid,ntorsion;



int main()
{
char coordinate_input[80];
char ammp_output[80];
char dictionary[80];
char  work[1000],keep[1000]; 
char atype[20],aname[20];
char *fgets();

float sigma,emin;

int ifile,i,imreading,ii,myres;

FILE *in1,*dict,*output,*fopen();

	printf(" atom definition file:>\n");
/* strip out the file name and put it in the right place for opening */
        fgets( work,80,stdin );
        for(i= 0; i<80; i++)
                if( work[i] != ' ') break;
        for( ifile = i; ifile < 80 ; ifile++)
                {
                if( work[ifile] == ' ' ) break;
                if( work[ifile] == '\0' ) break;
                if( work[ifile] == '\n' ) break;
                dictionary[ifile -i ] = work[ifile];
                dictionary[ifile -i +1 ] = '\0';
                }
/* now read in the atoms */
	dict = fopen( dictionary, "r");
	if( dict == NULL)
	{ fprintf(stderr," sorry i can't open %s\n",dictionary); exit(0);}
	inkinds = 0;
	while ( fgets( work,256, dict) != NULL)
	{
	sscanf( work,"%s %f %f %f %f %f %f %f %f %f %f %f %f %f \n",
		&kinds[inkinds].type[0],
		&kinds[inkinds].r ,
		&kinds[inkinds].theta ,
		&sigma ,
		&emin ,
		&kinds[inkinds].Z ,
		&kinds[inkinds].angle_inc,
		&kinds[inkinds].X ,
		&kinds[inkinds].jaa ,
		&kinds[inkinds].mass ,
		&kinds[inkinds].charge ,
		&kinds[inkinds].V, &kinds[inkinds].U ,
		&kinds[inkinds].hybrid
		); 
/*		kinds[inkinds].theta *= 3.14159265/180.; 
*/
		sigma = sigma*sigma;
		sigma = sigma*sigma*sigma;
/*		kinds[inkinds].b = sqrt( 4.*emin*sigma*sigma);
		kinds[inkinds].a = sqrt( 4.*emin*sigma);
*/
		kinds[inkinds].b = sqrt( emin*sigma*sigma);
		kinds[inkinds].a = sqrt( 2*emin*sigma);
		inkinds ++;
		if( inkinds > MAX_ATOM) 
		{ fprintf( stderr," too many kinds of atoms in dictionary\n");
		  exit(0); }
	}			
	fclose( dict);

	printf(" Residue dictionary directory:>\n");
/* strip out the file name and put it in the right place for opening */
        fgets( work,80,stdin );
        for(i= 0; i<80; i++)
                if( work[i] != ' ') break;
        for( ifile = i; ifile < 80 ; ifile++)
                {
                if( work[ifile] == ' ' ) break;
                if( work[ifile] == '\0' ) break;
                if( work[ifile] == '\n' ) break;
                dictionary[ifile -i ] = work[ifile];
                dictionary[ifile -i +1 ] = '\0';
                }
	ifile = ifile -i ;
	if( ifile == 0 ){ dictionary[ifile++] = '.'; 
			dictionary[ifile] = '\0'; }
	if( dictionary[ifile-1] != '/') 
	{ ifile++; dictionary[ifile-1] = '/'; dictionary[ifile] = '\0';}
	printf(" atom input coordinate file:>\n");
/* strip out the file name and put it in the right place for opening */
        fgets( work,80,stdin );
        for(i= 0; i<80; i++)
                if( work[i] != ' ') break;
        for( ifile = i; ifile < 80 ; ifile++)
                {
                if( work[ifile] == ' ' ) break;
                if( work[ifile] == '\0' ) break;
                if( work[ifile] == '\n' ) break;
                coordinate_input[ifile -i ] = work[ifile];
                coordinate_input[ifile -i +1 ] = '\0';
                }

	printf(" AMMP file:>\n");
/* strip out the file name and put it in the right place for opening */
        fgets( work,80,stdin );
        for(i= 0; i<80; i++)
                if( work[i] != ' ') break;
        for( ifile = i; ifile < 80 ; ifile++)
                {
                if( work[ifile] == ' ' ) break;
                if( work[ifile] == '\0' ) break;
                if( work[ifile] == '\n' ) break;
                ammp_output[ifile -i ] = work[ifile];
                ammp_output[ifile -i +1 ] = '\0';
                }
	in1 = fopen( coordinate_input, "r");
	output = fopen( ammp_output, "w");
	if (in1 == NULL) goto NO_IN1_FILE;

/* read in the first residue
*  stuff is read into work when the residue number changes then
* work is copied into keep, so the first atom of any residue is
* going to be in the keep buffer.  to initialize we read the
* first atom into keep
*/
	keep[0] = '\0';
	while( !(keep[0] ==  'A' && keep[1] == 'T' && keep[2] == 'O' && keep[3] == 'M')
	 && !(keep[0] == 'H'&& keep[1] == 'E'
		&& keep[1] == 'T'&& keep[1] == 'A'&& keep[1] == 'T' ) )
	{  if( fgets(keep,90,in1) == NULL) {
		fprintf(stderr," no atoms in input file ? \n") ;
NO_IN1_FILE:
		printf(" enter the dictionary name >\n");
		fgets( keep, 90, stdin);
	i=0; while( (keep[i] == ' ' || keep[i] == '\t') && keep[i] != '\0') i++;	
	for( imreading = 0; imreading < 3; imreading++)
	{ aname[imreading] = toupper( keep[i + imreading] ); }
	aname[3] = '\0';
		 goto NO_ATOMS; 
					   }
	}

	imreading = 0;
	while( imreading == 0)
	{/* begining of imreading loop */
	inres = 0;
	for(i=0; i< MAX_RES; i++)
		{ mykind[i] = -1;
		defined[i] = -1;
		name[i][0] = '\0';
		x[i] = 0.;
		y[i] = 0.;
		z[i] = 0.;
		}
	sscanf( &keep[22],"%d",&myres);
        sscanf(&keep[29],"%f %f %f",&x[inres],&y[inres],&z[inres]);
        sscanf(&keep[11],"%s",&atype[0]);
        sscanf(&keep[17],"%s",&aname[0]);
	serial[inres] = 100*myres + inres;
        ii = 0;
        for( i=0; i< 3; i++ )
        {
                if( aname[i] == '\0') break;
                name[inres][ii++] = (char)tolower(aname[i]);
                }
        name[inres][ii++] = '.';
        for(i=0; i< 4; i++)
        {
                if( atype[i] == '\0' ) break;
                name[inres][ii++] = (char)tolower(atype[i]);

        }
        name[inres][ii] = '\0';
	a[inres] = 0.; b[inres] = 0.; q[inres] = 0.; mass[inres] = 1.;
	defined[inres] = 0;
	for( inres = 1; inres < MAX_RES; inres ++)
	{ /* start of reading an atom loop */

	if( fgets( work,90,in1) == NULL) {
		imreading = 1; break;    }
	if( work[0] != 'A' && work[0] != 'H') { 
		imreading = 1; break;    }
	

/*
	printf("%s\n",work);
	printf("%d\n",inres);
*/

	sscanf( &work[22],"%d",&ii);
	if( ii != myres )
	{ii = 0; while( work[ii] != '\0'){ keep[ii] = work[ii] ; ii++;}
/*	inres += 1; */
	 break;
	}
        sscanf(&work[29],"%f %f %f",&x[inres],&y[inres],&z[inres]);
        sscanf(&work[11],"%s",&atype[0]);
        sscanf(&work[17],"%s",&aname[0]);
	serial[inres] = 100*myres + inres;
        ii = 0;
        for( i=0; i< 3; i++ )
        {
                if( aname[i] == '\0') break;
                name[inres][ii++] = (char)tolower(aname[i]);
                }
        name[inres][ii++] = '.';
        for(i=0; i< 4; i++)
        {
                if( atype[i] == '\0' ) break;
                name[inres][ii++] = (char)tolower(atype[i]);

        }
        name[inres][ii] = '\0';
	a[inres] = 0.; b[inres] = 0.; q[inres] = 0.; mass[inres] = 1.;
	defined[inres] = 0;
	}/* end of reading an atom loop */
/*
	for( i=0; i< inres; i++)
	printf(">%s<\n",&name[i][0]); 
*/

/* if here then all of the atoms in the residue are read in */
/* now we have to open and read the dictionary file 
*  there could be atoms in either the residue or dictionary which 
*  are not there i.e. oxt for a residue or a missing H */
/* prepare the filename */
NO_ATOMS:
	i = 0;
	while(dictionary[i] != '\0') { work[i] = dictionary[i]; i++;}
	ii = 0;
	while(aname[ii] != '\0') {work[i] = aname[ii]; i++; ii++;}
	work[i++] = '\0';
	dict = fopen(work,"r");
	if( dict == NULL ){
		fprintf(stderr," dictionary %s not found \n", work);
		write_atoms(output);
		goto DONE;
			}
/* read the number of atoms */
	if(fgets( work,80,dict)==NULL)
	{ fprintf(stderr,"BAD DICTIONARY %s\n",aname); exit(0);}
	sscanf(work,"%d",&natoms);
	read_atoms( dict,natoms );
	write_atoms(output);
	if(fgets( work,80,dict)==NULL)
	{ fprintf(stderr,"BAD DICTIONARY %s\n",aname); exit(0);}
	sscanf(work,"%d",&nbond);
	do_bonds( dict,nbond,output);
	inbondlist = nbond;
	nangle = 0;
	if(fgets( work,80,dict)!=NULL)
	sscanf(work,"%d",&nangle);
	do_angle(dict,nangle,output);
	nhybrid = 0;
	if(fgets( work,80,dict)!=NULL)
	sscanf(work,"%d",&nhybrid);
	do_hybrid(dict,nhybrid,output);
	ntorsion = 0;
	if(fgets( work,80,dict)!=NULL)
	sscanf(work,"%d",&ntorsion);
	do_torsion(dict,ntorsion,output);


DONE:  
	fclose( dict );
	}/* end of imreading loop */
}/* end of module main */


write_atoms( fp )
FILE *fp;
{
/* write an AMMP atom command */
	int i;
	for( i=0; i< inres; i++)
	{
/*
*	printf( "atom %f %f %f %d %s %f %f %f %f ; \n",
*	x[i],y[i],z[i],serial[i], &name[i][0],
*	q[i],a[i],b[i],mass[i] ); 
*/
	fprintf( fp,"atom %f %f %f %d %s %f %f %f %f ; \n",
	x[i],y[i],z[i],serial[i], &name[i][0],
	q[i],a[i],b[i],mass[i] ); 
	if( mykind[i] > -1)
	if( kinds[mykind[i]].X > 0. && kinds[mykind[i]].jaa > 0.) 
	fprintf( fp,"mompar %d %f %f ;\n",
	serial[i], kinds[mykind[i]].X, kinds[mykind[i]].jaa);
	}
}/*end of routine write_atoms */

read_atoms( fp, n)
int n;
FILE *fp;
{
int i,j,k,l; 
char work[257];
char strings[21][40];
char *fgets();
float sigma, emin;
int indic[MAX_RES];
int theatom;

int inkey;
int found[9];
char *keyword[9];

keyword[0] = "name";
keyword[1] = "type";
keyword[2] = "mass";
keyword[3] = "charge";
keyword[4] = "a";
keyword[5] = "b";
keyword[6] = "sigma";
keyword[7] = "e";

inkey = 8;

for( j=0; j< MAX_RES; j++)
 indic[j] = -1;

for(j=0; j< n; j++)
{ /*j*/
	for(i=0; i< inkey ; i++) found[i] = -1;
	for( i=0; i< 20; i++)
	strings[i][0] = '\0';

	if( fgets(work,256,fp) == NULL )
	{  fprintf(stderr," corrupt dictionary - not enough atoms \n");
	   exit(0);
	}
/* lex out the strings */
	k = 0;
	l = 0;
	for( i=0; i< 255; i++) /* clean up first */
	{
		if(work[i] == ',') work[i] = ' ' ;
		if(work[i] == '\n') {k++; break ;}
		if(work[i] == '\0') {k++; break ;}
		work[i] = tolower(work[i]);
		if( l < 39 ) 
	{	strings[k][l] = work[i]; l+=1; strings[k][l] = '\0'; }
/*	if( work[i] == ' ') {printf("%s\n",&strings[k][0]);} */

		if( work[i] == ' ' ) 
		{ while(work[i] == ' ') {i++;} i-=1;
		strings[k][l-1] = '\0'; k++; l = 0;}
	}
/* now have the strings and we know how many (k) */
	for( l= 0; l < inkey ;l++)	
	{
	for( i=0; i< k; i++)
	{
	if( strcmp( keyword[l], &strings[i][0]) == 0)
	{ /*then we've got a hit */
	found[l] = i+1;
	break;
	}
	}/* i */
/*
	if( found[l] > 0)
	printf("%d %d %s\n",l,found[l],&strings[found[l]][0]);
	else
	printf("%d %d \n",l,found[l]);
*/

	}/* l */
/* now have the values associated with each keyword */
	if( found[0] < 0 || found[1] < 0) 
{ fprintf(stderr," corrupt dictionary  name and type fields required\n");
	  exit(0);}
/* check to see if the atom is in the file */ 
	for( i=0; i< inres; i++)
	{
	if( strcmp( &strings[found[0]][0], &name[i][0]) == 0)
	{
/*	printf("%d %s %s\n",i,&name[i][0],&strings[found[0]][0]);
*/
	indic[i] = 1000;
	theatom = i;
	break;	
	}
	}
/* if we're here and the following is true the atom has not been found */
	if( i == inres )
	{
	strcpy( &name[inres][0] , &strings[found[0]][0] );
	serial[inres] = serial[inres-1] + 1;
	theatom = inres;
	inres += 1;
	}
/* now find the kind out */
	for( i=0; i< inkinds ; i++)
	{
	if( strcmp( &strings[found[1]][0] ,  &kinds[i].type[0]) == 0 )
	{
/*
	printf("%d %s %s\n",i,&kinds[i].type[0],&strings[found[1]][0]);
*/
	mykind[theatom] = i;
	break;
	}
	}
	if(i == inkinds )
	{/* if here we havn't found it */
	fprintf(stderr,"atom %s kind %s is undefined\n", 
		&strings[found[0]][0],&strings[found[1]][0]);
	fprintf(stderr," sorry - must exit ");
	exit(0);
	}
/* put in the dictionary values */
	a[theatom] = kinds[mykind[theatom]].a;
	b[theatom] = kinds[mykind[theatom]].b;
	q[theatom] = kinds[mykind[theatom]].charge;
	mass[theatom] = kinds[mykind[theatom]].mass;
/* now we check the optional parameters */
/* keyword[0] = "name";
*  keyword[1] = "type";
*  keyword[2] = "mass";
*  keyword[3] = "charge";
*  keyword[4] = "a";
*  keyword[5] = "b";
*  keyword[6] = "sigma";
*  keyword[7] = "e";
*/
	if( found[3] > 0) 
		sscanf(&strings[found[3]][0],"%f",&q[theatom]);
	if( found[2] > 0) 
		sscanf(&strings[found[2]][0],"%f",&mass[theatom]);
/* sigma,e representation of vdw ( 6,7) */
	if( found[6] > 0 && found[7] > 0)
	{
		sscanf(&strings[found[6]][0],"%f",&sigma);
		sscanf(&strings[found[7]][0],"%f",&emin);
		sigma = sigma*sigma;
		sigma = sigma*sigma*sigma;
		b[theatom] = sqrt( 4.*emin*sigma*sigma);
		a[theatom] = sqrt( 4.*emin*sigma);
	}
/* a,b values */
	if( found[4] > 0) 
		sscanf(&strings[found[4]][0],"%f",&a[theatom]);
	if( found[5] > 0) 
		sscanf(&strings[found[5]][0],"%f",&b[theatom]);

} /*j*/
/* before we return check for undefined atoms 
	either in the file and not defined in the dictionary
	or     in the dictionary and not in the file 
  neither of these is a lethal
*/
	for(j=0; j< inres; j++)
	{
/*
	printf("%d %d %s\n", defined[j],indic[j],&name[j][0]);
*/
	if( defined[j] ==  0  && indic[j] < 0)
	{ fprintf(stderr,"Serious Warning atom %s in coordinates but not dictionary \n",
	&name[j][0]);
	}
	if( defined[j] <  0  && indic[j] == 0)
	{ fprintf(stderr,"Mild Warning atom %s in dictionary but not coordinates \n",
	&name[j][0]);
	}
	}

} 
