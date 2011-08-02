/* 
* ghost.c
*
*  ghost allows the easy definition of sets of mutually invisible atoms
*
*  for estimating binding energies.  Basically put waters where the ligand is
*  and then analyze each ghost bin
*
*
*  ghost  ghostid  low_atom high_atom;  
*  buster  ghostid;
*  dump ghost;
*
*/

/*
*  copyright 1992-2001 Robert W. Harrison
*  
*  This notice may not be removed
*  This program may be copied for scientific use
*  It may not be sold for profit without explicit
*  permission of the author(s) who retain any
*  commercial rights including the right to modify 
*  this notice
*/

#define ANSI 1
/* misc includes - ANSI and some are just to be safe */
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#ifdef ANSI
#include <stdlib.h>
#endif
#include "ammp.h"
/* ATOM structure contains a serial number for indexing into
* arrays and the like (a Hessian)
* but otherwise is self-contained. Note the hooks for Non-debyeed potentials
*/



typedef struct {
	int low,high; /* serial numbers */
	int serial ;  /* ghost bin */
	void *next_ghost;
}GHOST;

#define GLONG sizeof(GHOST)

static GHOST *firstGHOST = NULL;

int aaerror(char *);
/* the analog of the constructor */

int ghost( int serial, int low, int high)
{
	int tailor_exclude(int,int); /* defined in tailor.c */

	GHOST *gp;

	gp = (GHOST *)malloc( GLONG);
	if( gp == NULL) {aaerror("memory allocation error in ghost\n"); exit(0);}

	gp->next_ghost = firstGHOST;
	firstGHOST = gp;  

	if( low > high) { int i = low; low = high; high = i;}
	gp->low = low;
	gp->high = high;
	gp->serial = serial;

	while( gp->next_ghost != NULL)
	{
		int i,j;
		gp = gp->next_ghost;
		if( gp->serial == serial) continue; /* this allows multiple regions */
		for( i= low; i<= high; i++)
			for( j= gp->low; j <= gp->high; j++)
				tailor_exclude(i,j);
	}
	return 0;
}/* end of ghost */


void buster( int (*vfs[])(float *, float), int nfs, int serial, FILE *op)
{
	GHOST *gp;
	int analyze( int (*[])(float* ,float), int, int, int, FILE *);

	gp = firstGHOST;

	while( gp != NULL)
	{
		if( gp->serial == serial)
		{analyze( vfs,nfs, gp->low, gp->high, op); }
		gp = gp->next_ghost;
	}
	
}/* end of buster */


int dump_ghost( FILE *op)
{
	GHOST *gp;
	gp = firstGHOST;

	while(gp != NULL)
	{
		fprintf(op,"ghost %d %d %d ;\n", gp->serial, gp->low, gp->high);
		gp = gp->next_ghost;
	}
	return 0;
}/* end of dump_ghost*/
