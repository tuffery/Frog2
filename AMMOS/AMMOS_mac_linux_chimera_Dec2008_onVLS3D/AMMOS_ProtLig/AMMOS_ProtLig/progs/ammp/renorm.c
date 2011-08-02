/* renorm.c
*
*  find the closest pair with a specified type
*  the pair cannot be in a covalent geometry term (i.e. bond angle torsion hybrid)
*
*  find the average distance for these
*  
*  then set up a restrain-like term for them
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
#include "ammp.h"
/* ATOM structure contains a serial number for indexing into
* arrays and the like (a Hessian)
* but otherwise is self-contained. Note the hooks for Non-restrained potentials
*/
typedef struct{
	ATOM *atom1,*atom2;
	float length,k;
	char active; 
	void *next;
	}  RESTRAIN;
typedef struct{
	char a1[6],a2[6];
	float mean, k;
	void *next;
	RESTRAIN *root; 
} RENORM;

#define RLONG sizeof(RESTRAIN)


static RENORM *renorm_first = NULL;
static RESTRAIN *empty_list = NULL;

//int renorm( char *at1, char *at2, float fk)
int renorm( at1,at2,fk)
char *at1,*at2;
float fk;
{
	RENORM *rp,*newrp;

	newrp = NULL;
	rp = renorm_first;
	while( 1==1)
		{
			if( rp == NULL) break;
		if( strncmp(at1,rp->a1,6) == 0 &&
			strncmp(at2,rp->a2,6) == 0)
			{newrp = rp; break;}
		    rp = rp->next;
		} 

	if( newrp == NULL)
	{
		rp = renorm_first;
		newrp = (RENORM *)malloc( sizeof(RENORM));
		if( newrp == NULL) {aaerror("cannot initialize memory in renorm"); return 0;}
		newrp->next = NULL;
		newrp->mean = -1;
		newrp->root = NULL;
		if( renorm_first == NULL) renorm_first = newrp;
		else
		{while( rp->next != NULL ) rp = rp->next; 
		   rp->next = newrp; }
	}

	strcpy(newrp->a1,at1); strcpy(newrp->a2,at2);
	newrp->k = fk;
	newrp->mean = -1;


}/* renorm definition */

int dump_renorm( op)
FILE *op;
{
	RENORM *rp;
	rp = renorm_first;
	while( 1==1)
	{
		if( rp == NULL) break;
		fprintf(op,"renorm %s %s %f;\n", rp->a1, rp->a2, rp->k);
		rp = rp->next;
	}
}

/* force_renorm_update()
*
* will set all of the mean values to -1;
*
* this will then flag the program to update the
* lists
*
*  it will depend on the mxdq control
*  and get updated when the non-bonded terms are updated
*
*/
void force_renorm_update()
{
	RENORM *rp;
	rp = renorm_first;
	while(1==1) {if( rp== NULL) break;rp->mean = -1.; rp = rp->next; }
}


static int check_renorm_update()
{
	RENORM *rp;
	rp = renorm_first;
	if( rp == NULL) return 0;
	while( 1==1)
	{
		if( rp->mean < 0.) renorm_update(rp);
		rp = rp->next;
		if( rp == NULL ) break;
	}
	return 0;
}

static int renorm_update( RENORM *rp)
{
	RESTRAIN *resp;
	float r,x,y,z,rmin, rminbar;
	ATOM *ap,*bp,*cp, *a_next();
	int numatm, a_number();
	int i,j, inbar;
	int math_match_atom(); /* (who == atomtype, ATOM *ap) returns 0 if not a match*/

	numatm = a_number();
	if( numatm < 2) return 0; /* cannot work if there isn't at least one pair of atoms */

	resp = empty_list;
	if( rp->root == NULL) goto DONT_WALK_THE_LIST;
	if( resp == NULL ){ empty_list = rp->root; resp = rp->root; resp= resp->next;}
	
	while( 1==1)
	{
		if( resp == NULL) break;
		resp = resp->next;
		empty_list->next = resp;
	}

	rp->root = NULL;

DONT_WALK_THE_LIST: ;

	inbar = 0;
	rminbar = 0.;
	ap = a_next(0);
	for( i=0; i< numatm ; i++)
	{
		if(math_match_atom( rp->a1, ap) == 0 ){ap= ap->next; continue;}
	
			if( empty_list == NULL)/* put something on the list so we can take it off */
				{  empty_list = (RESTRAIN *)malloc( sizeof( RESTRAIN));
				if( empty_list == NULL){ aaerror("memory error") ; return;}
				empty_list->next = NULL;
				}
			

		bp = ap->next;
		cp = NULL;
		rmin = 10.e10;
		for( j=i; j < numatm; j++)
		{
			int k,l;
			ATOM *gp;

			/* skip all atoms that aren't candidates */
			if( math_match_atom(rp->a2, bp) == 0){bp = bp->next; continue;}

			
			/* the double loop will catch 1-2,1-3,1-4 interactions */
			for( k=0;k< ap->dontuse; k++)
			{if( ap->excluded[k] == bp) goto SKIP; /* excluded list is symmetric */
					gp = ap->excluded[k];
					for( l=0; l < gp->dontuse; l++)
						if( gp->excluded[l] == bp) goto SKIP; 
			}
			x = ap->x - bp->x;
			y = ap->y - bp->y;
			z = ap->z - bp->z;
			r = x*x + y*y + z*z;
			if( r < rmin){ rmin = r; cp = bp;}
			
SKIP: ;
			bp = bp->next;
		}
		if( cp!= NULL) {rmin = sqrt(rmin); rminbar += rmin; inbar += 1;
				empty_list->atom1 = ap;
				empty_list->atom2 = cp;
				empty_list->length = rmin;
				empty_list->k = rp->k;
				resp = empty_list; empty_list = resp->next;
				resp->next = rp->root;
				rp->root = resp;
				}

		ap = ap->next;
	}

	if(inbar != 0)
	{
		rminbar = rminbar/inbar;
		printf("%f \n", rminbar);
		resp = rp->root;
		x = 0.; y = 0.;
		while( 1==1)
		{
			resp->length = rminbar *0.95;
//		resp->length += rminbar;
//		resp->length *= 0.5;
		resp = resp->next;
		if( resp == NULL) break;
		}

	}
	
}

static RESTRAIN *restrain_first;  /* should overload the value in restrain */

int v_renorm( V, lambda)
float *V, lambda;
{
	RENORM *rp;
	check_renorm_update();

	rp =renorm_first;
	while( 1==1)
	{  if( rp == NULL) break;
		restrain_first = rp->root;
		if( restrain_first != NULL) 
		v_renorm_restrain(V,lambda);
		rp = rp->next;
	}
}

int f_renorm( lambda)
float lambda;
{
	RENORM *rp;
	check_renorm_update();
	rp = renorm_first;
	while( 1==1)
	{  if( rp == NULL) break;
		restrain_first = rp->root;
		if( restrain_first != NULL)
		f_renorm_restrain(lambda);
		rp = rp->next;
	}
}

/* so now we're ready to do the force and potential calculations */
/* basically re-use restrain */

/* v_restrain()
* this function sums up the potentials
* for the atoms defined in the RESTRAIN data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
static int v_renorm_restrain( V, lambda )
	float *V,lambda;
{
	RESTRAIN *bp;
	float r,xt,yt,zt;
	ATOM *a1,*a2;


	bp = restrain_first;
       if( bp == NULL ) return 1;
       while(1)
       {
	if( bp == NULL) return 0;
	a1 = bp->atom1; a2 = bp->atom2;
	if( lambda == 0.)
	{
	r = (a1->x - a2->x)*(a1->x - a2->x);
	r = r + (a1->y - a2->y)*(a1->y - a2->y);
	r = r + (a1->z - a2->z)*(a1->z - a2->z);
	} else
	{
	xt = (a1->x -a2->x +lambda*(a1->dx-a2->dx));
	yt = (a1->y -a2->y +lambda*(a1->dy-a2->dy));
	zt = (a1->z -a2->z +lambda*(a1->dz-a2->dz));
	r = xt*xt+yt*yt+zt*zt;
	}
	r = sqrt(r);
	if( a1->active || a2->active)
	 *V += bp->k*( r - bp->length)*(r - bp->length);

	if( bp == bp->next ) return 1;
	bp = bp->next;
       }
}
/* f_restrain()
*
* f_restrain increments the forces in the atom structures by the force
* due to the restrain components.  NOTE THE WORD increment.
* the forces should first be zero'd.
* if not then this code will be invalid.  THIS IS DELIBERATE.
* on bigger (and better?) machines the different potential terms
* may be updated at random or in parrellel, if we assume that this routine
* will initialize the forces then we can't do this.
*/
static int f_renorm_restrain(lambda)
float lambda;
/*  returns 0 if error, 1 if OK */
{
	RESTRAIN *bp;
	float r,k,ux,uy,uz;
	ATOM *a1,*a2;


	bp = restrain_first;
       if( bp == NULL ) return 1;
       while(1)
       {
	if( bp == NULL) return 0;
	k = bp->k;
	a1 = bp->atom1; a2 = bp->atom2;
	if( lambda == 0.)
	{
	ux = (a2->x - a1->x);
	uy = (a2->y - a1->y);
	uz = (a2->z - a1->z);
	}else{
	ux = (a2->x -a1->x +lambda*(a2->dx-a1->dx));
	uy = (a2->y -a1->y +lambda*(a2->dy-a1->dy));
	uz = (a2->z -a1->z +lambda*(a2->dz-a1->dz));
	}
	r = ux*ux + uy*uy + uz*uz;
	 /* watch for FP errors*/
	 if( r <= 1.e-5)
	 { r = 0; ux = 1.; uy = 0.; uz = 0.; }else{
	r = sqrt(r); ux = ux/r; uy = uy/r; uz = uz/r;
	}
	ux = 2*k*(r-bp->length)*ux; 
	uy = 2*k*(r-bp->length)*uy; 
	uz = 2*k*(r-bp->length)*uz;
	if( a1->active){
	a1->fx += ux; 
	a1->fy += uy; 
	a1->fz += uz; 
		}
	if( a2->active){
	a2->fx -= ux; 
	a2->fy -= uy; 
	a2->fz -= uz; 
	}
	if( bp == bp->next ) return 1;
	bp = bp->next;
       }
}