/* dill.c
*
*   basic form is
*
*  v = k (1./( 1. + exp( (r-d0)/dt) )
*  k < 0 for attractive
*  k > 0 for repulsive
*
*  dill uses dt = 2.5, d0 = 6.5 k = 1.for hydrophobic attractive
*            dt = 0.1  d0 = 3.6 ca, 3.2 cb_centriod. k = 5 for repulsion
*
*
* POOP (Poor-mans Object Oriented Programming) using scope rules
*
* these routines hold a data base (in terms of array indeces)
* of DILLs, with the associated length and force constant
* These are updateable - unlike bonds which are "permanent"
*
* (this could be table driven but what the hell memories cheap)
*
* the routines for potential value, force and (eventually) second
* derivatives are here also
*
* force and 2nd derivative routines assume zero'd arrays for output
* this allows for parralellization if needed (on a PC?)
*
* forces are symmetric 
*/
/*
*  copyright 1992,1994,2000 Robert W. Harrison
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
* but otherwise is self-contained. Note the hooks for Non-bonded potentials
*/
typedef struct{
    ATOM *atom1,*atom2;
    float k,d0,dt;
    void *next;
}  DILL;
#define DILLLONG sizeof(DILL)

DILL *DILL_first = NULL;
DILL *DILL_last = NULL;
/* function DILL adds a DILL to the DILL list
* returns 1 if ok
* returns 0 if not
*  is passed the atom serial numbers, length and constant
* allocates the new memory, initializes it and
* returns
*/
int dill( p1,p2,k,d0,dt)
int p1,p2;
float k,d0,dt ;
{
    ATOM *ap1,*ap2,*a_m_serial();
    DILL *new;
    float x;
    /*	char line[80]; */
    /* get the atom pointers for the two serial numbers */
    ap1 = a_m_serial( p1 );
    ap2 = a_m_serial( p2 );
    if( (ap1 == NULL) || (ap2 == NULL) )
    {
        /*      sprintf( line,"undefined atom in DILL %d %d \0",p1,p2);
        	aaerror( line );
        */
        return 1;
    }
    /* check to see if a DILLt is already defined */
#ifdef replace_DILL
    new = DILL_first;
    if( new != NULL)
    {
        while(1)
        {
            if( new == NULL) break;
            if( (new->atom1 == ap1 && new->atom2 == ap2) ||
                    (new->atom1 == ap2 && new->atom2 == ap1) )
            {
                new->dt = dt;
                new->d0 = d0;
                new->k = k;
                return 1;
            }
            if( new == new->next) break;
            new = new->next;
        }
    }
#endif
    if( ( new = malloc( DILLLONG ) ) == NULL)
    {
        return 0;
    }
    /* initialize the pointers */
    if( DILL_first == NULL) DILL_first = new;
    if( DILL_last == NULL) DILL_last = new;
    new -> atom1 = ap1;
    new -> atom2 = ap2;
    new->dt = dt;
    new->d0 = d0;
    new->k = k;
    new -> next = new;
    DILL_last -> next = new;
    DILL_last = new;
    return 1;
}




/* v_dill()
* this function sums up the potentials
* for the atoms defined in the DILL data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int v_dill( V, lambda )
float *V,lambda;
{
    DILL *bp;
    float r,xt,yt,zt;
    ATOM *a1,*a2;


    bp = DILL_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1; a2 = bp->atom2;
        if( a1->active || a2->active ){
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
            r = (r-bp->d0)/bp->dt;
            if( r <10.)
                *V += bp->k*(1./(1.+ exp( r)));
        }
        if( bp == bp->next ) return 1;
        bp = bp->next;

    }
}
/* f_dill()
*
* f_dill increments the forces in the atom structures by the force
* due to the DILL components.  NOTE THE WORD increment.
* the forces should first be zero'd.
* if not then this code will be invalid.  THIS IS DELIBERATE.
* on bigger (and better?) machines the different potential terms
* may be updated at random or in parrellel, if we assume that this routine
* will initialize the forces then we can't do this.
*/
int f_dill(lambda)
float lambda;
/*  returns 0 if error, 1 if OK */
{
    DILL *bp;
    float r,t,ux,uy,uz;
    ATOM *a1,*a2;


    bp = DILL_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1; a2 = bp->atom2;
        if( a1->active || a2->active){
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
            r = ux*ux + uy*uy + uz*uz ;
            /* watch for FP errors*/
            if( r <= 1.e-5)
            { r = 0; ux = 1.; uy = 0.; uz = 0.;}else{
                r = sqrt(r); t = 1/r; ux = ux*t; uy = uy*t; uz = uz*t;
            }
            /*		*V += bp->k*(1./(1.+ exp( (r-bp->d0)/bp->dt));
            */
            r = (r-bp->d0)/bp->dt;
            if( r < 10.)
            {
                t = exp( r);
                r = -t/((1.+t)*(1.+t)*bp->dt);
            } else {
                r = 0.;
            }
            ux *= r;
            uy *= r;
            uz *= r;
            if( a1->active){
                a1->fx += ux;
                a1->fy += uy;
                a1->fz += uz;
            }

            if( a2->active) {
                a2->fx -= ux;
                a2->fy -= uy;
                a2->fz -= uz;
            }
        }
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* function get_dill( a1,DILLed,10,inDILL);
* check the DILLS list for atoms DILLed to a1
*/
void get_dill( a1,DILLed,mDILL,inDILL)
ATOM *a1, *DILLed[];
int mDILL,*inDILL ;
{
    DILL *mine;
    mine = DILL_first;
    *inDILL = 0;
    while(1)
    {
        if( (mine == NULL) )
        {
            return;
        }
        if( mine->atom1 == a1)
        {
            DILLed[(*inDILL)++] = mine->atom2;
        }
        if( mine->atom2 == a1)
        {
            DILLed[(*inDILL)++] = mine->atom1;
        }
        if( mine == mine->next) return;
        mine = mine->next;
        if( *inDILL == mDILL ) return;
    }
}
/* function get_dill( a1,DILLed,10,inDILL);
* check the DILLS list for atoms DILLed to a1
*/
void get_dill_and_length( a1,DILLed,r,mDILL,inDILL)
ATOM *a1, *DILLed[];
int mDILL,*inDILL ;
float r[];
{
    DILL *mine;
    mine = DILL_first;
    *inDILL = 0;
    while(1)
    {
        if( (mine == NULL) )
        {
            return;
        }
        if( mine->atom1 == a1)
        {
            r[*inDILL] = mine->d0;
            DILLed[(*inDILL)++] = mine->atom2;
        }
        if( mine->atom2 == a1)
        {
            r[*inDILL] = mine->d0;
            DILLed[(*inDILL)++] = mine->atom1;
        }
        if( mine == mine->next) return;
        mine = mine->next;
        if( *inDILL == mDILL ) return;
    }
}
/* DILL_next()
*  like bond_next() but for DILL structures
*/
int DILL_next( int i, ATOM **n1, ATOM **n2  )
{
    static DILL *np = NULL ;
    *n1 = NULL ; *n2 = NULL;
    if( DILL_first == NULL ) return (1==0);
if( np == NULL || i <= 0 ) { np = DILL_first;}
    else{ np = np->next; }
    *n1 = np->atom1; *n2 = np->atom2;
    if( np->next != np)return (1==1);
    return ( 1== 0);
}
/* routine dump_dills
* this function outputs the DILL parameters
* and does it in a simple form
* DILL ser1,ser2,k,req
* the rest is just free format
*/
void dump_dills( where )
FILE *where;
{
    DILL *b;
    ATOM *a1,*a2;
    b = DILL_first;
    if( b == NULL ) return;
    while( (b->next != b) )
    {
        if( b->next == NULL) return;
        a1 = b->atom1; a2 = b->atom2;
        fprintf( where,"dill %d %d %f %f %f ;\n",a1->serial,a2->serial,
                 b->k,b->d0,b->dt);
        b = b->next;
    }
    if( b->next == NULL) return;
    a1 = b->atom1; a2 = b->atom2;
    fprintf( where,"dill %d %d %f %f %f ;\n",a1->serial,a2->serial,
             b->k,b->d0,b->dt);
}

/* a_dill()
* this function sums up the potentials
* for the atoms defined in the DILL data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int a_dill( V, lambda,ilow,ihigh,op )
float *V,lambda;
int ilow,ihigh;
FILE *op;
{
    DILL *bp;
    float r,xt,yt,zt;
    ATOM *a1,*a2;


    bp = DILL_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1; a2 = bp->atom2;
        if(( a1->serial >= ilow && a1->serial <=ihigh)
                ||( a2->serial >= ilow && a2->serial <=ihigh))
        {
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
            zt = 0;
            zt =  bp->k*(1./(1.+exp( (r-bp->d0)/bp->dt) ));
            *V += zt;
            fprintf(op,"dill %s %d %s %d E %f value %f\n"
                    ,a1->name,a1->serial,a2->name,a2->serial,zt,r);
        }
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* gsdg_dill( ATOM *ap )
*  
* setup the distances for DILL terms
*/
int gsdg_dill( ap)
ATOM *ap;
{
    ATOM *bp;
    DILL *np;

    np = DILL_first;
    while(1)
    { if( np == NULL ) return 0;
        if( np->atom1 == ap )
        {  bp = np->atom2; bp->vx = (np->d0*np->d0 );
            bp->vy = np->k; }
        if( np->atom2 == ap )
        {  bp = np->atom1; bp->vx = (np->d0*np->d0 );
            bp->vy = np->k; }

        if( np == np->next ) return 0;
        np = np->next;
    }
}
