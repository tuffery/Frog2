/* restrain.c
*
* collection of routines to service restrain potentials
*
* POOP (Poor-mans Object Oriented Programming) using scope rules
*
* these routines hold a data base (in terms of array indeces)
* of restraints, with the associated length and force constant
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
*  copyright 1992 Robert W. Harrison
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
* but otherwise is self-contained. Note the hooks for Non-restrained potentials
*/
typedef struct{
    ATOM *atom1,*atom2;
    float length,k;
    void *next;
}  RESTRAIN;
#define RLONG sizeof(RESTRAIN)

RESTRAIN *restrain_first = NULL;
RESTRAIN *restrain_last = NULL;
/* function restrain adds a restrain to the restrain list
* returns 1 if ok
* returns 0 if not
*  is passed the atom serial numbers, length and constant
* allocates the new memory, initializes it and
* returns
*/
int restrain( p1,p2,bl,fk)
int p1,p2;
float bl,fk ;
{
    ATOM *ap1,*ap2,*a_m_serial();
    RESTRAIN *new;
    char line[BUFSIZ];
    /* get the atom pointers for the two serial numbers */
    ap1 = a_m_serial( p1 );
    ap2 = a_m_serial( p2 );
    if( (ap1 == NULL) || (ap2 == NULL) )
    {
        sprintf( line,"undefined atom in restrain %d %d \0",p1,p2);
        aaerror( line );
        return 0;
    }
    /* check to see if a restraint is already defined */
    new = restrain_first;
    if( new != NULL)
    {
        while(1)
        {
            if( new == NULL) break;
            if( (new->atom1 == ap1 && new->atom2 == ap2) ||
                    (new->atom1 == ap2 && new->atom2 == ap1) )
            {
                new->length = bl; new->k = fk; return 1;
            }
            if( new == new->next) break;
            new = new->next;
        }
    }
    if( ( new = malloc( RLONG ) ) == NULL)
    {
        return 0;
    }
    /* initialize the pointers */
    if( restrain_first == NULL) restrain_first = new;
    if( restrain_last == NULL) restrain_last = new;
    new -> atom1 = ap1;
    new -> atom2 = ap2;
    new -> length = bl;
    new -> k = fk;
    new -> next = new;
    restrain_last -> next = new;
    restrain_last = new;
    return 1;
}


/* v_restrain()
* this function sums up the potentials
* for the atoms defined in the RESTRAIN data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int v_restrain( V, lambda )
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
int f_restrain(lambda)
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
/* function get_restrain( a1,restrained,10,inrestrain);
* check the RESTRAINS list for atoms restrained to a1
*/
void get_restrain( a1,restrained,mrestrain,inrestrain)
ATOM *a1, *restrained[];
int mrestrain,*inrestrain ;
{
    RESTRAIN *mine;
    mine = restrain_first;
    *inrestrain = 0;
    while(1)
    {
        if( (mine == NULL) )
        {
            return;
        }
        if( mine->atom1 == a1)
        {
            restrained[(*inrestrain)++] = mine->atom2;
        }
        if( mine->atom2 == a1)
        {
            restrained[(*inrestrain)++] = mine->atom1;
        }
        if( mine == mine->next) return;
        mine = mine->next;
        if( *inrestrain == mrestrain ) return;
    }
}
/* routine dump_restrains
* this function outputs the restrain parameters
* and does it in a simple form
* restrain ser1,ser2,k,req
* the rest is just free format
*/
void dump_restrains( where )
FILE *where;
{
    RESTRAIN *b;
    ATOM *a1,*a2;
    b = restrain_first;
    if( b == NULL ) return;
    while( (b->next != b) )
    {
        if( b->next == NULL) return;
        a1 = b->atom1; a2 = b->atom2;
        fprintf( where,"restrain %d %d %f %f \;\n",a1->serial,a2->serial,
                 b->length,b->k);
        b = b->next;
    }
    if( b->next == NULL) return;
    a1 = b->atom1; a2 = b->atom2;
    fprintf( where,"restrain %d %d %f %f \;\n",a1->serial,a2->serial,
             b->length,b->k);
}

/* a_restrain()
* this function sums up the potentials
* for the atoms defined in the RESTRAIN data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int a_restrain( V, lambda,ilow,ihigh,op )
float *V,lambda;
int ilow,ihigh;
FILE *op;
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
            r = sqrt(r); zt= bp->k*( r - bp->length)*(r - bp->length);
            *V += zt;
            fprintf(op,"Restrain %d %d E %f value %f error %f\n"
                    ,a1->serial,a2->serial,zt,r,r-bp->length);
        }
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
