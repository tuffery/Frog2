
/* central.c
*
* collection of routines to service central mean force  potentials
*
* POOP (Poor-mans Object Oriented Programming) using scope rules
*
* these routines hold a data base (in terms of array indeces)
* of central bonds, with the associated length and force constants
*
* (this could be table driven but what the hell memories cheap)
*
* the routines for potential value, force and (eventually) second
* derivatives are here also
*
* force and 2nd derivative routines assume zero'd arrays for output
* this allows for parralellization if needed (on a PC?)
*
* forces are bond wise symmetric - so we don't have to fuck around with
* s matrices and the like.
*/
/*
*  copyright 1999 Robert W. Harrison
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
* but otherwise is self-contained. Note the hooks for Non-centraled potentials
*
*/
typedef struct{
    ATOM *atom;
    float weight,exp;
    void *next;
}  CENTRAL;
#define CLONG sizeof(CENTRAL)

CENTRAL *central_first = NULL;
CENTRAL *central_last = NULL;
/* function central adds a central to the central list
* returns 1 if ok
* returns 0 if not
*  is passed the atom serial numbers, length and constant
* allocates the new memory, initializes it and
* returns
*/
int central( p1, weight,exp)
int p1;
float weight,exp;
{
    ATOM *ap1,*a_m_serial();
    CENTRAL *new;
    char line[BUFSIZ];
    /* get the atom pointers for the two serial numbers */
    ap1 = a_m_serial( p1 );
    if( (ap1 == NULL)  )
    {
        sprintf( line,"undefined atom in central %d \0",p1);
        aaerror( line );
        return 0;
    }
    new = central_first;
    if( central_first != NULL)
    {
        while( 1 )
        {
            if( new->atom == ap1)
            {
                new -> weight = weight;
                new -> exp = exp;
                return 1;
            }
            if( new->next == new) break;
            if( new->next == NULL) break;
            new = new->next;
        }
    }

    if( ( new = malloc( CLONG ) ) == NULL) return 0;
    new -> next = new;
    if( central_first == NULL) central_first = new;
    if( central_last == NULL) central_last = new;
    central_last -> next = new;
    central_last = new;
    /* initialize the pointers */
    new -> atom = ap1;
    new -> weight = weight;
    new -> exp = exp;
    return 1;
}


/* v_central()
* this function sums up the potentials
* for the atoms defined in the CENTRAL data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int v_central( V, lambda )
float *V,lambda;
{
    CENTRAL *bp;
    float r,xt,yt,zt;
    float cx,cy,cz;
    int numatom,a_number(),i;
    ATOM *a1,*a_next();;

    numatom = a_number();
    if( numatom <= 2) return 0;

    cx = 0.; cy = 0.; cz = 0.;
    for( i=0; i< numatom; i++)
    {
        a1 = a_next(i);
        cx += a1->x + lambda*a1->dx;
        cy += a1->y + lambda*a1->dy;
        cz += a1->z + lambda*a1->dz;
    }
    cx /= numatom;
    cy /= numatom;
    cz /= numatom;

    bp = central_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom;
        if( a1->active ){
            xt = a1->x + lambda*a1->dx -cx;
            yt = a1->y + lambda*a1->dy -cy;
            zt = a1->z + lambda*a1->dz -cz;
            xt = sqrt(xt*xt + yt*yt + zt*zt);
            *V += bp->weight*exp(bp->exp*xt);
        }
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* f_central()
*
* f_central increments the forces in the atom structures by the force
* due to the central components.  NOTE THE WORD increment.
* the forces should first be zero'd.
* if not then this code will be invalid.  THIS IS DELIBERATE.
* on bigger (and better?) machines the different potential terms
* may be updated at random or in parrellel, if we assume that this routine
* will initialize the forces then we can't do this.
*/
int f_central(lambda)
float lambda;
/*  returns 0 if error, 1 if OK */
{
    CENTRAL *bp;
    float ux,uy,uz;
    float r,xt,yt,zt;
    float cx,cy,cz;
    float (*cdx)[],(*cdy)[],(*cdz)[]; /* needed for center corrections */
    int numatom,a_number(),i;
    ATOM *a1,*a_next();
    ATOM *(*aps)[];

    numatom = a_number();
    if( numatom <= 2) return 0;
    bp = central_first;
if( bp == NULL ){ return 1;}

    cdx = malloc( numatom*sizeof(float));
    cdy = malloc( numatom*sizeof(float));
    cdz = malloc( numatom*sizeof(float));
    aps = malloc( numatom*sizeof(ATOM*));
    if( cdx == NULL || cdy == NULL || cdz == NULL || aps == NULL)
    {aaerror("cannot allocate memory in f_central"); return 0;}

    cx = 0.; cy = 0.; cz = 0.;
    for( i=0; i< numatom; i++)
    {
        a1 = a_next(i);
        (*aps)[i] = a1;
        cx += a1->x + lambda*a1->dx;
        cy += a1->y + lambda*a1->dy;
        cz += a1->z + lambda*a1->dz;
    }
    cx /= numatom;
    cy /= numatom;
    cz /= numatom;
    /* calculate the unit vectors to the center */
    for( i=0; i< numatom; i++)
    {
        a1 = (*aps)[i];
        xt = a1->x + lambda*a1->dx - cx;
        yt = a1->y + lambda*a1->dy - cy;
        zt = a1->z + lambda*a1->dz - cz;
        r = sqrt( xt*xt + yt*yt + zt*zt);
        if( r > 1.e-7) r = 1./r;
        (*cdx)[i] = xt*r;
        (*cdy)[i] = yt*r;
        (*cdz)[i] = zt*r;
    }

    while(1)
    {
        if( bp == NULL) return 0;
    if( bp == NULL ){free(aps); free(cdz); free(cdy); free(cdx); return 1;}
        a1 = bp->atom;
        if( a1->active ){
            xt = a1->x + lambda*a1->dx -cx;
            yt = a1->y + lambda*a1->dy -cy;
            zt = a1->z + lambda*a1->dz -cz;
            xt = sqrt(xt*xt + yt*yt + zt*zt);
            /* *V += bp->weight*exp(bp->exp*xt); */
            /* f = -dV/dx */
            yt = -bp->weight*exp(bp->exp*xt)*bp->exp;
            zt = yt/(numatom-1);
            for( i=0; i< numatom; i++)
            {
                a1 = (*aps)[i];
                if( a1->active){
                    if( a1 == bp->atom)
                    {
                        a1->fx += (*cdx)[i]*yt;
                        a1->fy += (*cdy)[i]*yt;
                        a1->fz += (*cdz)[i]*yt;
                    }else{
                        a1->fx -= (*cdx)[i]*zt;
                        a1->fy -= (*cdy)[i]*zt;
                        a1->fz -= (*cdz)[i]*zt;
                    }
                }}/* for active and i */
        }
        if( bp == bp->next ){free(aps); free(cdz); free(cdy); free(cdx); return 1;}
        bp = bp->next;
    }
    /* never get here */
}
/* routine dump_centrals
* this function outputs the central parameters
* and does it in a simple form
* central ser1,ser2,k,req
* the rest is just free format
*/
void dump_central( where )
FILE *where;
{
    CENTRAL *b;
    ATOM *a1;
    b = central_first;
    if( b == NULL ) return;
    while( (b->next != b) )
    {
        if( b->next == NULL) return;
        a1 = b->atom;
        fprintf( where,"central %d %f %f  ;\n",a1->serial,
                 b->weight,b->exp);
        b = b->next;
    }
    if( b->next == NULL) return;
    a1 = b->atom;
    fprintf( where,"central %d %f %f  ;\n",a1->serial,
             b->weight,b->exp);
}
