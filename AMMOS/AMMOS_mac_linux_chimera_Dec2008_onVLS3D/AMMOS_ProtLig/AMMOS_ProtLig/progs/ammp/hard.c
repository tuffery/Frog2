/* Hard.c
*
*  non-bonded potentials with 
*   very short fixed cutoff
*   only repulsive terms (although using the attractive ones) 
*
*   primarily for generating big structures 
*/
/*
*  copyright 1992,1993,1994 Robert W. Harrison
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



v_hard( V,lambda)
float *V,lambda;
{
    float x,y,z;
    float r, r2,r6,r12;
    float r13,r7,dr;
    ATOM *ap,*bp,*a_next();
    int i, j,k,numatm,a_number();
    float dielec,get_f_variable();

#define CUTOFF 4.
#define LETOFF -4.

    dielec = get_f_variable("dielec");
    if( dielec < 1.) dielec = 1.;
    dielec = 332.17752/dielec;

    numatm = a_number();

    /* make ap point to the second atom */
    ap = a_next(-1);
    for( i=1; i< numatm; i++)
    {  ap = ap->next;  if( ap == NULL ) break;
        for( j=0; j< i; j++)
        { bp = a_next(j);
            /* check for exclusions (use getbond ... for dos version */
            for( k=0; k< ap->dontuse; k++)
            { if( bp == ap->excluded[k]) goto SKIP ; }
            /* if here lets go */
            if( lambda > 0.) {
                x = ap->x - bp->x + lambda*(ap->dx -bp->dx);
                if( x > CUTOFF || x <  LETOFF) goto SKIP;
                y = ap->y - bp->y + lambda*(ap->dy -bp->dy);
                if( y > CUTOFF || y <  LETOFF) goto SKIP;
                z = ap->z - bp->z + lambda*(ap->dz -bp->dz);
                if( z > CUTOFF || z <  LETOFF) goto SKIP;
            } else {
                x = ap->x - bp->x ;
                if( x > CUTOFF || x <  LETOFF) goto SKIP;
                y = ap->y - bp->y ;
                if( y > CUTOFF || y <  LETOFF) goto SKIP;
                z = ap->z - bp->z ;
                if( z > CUTOFF || z <  LETOFF) goto SKIP;
            }

            r2 = x*x + y*y + z*z;
            if( r2 < 1. ) r2 = 1.;
            r6 = r2*r2*r2; r12 = r6*r6;
            r = sqrt(r2);
            r7 = r6*r; r13 = r12*r;
            dr = -dielec*ap->q*bp->q/r2 + 6*ap->a*bp->a/r7 - 12*ap->b*bp->b/r13;
            /*
            	printf(" hard> %d %d %f %f\n", ap->serial,bp->serial,r,dr);
            */
            if( dr < 0.)
                *V += dielec*ap->q*bp->q/r - ap->a*bp->a/r6 + ap->b*bp->b/r12;


SKIP:    k = j;
        }

        if( ap == ap->next ) break;
    }

}/* end of routine */



f_hard( lambda)
float lambda;
{
    float x,y,z;
    float ux,uy,uz;
    float r, r2,r6,r12;
    float r13,r7,dr;
    ATOM *ap,*bp,*a_next();
    int i, j,k,numatm,a_number();
    float dielec,get_f_variable();

#define CUTOFF 4.
#define LETOFF -4.

    dielec = get_f_variable("dielec");
    if( dielec < 1.) dielec = 1.;
    dielec = 332.17752/dielec;

    numatm = a_number();

    /* make ap point to the second atom */
    ap = a_next(-1);
    for( i=1; i< numatm; i++)
    {  ap = ap->next;  if( ap == NULL ) break;
        for( j=0; j< i; j++)
        { bp = a_next(j);
            /* check for exclusions (use getbond ... for dos version */
            for( k=0; k< ap->dontuse; k++)
            { if( bp == ap->excluded[k]) goto SKIP ; }
            /* if here lets go */
            if( lambda > 0.) {
                x = ap->x - bp->x + lambda*(ap->dx -bp->dx);
                if( x > CUTOFF || x <  LETOFF) goto SKIP;
                y = ap->y - bp->y + lambda*(ap->dy -bp->dy);
                if( y > CUTOFF || y <  LETOFF) goto SKIP;
                z = ap->z - bp->z + lambda*(ap->dz -bp->dz);
                if( z > CUTOFF || z <  LETOFF) goto SKIP;
            } else {
                x = ap->x - bp->x ;
                if( x > CUTOFF || x <  LETOFF) goto SKIP;
                y = ap->y - bp->y ;
                if( y > CUTOFF || y <  LETOFF) goto SKIP;
                z = ap->z - bp->z ;
                if( z > CUTOFF || z <  LETOFF) goto SKIP;
            }

            r2 = x*x + y*y + z*z;
            if( r2 < 1.e-5) goto SKIP;
            r = sqrt(r2);
            ux = x/r; uy = y/r; uz = z/r;
            if( r2 < 1. ) r2 = 1.;
            r6 = r2*r2*r2; r12 = r6*r6;
            r = sqrt(r2);
            r7 = r6*r; r13 = r12*r;
            dr = -dielec*ap->q*bp->q/r2 + 6*ap->a*bp->a/r7 - 12*ap->b*bp->b/r13;
            if( dr < 0.)
            {
                ux *= dr;
                uy *= dr;
                uz *= dr;
                ap->fx -= ux;
                ap->fy -= uy;
                ap->fz -= uz;
                bp->fx += ux;
                bp->fy += uy;
                bp->fz += uz;
            }


SKIP:    k = j;
        }

        if( ap == ap->next ) break;
    }

}/* end of routine */

