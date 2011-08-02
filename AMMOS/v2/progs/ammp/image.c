/* image charge solvent model
*
*  1) find the center of charge and total charge
*  2) apply a -(total_charge) field form the center of charge
*  
*  use a guassian model for the non-sharp boundary
*/

/*
*  copyright 1997 Robert W. Harrison
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


/*
#define ZETA 5.
#define ZETA2 25.
#define ZETA 1.
#define ZETA2 1.
*/
#define ZETA 2.
#define ZETA2 4.


/* zeta and zeta2 are the zeta layer depth and the square of it. */

/* these are defined in orbit.c */
/*  the next lines can be uncommented when prototypes are used */
/*float Fzero(float) ;  special integral for integral( exp(-x*x)/r)
* float Fone(float) ;  -d/dx(Fzero) 
*/
float image_get_center( float *,float *, float *,float);

int v_image( V, lambda)
float *V,lambda;
{
    ATOM *ap,*a_next();
    int i,na,a_number();
    float qt,x,y,z,r,t,vt;
    float Fzero();

    qt = image_get_center(&x,&y,&z,lambda);
    if( qt < 1.e-7 && qt > -1.e-7) return;
    na = a_number();
    vt = 0.;
    for( i=0; i< na; i++)
    {
        ap = a_next(i);
        t = x - (ap->x +lambda*ap->dx);
        r = t*t;
        t = y - (ap->y +lambda*ap->dy);
        r += t*t;
        t = z - (ap->z +lambda*ap->dz);
        r += t*t;
        vt += ap->q*Fzero( r*ZETA2);

    }
    vt = vt * qt* 332.17752* TWOPI/ZETA2	;
    *V -= vt;
}/* end of v_image */


int f_image( lambda)
float lambda;
{
    ATOM *ap,*a_next();
    int i,na,a_number();
    float qt,x,y,z,r,t,dr;
    float Fone();
    qt = image_get_center(&x,&y,&z,lambda);
    if( qt < 1.e-7 && qt > -1.e-7) return;
    na = a_number();
    for( i=0; i< na; i++)
    {
        ap = a_next(i);
        t = x - (ap->x +lambda*ap->dx);
        r = t*t;
        t = y - (ap->y +lambda*ap->dy);
        r += t*t;
        t = z - (ap->z +lambda*ap->dz);
        r += t*t;
        /*		*V += 332.17752*ap->q*qt*TWOPI/ZETA2*Fzero( ZETA2*r);
        */
        t =  332.17752*ap->q*qt*TWOPI/ZETA2*Fone( ZETA2*r);
        r = sqrt(r);
        if( r > 1.e-5) {
            dr = (x - (ap->x +lambda*ap->dx))/r;
            ap->fx += t*dr;
            dr = (y - (ap->y +lambda*ap->dy))/r;
            ap->fy += t*dr;
            dr = (z - (ap->z +lambda*ap->dz))/r;
            ap->fz += t*dr;

        }/* r> 1.e-5 */
    }
}/* end of f_image */

/* get the center of charge and the total charge */
/*
float image_get_center( x,y,z,lambda)
float *x,*y,*z,lambda;
*/
float image_get_center(float *x, float *y, float *z,float lambda)
{
    ATOM *ap,*a_next();
    int na,a_number();
    int i;
    float qt;

    na = a_number();
    if( na <= 1) return 0.;


    qt = 0.;
    *x = 0.;
    *y = 0.;
    *z = 0.;
    for( i=0; i< na; i++)
    {
        ap = a_next(i);
        qt +=  ap->q;
        *x += (ap->x+lambda*ap->dx)*ap->q;
        *y += (ap->y+lambda*ap->dy)*ap->q;
        *z += (ap->z+lambda*ap->dz)*ap->q;
    }
    *x /= na;
    *y /= na;
    *z /= na;
    return qt;
}

int a_image( V,lambda,ilow,ihigh,op)
float *V,lambda;
int ilow,ihigh;
FILE *op;
{
    ATOM *ap,*a_next();
    int i,na,a_number();
    float qt,x,y,z,r,t,vt,vtt;
    float Fzero();

    qt = image_get_center(&x,&y,&z,lambda);
    na = a_number();
    if( na == 0) return;
    for( i=0; i< na; i++)
    {
        ap = a_next(i);
        if( ap->serial <= ihigh && ap->serial>=ilow)
        {
            t = x - (ap->x +lambda*ap->dx);
            r = t*t;
            t = y - (ap->y +lambda*ap->dy);
            r += t*t;
            t = z - (ap->z +lambda*ap->dz);
            r += t*t;
            vt = 332.17752*ap->q*qt*TWOPI/ZETA2*Fzero( ZETA2*r);
            vtt += vt;
            fprintf(op,"Image Charge energy %s %d %f\n",
                    ap->name,ap->serial,vt);
        } /* if( ap->serial */
    }/* for i */
    fprintf(op,"Total image charge energy %f\n",vtt);
    *V += vtt;
}
