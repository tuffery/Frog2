/*  generate a reflecting sphere
*
*  if the atom is beyond a specified radius
*  AND the radial component of the velocity is
* positive THEN invert the radial component.
*  also make sure the modulus of the velocity is not
* changed
*
*  this acts like a perfect reflecting sphere
*  (of Neutronium??)  and enforces average transfer
*  balance thereby inhibiting evaporation
*
*  sort of a poor man's substitute for 
* periodic boundary conditions
*  but maybe actually better.
*/
/* will depend on a variable
*  bbox which must be set
*  it will do nothing if this is not
*  defined
*/
/*
*  copyright 1998 Robert W. Harrison
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

/* v_box is a NULL function */
int v_box( V,lambda)
float *V,lambda;
{
    return 0;
}


int f_box( lambda)
float lambda;
{
    /* lambda is ignored and the dx == vx in pac et al are
    * also updated */
    float bbox,get_f_variable();
    float xc,yc,zc;
    float rx,ry,rz;
    float rad, vmod,vf,vdot;
    ATOM *ap,*a_next();
    int natom,a_number();
    int i;

    /* check for silly calls */
    natom = a_number();
    if( natom < 1) return 0;
    bbox = get_f_variable("bbox");
    if( bbox < 1.) return 0;
    /* now find the center */

    xc = 0.; yc = 0.; zc = 0.;
    for( i=0;  i< natom; i++)
    {
        ap =  a_next(i);
        xc += ap->x ;
        yc += ap->y ;
        zc += ap->z ;
    }
    xc  /= natom;
    yc  /= natom;
    zc  /= natom;
    /* and now check for inverting the velocities */

    bbox = bbox*bbox;
    for( i=0; i< natom; i++)
    {
        ap = a_next(i);
        rx = ap->x -xc ;
        ry = ap->y -yc ;
        rz = ap->z -zc ;
        rad = rx*rx + ry*ry + rz*rz;
        if( rad > bbox )
        {
            rad = 1./sqrt(rad);
            rx *= rad;
            ry *= rad;
            rz *= rad;
            vmod = ap->vx*ap->vx + ap->vy*ap->vy + ap->vz*ap->vz;
            vdot = ap->vx*rx + ap->vy*ry + ap->vz*rz;
            if( vdot > 0. ){
                ap->vx -= vdot*rx;
                ap->vx -= vdot*rx;
                ap->vy -= vdot*ry;
                ap->vy -= vdot*ry;
                ap->vz -= vdot*rz;
                ap->vz -= vdot*rz;
                vf = ap->vx*ap->vx + ap->vy*ap->vy + ap->vz*ap->vz;
                if( vf > 0.) vmod = sqrt(vmod/vf);

                ap->vx *= vmod;
                ap->vy *= vmod;
                ap->vz *= vmod;
                ap->dx = ap->vx;
                ap->dy = ap->vy;
                ap->dz = ap->vz;

            } /* if vdot > 0. */

        }/* if rad > bbox */
    }/* for( i */

    return 0;
}/* endof f_box */
