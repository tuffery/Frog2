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
/* will depend on  variables
*  acube, bcube, ccube  which must be set
*  this will force the coordinates to lie within a 
*  orthorhombic box
* these are the half sizes of the box
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

/* v_cube is a NULL function */
int v_cube( V,lambda)
float *V,lambda;
{
    return 0;
}


int f_cube( lambda)
float lambda;
{
    /* lambda is ignored and the dx == vx in pac et al are
    * also updated */
    float acube,bcube,ccube,get_f_variable();
    float xc,yc,zc;
    float rx,ry,rz;
    float rad, vmod,vf,vdot;
    ATOM *ap,*a_next();
    int natom,a_number();
    int i;

    /* check for silly calls */
    natom = a_number();
    if( natom < 1) return 0;
    acube = get_f_variable("acube");
    if( acube < 1.) return 0;
    bcube = get_f_variable("bcube");
    if( bcube < 1.) return 0;
    ccube = get_f_variable("ccube");
    if( ccube < 1.) return 0;
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

    for( i=0; i< natom; i++)
    {
        ap = a_next(i);
        rx = (ap->x -xc) ;
        ry = (ap->y -yc) ;
        rz = (ap->z -zc) ;
        if( fabs(rx) > acube )
        {
            if( (ap->vx > 0 && rx > 0) || (ap->vx<0 &&rx<0) ){
                ap->vx = -ap->vx;
                ap->dx = -ap->dx;
            }
        }/* if rx > acube */
        if( fabs(ry) > bcube )
        {
            if( (ap->vy > 0 && ry > 0) || (ap->vy<0 &&ry<0) ){
                ap->vy = -ap->vy;
                ap->dy = -ap->dy;
            }
        }/* if ry > bcube */
        if( fabs(rz) > ccube )
        {
            if( (ap->vz > 0 && rz > 0) || (ap->vz<0 &&rz<0) ){
                ap->vz = -ap->vz;
                ap->dz = -ap->dz;
            }
        }/* if rz > ccube */
    }/* for( i */

    return 0;
}/* endof f_cube */
