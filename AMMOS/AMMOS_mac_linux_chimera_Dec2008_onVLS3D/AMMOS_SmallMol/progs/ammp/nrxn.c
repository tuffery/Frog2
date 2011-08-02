/* rxnfield.c
*
* expand the distal feild in a multipole expansion
* then use the conformal map a/r to generate the reaction feild
* with the (more or less) correct result for the reflection at a
* boundary.
*
*
* code taken from rectmm.c
*
*
*  rectangular multipole expansion to the r^-6 order (5th order expansion)
*
*  from Eyges "The Classical Electromagnetic Field"
*  note that we use the oposite convention for sign of
*  expansion so use + for all the cumulants, while he uses
*   -1^n.   This is solely due to choice of origin and for
*   other applications (rxn feild ) the -1^n is correct.
*
*
* collection of routines to service nonbonded potentials
*
* POOP (Poor-mans Object Oriented Programming) using scope rules
*
* the routines for potential value, force and (eventually) second
* derivatives are here also
*
* force and 2nd derivative routines assume zero'd arrays for output
* this allows for parralellization if needed (on a PC?)
*
* forces are symmetric - so we don't have to fuck around with
* s matrices and the like.
*
* note that the non-bonded information is in the ATOM structures 
*
*
* attempts at vectorization
*/
/*
*  copyright 1992, 1993, 1994, 1995 Robert W. Harrison
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
* but otherwise is self-contained. Note the hooks for Non-nonboned potentials
*/

/*
#define FOURTH
#define  FIFTH
#ifdef FIFTH
#define FOURTH
#endif
*/
typedef struct {
    float xc,yc,zc;
    float q000;
    float q100,q010,q001;
    float q200,q020,q002,q110,q101,q011;
    float q300,q030,q003,q210,q201,q120,q021,q102,q012,q111;
#ifdef FOURTH
    float q400,q040,q004,q310,q301,q130,q031,q103,q013,q220,q202,q022,q211,q121,q112;
#endif
#ifdef FIFTH
    float q500,q050,q005,q410,q401,q140,q041,q104,q014,q320,q230,q302,q203,q032,q023,q311,q131,q113,q221,q212,q122;
#endif
} MULTIPOLE;
int rxn_common_routine( mmp ,lambda)
MULTIPOLE *mmp;
float lambda;
{
    float r,r0,xt,yt,zt;
    float xc,yc,zc;
    float xt2,xt3,xt4,yt2,yt3,yt4,zt2,zt3,zt4;
    float k,k1,k2,k3,k4,k5;
    float c1,c2,c3,c4,c5; /* constants for the mm expansion */
    float get_f_variable();
    ATOM *ap;
    ATOM *a_next( ); /* returns first ATOM when called with -1 */
    int a_number();
    int j,imax,inclose;
    char line[BUFSIZ];



    imax = a_number();
    if( imax <= 1 ) return

            /* get the center
            *  we always expand around it */
            xc = 0.; yc = 0.; zc = 0.;
    for( i=0; i< imax; i++)
    {
        ap = a_next(i);
        xc +=  ap->x + lambda*ap->dx;
        yc +=  ap->y + lambda*ap->dy;
        zc +=  ap->z + lambda*ap->dz;
    }
    xc /= imax;
    yc /= imax;
    zc /= imax;
    (*mmp).xc = xc;
    (*mmp).yc = yc;
    (*mmp).zc = zc;
    /* initialze the data */
    (*mmp).q000 = 0.;
    (*mmp).q100 = 0.;
    (*mmp).q010 = 0.;
    (*mmp).q001 = 0.;
    (*mmp).q200 = 0.;
    (*mmp).q020 = 0.;
    (*mmp).q002 = 0.;
    (*mmp).q101 = 0.;
    (*mmp).q110 = 0.;
    (*mmp).q011 = 0.;
    (*mmp).q300 = 0.;
    (*mmp).q030 = 0.;
    (*mmp).q003 = 0.;
    (*mmp).q210 = 0.;
    (*mmp).q120 = 0.;
    (*mmp).q201 = 0.;
    (*mmp).q102 = 0.;
    (*mmp).q021 = 0.;
    (*mmp).q012 = 0.;
    (*mmp).q111 = 0.;
#ifdef FOURTH
    (*mmp).q400 = 0.;
    (*mmp).q040 = 0.;
    (*mmp).q004 = 0.;
    (*mmp).q310 = 0.;
    (*mmp).q130 = 0.;
    (*mmp).q301 = 0.;
    (*mmp).q103 = 0.;
    (*mmp).q031 = 0.;
    (*mmp).q013 = 0.;
    (*mmp).q220 = 0.;
    (*mmp).q202 = 0.;
    (*mmp).q022 = 0.;
    (*mmp).q211 = 0.;
    (*mmp).q121 = 0.;
    (*mmp).q112 = 0.;
#endif
#ifdef FIFTH
    (*mmp).q500 = 0.;
    (*mmp).q050 = 0.;
    (*mmp).q005 = 0.;
    (*mmp).q410 = 0.;
    (*mmp).q140 = 0.;
    (*mmp).q401 = 0.;
    (*mmp).q104 = 0.;
    (*mmp).q041 = 0.;
    (*mmp).q014 = 0.;
    (*mmp).q320 = 0.;
    (*mmp).q230 = 0.;
    (*mmp).q302 = 0.;
    (*mmp).q203 = 0.;
    (*mmp).q032 = 0.;
    (*mmp).q023 = 0.;
    (*mmp).q221 = 0.;
    (*mmp).q212 = 0.;
    (*mmp).q122 = 0.;
    (*mmp).q311 = 0.;
    (*mmp).q131 = 0.;
    (*mmp).q113 = 0.;
#endif
    for( i=0; i< imax; i++)
    {
        ap = a_next(i);
        xt = ap->x + lambda *ap->dx -xc;
        yt = ap->y + lambda *ap->dy -yc;
        zt = ap->z + lambda *ap->dz -zc;
        (*mmp).q000 +=  ap->q;
        xt2 = xt*xt;
        xt3 = xt2*xt;
        xt4 = xt3*xt;
        yt2 = yt*yt;
        yt3 = yt2*yt;
        yt4 = yt3*yt;
        zt2 = zt*zt;
        zt3 = zt2*zt;
        zt4 = zt3*zt;
        (*mmp).q100 += ap->q*xt;
        (*mmp).q010 += ap->q*yt;
        (*mmp).q001 += ap->q*zt;
        (*mmp).q200 += ap->q*xt2;
        (*mmp).q020 += ap->q*yt2;
        (*mmp).q002 += ap->q*zt2;
        (*mmp).q101 += ap->q*xt*zt;
        (*mmp).q110 += ap->q*xt*yt;
        (*mmp).q011 += ap->q*yt*zt;
        (*mmp).q300 += ap->q*xt3;
        (*mmp).q030 += ap->q*yt3;
        (*mmp).q003 += ap->q*zt3;
        (*mmp).q210 += ap->q*xt2*yt;
        (*mmp).q120 += ap->q*xt*yt2;
        (*mmp).q201 += ap->q*xt2*zt;
        (*mmp).q102 += ap->q*xt*zt2;
        (*mmp).q021 += ap->q*yt2*zt;
        (*mmp).q012 += ap->q*yt*zt2;
        (*mmp).q111 += ap->q*xt*yt*zt;
#ifdef FOURTH
        (*mmp).q400 += ap->q*xt4;
        (*mmp).q040 += ap->q*yt4;
        (*mmp).q004 += ap->q*zt4;
        (*mmp).q310 += ap->q*xt3*yt;
        (*mmp).q130 += ap->q*xt*yt3;
        (*mmp).q301 += ap->q*xt3*zt;
        (*mmp).q103 += ap->q*xt*zt3;
        (*mmp).q031 += ap->q*yt3*zt;
        (*mmp).q013 += ap->q*yt*zt3;
        (*mmp).q220 += ap->q*xt2*yt2;
        (*mmp).q202 += ap->q*xt2*zt2;
        (*mmp).q022 += ap->q*yt2*zt2;
        (*mmp).q211 += ap->q*xt2*yt*zt;
        (*mmp).q121 += ap->q*xt*yt2*zt;
        (*mmp).q112 += ap->q*xt*yt*zt2;
#endif
#ifdef FIFTH
        (*mmp).q500 += ap->q*xt4*xt;
        (*mmp).q050 += ap->q*yt4*yt;
        (*mmp).q005 += ap->q*zt4*zt;
        (*mmp).q410 += ap->q*xt4*yt;
        (*mmp).q140 += ap->q*yt4*xt;
        (*mmp).q401 += ap->q*xt4*zt;
        (*mmp).q104 += ap->q*zt4*xt;
        (*mmp).q041 += ap->q*yt4*zt;
        (*mmp).q014 += ap->q*zt4*yt;
        (*mmp).q320 += ap->q*xt3*yt2;
        (*mmp).q230 += ap->q*yt3*xt2;
        (*mmp).q302 += ap->q*xt3*zt2;
        (*mmp).q203 += ap->q*zt3*xt2;
        (*mmp).q032 += ap->q*yt3*zt2;
        (*mmp).q023 += ap->q*zt3*yt2;
        (*mmp).q221 += ap->q*xt2*yt2*zt;
        (*mmp).q212 += ap->q*xt2*yt*zt2;
        (*mmp).q122 += ap->q*xt*yt2*zt2;
        (*mmp).q311 += ap->q*xt3*yt*zt;
        (*mmp).q131 += ap->q*xt*yt3*zt;
        (*mmp).q113 += ap->q*xt*yt*zt3;
#endif
    }
    /* and now (almost done with the MM setup)
    * normalize the accumulated nodal data */
#pragma _CNX no_recurrence
    /* multiplied by .5 to correct for double counting */
    k =
        xt = .5/3.;
    yt = xt/4.;
    zt = yt/5.;
    {
        (*mmp).q000 *= k;
        (*mmp).q100 *= k;
        (*mmp).q010 *= k;
        (*mmp).q001 *= k;
        (*mmp).q200 *= .5*k;
        (*mmp).q020 *= .5*k;
        (*mmp).q002 *= .5*k;
        (*mmp).q101 *= k;
        (*mmp).q110 *= k;
        (*mmp).q011 *= k;
        (*mmp).q300 *= xt*k;
        (*mmp).q030 *= xt*k;
        (*mmp).q003 *= xt*k;
        (*mmp).q210 *= 0.5*k;
        (*mmp).q120 *= 0.5*k;
        (*mmp).q201 *= 0.5*k;
        (*mmp).q102 *= 0.5*k;
        (*mmp).q021 *= 0.5*k;
        (*mmp).q012 *= 0.5*k;
        (*mmp).q111 *= k;
#ifdef FOURTH
        (*mmp).q400 *= yt*k;
        (*mmp).q040 *= yt*k;
        (*mmp).q004 *= yt*k;
        (*mmp).q310 *= xt*k;
        (*mmp).q130 *= xt*k;
        (*mmp).q301 *= xt*k;
        (*mmp).q103 *= xt*k;
        (*mmp).q031 *= xt*k;
        (*mmp).q013 *= xt*k;
        (*mmp).q220 *= .25*k;
        (*mmp).q202 *= .25*k;
        (*mmp).q022 *= .25*k;
        (*mmp).q211 *= .5*k;
        (*mmp).q121 *= .5*k;
        (*mmp).q112 *= .5*k;
#endif
#ifdef FIFTH
        (*mmp).q500 *= zt*k;
        (*mmp).q050 *= zt*k;
        (*mmp).q005 *= zt*k;
        (*mmp).q410 *= yt*k;
        (*mmp).q140 *= yt*k;
        (*mmp).q401 *= yt*k;
        (*mmp).q104 *= yt*k;
        (*mmp).q041 *= yt*k;
        (*mmp).q014 *= yt*k;
        (*mmp).q320 *= .5*xt*k;
        (*mmp).q230 *= .5*xt*k;
        (*mmp).q302 *= .5*xt*k;
        (*mmp).q203 *= .5*xt*k;
        (*mmp).q032 *= .5*xt*k;
        (*mmp).q023 *= .5*xt*k;
        (*mmp).q221 *= .25*k;
        (*mmp).q212 *= .25*k;
        (*mmp).q122 *= .25*k;
        (*mmp).q311 *= xt*k;
        (*mmp).q131 *= xt*k;
        (*mmp).q113 *= xt*k;
#endif
    }

    /* initiallization of the multipole is done !!! */
    /* therefore return for the explicity use */

    return ;
}

/* v_react calculates the reaction field
*/

int v_react( V, lambda)
float *V,lambda;
{
    ATOM *ap,*a_next();
    float fx,fy,fz,vt;
    float fxt,fyt,fzt;
    float xt,yt,zt,r,r0;
    float c1,c2,c3,c4,c5;
    float xt2,yt2,zt2;
    float xt3,yt3,zt3;
    float xt4,yt4,zt4;
    float radius,reflect; /* used for the mapping */
    float get_f_variable();
    int i,numatm,a_number();
    int rxn_common_routine( );

    MULTIPOLE *mmp;


    numatm = a_number();
    if( numatm < 1) return 1;

    mmp = malloc( sizeof( MULTIPOLE));
    if( mmp == NULL )
    {
        aerror("cannot allocate memory in v_react");
        return 0;
    }

    radius = get_f_variable("radius");
    if( radius <= one)
    {
        radius = get_radius_of_gyration();
        radius = radius *3;
        set_f_variable( "radius",radius);
    }
    reflect = get_f_variable("reflct");
    if( reflect <= one )
    {
        reflect = 79./81.;  /* swr with dielectric = 80, and 1 */
    }

    rxn_common_routine( mmp,lambda);

    /* now do the multipole expansion for the electrostatic terms */
    /* note that dielectric is included in the multipole expansion */
    xc = (*mmp).xc;
    yc = (*mmp).yc;
    zc = (*mmp).zc;
    fx = 0.;
    fy = 0.;
    fz = 0.;
    for( i=0 ; i< numatm; i++ )
    {
        fxt = 0.;
        fyt = 0.;
        fzt = 0.;
        vt = 0.;
        ap = a_next(i);
        xt = ap->x + lambda*ap->dx -xc;
        yt = ap->y + lambda*ap->dy -yc;
        zt = ap->z + lambda*ap->dz -zc;
        if( fabs(xt) < 1.e-4) xt = 1.e-4
                                       if( fabs(yt) < 1.e-4) yt = 1.e-4
                                                                      if( fabs(zt) < 1.e-4) zt = 1.e-4
                                                                                                     /* the transformation */
                                                                                                     xt = radius/xt;
        yt = radius/yt;
        zt = radius/zt;

        r = one/(xt*xt + yt*yt + zt*zt);
        r0 = sqrt(r);
        c1 =  -r*r0;
        c2 = -three*c1*r;
        c3 = -five*c2*r;
        c4 = -seven*c3*r;
        c5 = -nine*c4*r;
        xt2 = xt*xt;
        xt3 = xt2*xt;
        xt4 = xt3*xt;
        yt2 = yt*yt;
        yt3 = yt2*yt;
        yt4 = yt3*yt;
        zt2 = zt*zt;
        zt3 = zt2*zt;
        zt4 = zt3*zt;
        vt += (*mmp).q000*ap->q*r0;
        k = c1*ap->q*xt;
        vt += k*(*mmp).q100;
        fxt += k*(*mmp).q000;
        k = c1*ap->q*yt;
        vt += k*(*mmp).q010;
        fyt += k*(*mmp).q000;
        k = c1*ap->q*zt;
        vt += k*(*mmp).q001;
        fzt += k*(*mmp).q000;
        /* n=2 */
        k = (c2*xt2 +c1)*ap->q;
        vt += k*(*mmp).q200;
        fxt += k*(*mmp).q100;
        k = (c2*yt2 +c1)*ap->q;
        vt += k*(*mmp).q020;
        fyt += k*(*mmp).q010;
        k = (c2*zt2 +c1)*ap->q;
        vt += k*(*mmp).q002;
        fzt += k*(*mmp).q001;
        k = c2*xt*yt*ap->q;
        vt += k*(*mmp).q110;
        fxt += k*(*mmp).q010;
        fyt += k*(*mmp).q100;
        k = c2*xt*zt*ap->q;
        vt += k*(*mmp).q101;
        fxt += k*(*mmp).q001;
        fzt += k*(*mmp).q100;
        k = c2*yt*zt*ap->q;
        vt += k*(*mmp).q011;
        fyt += k*(*mmp).q001;
        fzt += k*(*mmp).q010;
        /* n=3 */
        k = (c3*xt3 +3*c2*xt)*ap->q;
        vt += k*(*mmp).q300;
        fxt += k*(*mmp).q200;
        k = (c3*yt3 +3*c2*yt)*ap->q;
        vt += k*(*mmp).q030;
        fyt += k*(*mmp).q020;
        k = (c3*zt3 +3*c2*zt)*ap->q;
        vt += k*(*mmp).q003;
        fzt += k*(*mmp).q002;
        k = (c3*xt2*yt+c2*yt)*ap->q;
        vt += k*(*mmp).q210;
        fxt += k*(*mmp).q110;
        fyt += k*(*mmp).q200;
        k = (c3*yt2*xt+c2*xt)*ap->q;
        vt += k*(*mmp).q120;
        fxt += k*(*mmp).q020;
        fyt += k*(*mmp).q110;
        k = (c3*xt2*zt+c2*zt)*ap->q;
        vt += k*(*mmp).q201;
        fxt += k*(*mmp).q101;
        fzt += k*(*mmp).q200;
        k = (c3*zt2*xt+c2*xt)*ap->q;
        vt += k*(*mmp).q102;
        fxt += k*(*mmp).q002;
        fzt += k*(*mmp).q101;
        k = (c3*yt2*zt+c2*zt)*ap->q;
        vt += k*(*mmp).q021;
        fyt += k*(*mmp).q011;
        fzt += k*(*mmp).q020;
        k = (c3*zt2*yt+c2*yt)*ap->q;
        vt += k*(*mmp).q012;
        fyt += k*(*mmp).q002;
        fzt += k*(*mmp).q011;
        k = (c3*zt*yt*xt)*ap->q;
        vt += k*(*mmp).q111;
        fxt += k*(*mmp).q011;
        fyt += k*(*mmp).q101;
        fzt += k*(*mmp).q110;
        /* n=4 */
#ifdef FOURTH
        k = (c4*xt4 +six*c3*(xt2) +three*c2)*ap->q;
        vt += k*(*mmp).q400;
        fxt += k*(*mmp).q300;
        k = (c4*yt4 +six*c3*(yt2) +three*c2)*ap->q;
        vt += k*(*mmp).q040;
        fyt += k*(*mmp).q030;
        k = (c4*zt4 +six*c3*(zt2) +three*c2)*ap->q;
        vt += k*(*mmp).q004;
        fzt += k*(*mmp).q003;
        k = (c4*xt3*yt + three*c3*xt*yt)*ap->q;
        vt += k*(*mmp).q310;
        fxt += k*(*mmp).q210;
        fyt += k*(*mmp).q300;
        k = (c4*yt3*xt + three*c3*xt*yt)*ap->q;
        vt += k*(*mmp).q130;
        fxt += k*(*mmp).q030;
        fyt += k*(*mmp).q120;
        k = (c4*xt3*zt + three*c3*xt*zt)*ap->q;
        vt += k*(*mmp).q301;
        fxt += k*(*mmp).q201;
        fzt += k*(*mmp).q300;
        k = (c4*zt3*yt + three*c3*xt*yt)*ap->q;
        vt += k*(*mmp).q103;
        fxt += k*(*mmp).q003;
        fzt += k*(*mmp).q102;
        k = (c4*yt3*zt + three*c3*zt*yt)*ap->q;
        vt += k*(*mmp).q031;
        fzt += k*(*mmp).q030;
        fyt += k*(*mmp).q021;
        k = (c4*zt3*yt + three*c3*zt*yt)*ap->q;
        vt += k*(*mmp).q013;
        fzt += k*(*mmp).q012;
        fyt += k*(*mmp).q003;
        k = (c4*xt2*yt2 + c3*(xt2+yt2) +c2)*ap->q;
        vt += k*(*mmp).q220;
        fxt += k*(*mmp).q120;
        fyt += k*(*mmp).q210;
        k = (c4*xt2*zt2 + c3*(xt2+zt2) +c2)*ap->q;
        vt += k*(*mmp).q202;
        fxt += k*(*mmp).q102;
        fzt += k*(*mmp).q201;
        k = (c4*zt2*yt2 + c3*(zt2+yt2) +c2)*ap->q;
        vt += k*(*mmp).q022;
        fzt += k*(*mmp).q021;
        fyt += k*(*mmp).q012;
        k = (c4*xt2*yt*zt +c3*yt*zt)*ap->q;
        vt += k*(*mmp).q211;
        fzt += k*(*mmp).q210;
        fyt += k*(*mmp).q201;
        fxt += k*(*mmp).q111;
        k = (c4*xt*yt2*zt +c3*xt*zt)*ap->q;
        vt += k*(*mmp).q121;
        fzt += k*(*mmp).q120;
        fyt += k*(*mmp).q111;
        fxt += k*(*mmp).q021;
        k = (c4*xt*yt*zt2 +c3*yt*xt)*ap->q;
        vt += k*(*mmp).q112;
        fzt += k*(*mmp).q111;
        fyt += k*(*mmp).q102;
        fxt += k*(*mmp).q012;
#endif
        /* n=5 */
#ifdef FIFTH
        k = ((c5*xt+9*c4)*xt4  +15*c3*xt)*ap->q;
        vt += k*(*mmp).q500;
        fxt += k*(*mmp).q400;
        k = ((c5*yt+9*c4)*yt4  +15*c3*yt)*ap->q;
        vt += k*(*mmp).q050;
        fyt += k*(*mmp).q040;
        k = ((c5*zt+9*c4)*zt4  +15*c3*zt)*ap->q;
        vt += k*(*mmp).q005;
        fzt += k*(*mmp).q004;
        k = (c5*xt4+six*c4*xt2 +three*c3)*yt*ap->q;
        vt += k*(*mmp).q410;
        fxt += k*(*mmp).q310;
        fyt += k*(*mmp).q400;
        k = (c5*yt4+six*c4*yt2 +three*c3)*xt*ap->q;
        vt += k*(*mmp).q140;
        fxt += k*(*mmp).q040;
        fyt += k*(*mmp).q130;
        k = (c5*xt4+six*c4*xt2 +three*c3)*zt*ap->q;
        vt += k*(*mmp).q401;
        fxt += k*(*mmp).q301;
        fzt += k*(*mmp).q400;
        k = (c5*zt4+six*c4*zt2 +three*c3)*xt*ap->q;
        vt += k*(*mmp).q104;
        fxt += k*(*mmp).q004;
        fzt += k*(*mmp).q103;
        k = (c5*yt4+six*c4*yt2 +three*c3)*zt*ap->q;
        vt += k*(*mmp).q041;
        fyt += k*(*mmp).q031;
        fzt += k*(*mmp).q040;
        k = (c5*zt4+six*c4*zt2 +three*c3)*yt*ap->q;
        vt += k*(*mmp).q014;
        fyt += k*(*mmp).q004;
        fzt += k*(*mmp).q013;
        k = (c5*xt3*yt2 +c4*(three*xt*yt2-xt3) +three*c3*xt)*ap->q;
        vt += k*(*mmp).q320;
        fxt += k*(*mmp).q220;
        fyt += k*(*mmp).q310;
        k = (c5*yt3*xt2 +c4*(three*yt*xt2-yt3) +three*c3*yt)*ap->q;
        vt += k*(*mmp).q230;
        fxt += k*(*mmp).q130;
        fyt += k*(*mmp).q220;
        k = (c5*xt3*zt2 +c4*(three*xt*zt2-xt3) +three*c3*xt)*ap->q;
        vt += k*(*mmp).q302;
        fxt += k*(*mmp).q202;
        fzt += k*(*mmp).q301;
        k = (c5*zt3*xt2 +c4*(three*zt*xt2-zt3) +three*c3*zt)*ap->q;
        vt += k*(*mmp).q203;
        fxt += k*(*mmp).q103;
        fzt += k*(*mmp).q202;
        k = (c5*yt3*zt2 +c4*(three*yt*zt2-yt3) +three*c3*yt)*ap->q;
        vt += k*(*mmp).q032;
        fyt += k*(*mmp).q022;
        fzt += k*(*mmp).q031;
        k = (c5*zt3*yt2 +c4*(three*zt*yt2-zt3) +three*c3*zt)*ap->q;
        vt += k*(*mmp).q023;
        fyt += k*(*mmp).q013;
        fzt += k*(*mmp).q022;
        k = (c5*xt2*yt2 +c4*(xt2+yt2) +c3)*zt*ap->q;
        vt += k*(*mmp).q221;
        fxt += k*(*mmp).q121;
        fyt += k*(*mmp).q211;
        fzt += k*(*mmp).q220;
        k = (c5*xt2*zt2 +c4*(xt2+zt2) +c3)*yt*ap->q;
        vt += k*(*mmp).q212;
        fxt += k*(*mmp).q112;
        fyt += k*(*mmp).q202;
        fzt += k*(*mmp).q211;
        k = (c5*zt2*yt2 +c4*(zt2+yt2) +c3)*xt*ap->q;
        vt += k*(*mmp).q122;
        fxt += k*(*mmp).q022;
        fyt += k*(*mmp).q112;
        fzt += k*(*mmp).q121;
        k = (c5*xt3+three*c4*xt)*yt*zt*ap->q;
        vt += k*(*mmp).q311;
        fxt += k*(*mmp).q211;
        fyt += k*(*mmp).q301;
        fzt += k*(*mmp).q310;
        k = (c5*yt3+three*c4*yt)*xt*zt*ap->q;
        vt += k*(*mmp).q131;
        fxt += k*(*mmp).q031;
        fyt += k*(*mmp).q121;
        fzt += k*(*mmp).q130;
        k = (c5*zt3+three*c4*zt)*yt*xt*ap->q;
        vt += k*(*mmp).q113;
        fxt += k*(*mmp).q013;
        fyt += k*(*mmp).q103;
        fzt += k*(*mmp).q112;
#endif
        *V -= vt*reflect;
    }

    free( mmp );
    return 1;

}


