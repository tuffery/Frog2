/* perioldic.c
*
* collection of routines to service periodic  potentials
*  nonbonded potentials with a periodic unit cell
*  version 0. only orthorhombic cells 
*
* POOP (Poor-mans Object Oriented Programming) using scope rules
*
* the routines for potential value, force and (eventually) second
* derivatives are here also
*
* force and 2nd derivative routines assume zero'd arrays for output
* this allows for parralellization if needed (on a PC?)
*
*
* note that the non-bonded information is in the ATOM structures 
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
* but otherwise is self-contained. Note the hooks for Non-periodiced potentials
*/

/* v_periodic()
* this function sums up the potentials
* for the atoms defined in the periodic data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int v_periodic( V, lambda )
float *V,lambda;
{
    float r,r0,xt,yt,zt,temp;
    float lcutoff,cutoff,get_f_variable();
    float acell,bcell,ccell;
    float am,bm,cm;
    int inbond,inangle,i,ibonded;
    ATOM *a1,*a2,*bonded[10],*angled[10];
    ATOM *a_next( ); /* returns first ATOM when called with -1 */

    /* nonbonded potentials
    * do a double loop starting from the first atom to the 
    * last 
    * then from the second to the last 
    * etc
    *
    * also check to avoid bonded and 1-3 bonded atoms
    */
    cutoff = get_f_variable("cutoff");
    acell = get_f_variable("acell");
    bcell = get_f_variable("bcell");
    ccell = get_f_variable("ccell");
    if( (acell<1.) || (bcell<1.) || (ccell<1.))
    {
        aaerror("periodic - unit cell must be defined ");
        aaerror("use setf acell <value(>1.)> ");
        aaerror("use setf bcell <value(>1.)> ");
        aaerror("use setf ccell <value(>1.)> ");
        return 0;
    }
    am = acell/2;bm = bcell/2; cm = ccell/2;
    if( cutoff < 1.) cutoff = 1.e10;
    lcutoff = -cutoff;
    a1 = a_next(-1);
    while( (a1->next != a1) && (a1->next != NULL))
    {
        a1 = a1->next;
        get_bond( a1,bonded,10,&inbond);
        get_angle( a1,angled,10,&inangle);
        /*	a2 = a1->next; */
        a2 = a_next(-1);
        while(  (a2->next != a2) && (a2->next != NULL))
        {
            ibonded = 0;
            if( a1 == a2) ibonded = 1;
            for(i = 0; i< inbond; i++)
                if( a2 == bonded[i]) ibonded++;
            for(i = 0; i< inangle; i++)
                if( a2 == angled[i]) ibonded++;
            /* when bonded the closest image atom should be used  */
            if( ibonded > 0)
            {
                xt = (a1->x - a2->x +lambda*(a1->dx - a2->dx)+acell);
                temp = (a1->x - a2->x +lambda*(a1->dx - a2->dx)-acell);
                if( abs(temp) < abs(xt) ) xt = temp;
                yt = (a1->y - a2->y +lambda*(a1->dy - a2->dy)+bcell);
                temp = (a1->y - a2->y +lambda*(a1->dy - a2->dy)-bcell);
                if( abs(temp) < abs(yt) ) yt = temp;
                zt = (a1->z - a2->z +lambda*(a1->dz - a2->dz)+ccell);
                temp = (a1->z - a2->z +lambda*(a1->dz - a2->dz)-ccell);
                if( abs(temp) < abs(zt) ) zt = temp;
            }else
                /* otherwise just use the closest atom or image atom  */
            {
                if( lambda == 0.)
                {
                    xt =  (a1->x - a2->x);
                    yt =  (a1->y - a2->y);
                    zt =  (a1->z - a2->z);
                } else
                {
                    xt = (a1->x - a2->x +lambda*(a1->dx - a2->dx));
                    yt = (a1->y - a2->y +lambda*(a1->dy - a2->dy));
                    zt = (a1->z - a2->z +lambda*(a1->dz - a2->dz));
                }
CHECK:
                if( xt > acell) { xt = xt-acell; goto CHECK;}
                if( yt > bcell) { yt = yt-bcell; goto CHECK;}
                if( zt > ccell) { zt = zt-ccell; goto CHECK;}
                if( xt < 0.) { xt = xt+acell; goto CHECK;}
                if( yt < 0.) { yt = yt+bcell; goto CHECK;}
                if( zt < 0.) { zt = zt+ccell; goto CHECK;}
                if( xt > am) {xt = - acell + xt; }
                if( yt > bm) {yt = - bcell + yt; }
                if( zt > cm) {zt = - ccell + zt; }
            }
            if( (xt > cutoff) || (xt < lcutoff) ) goto SKIP;
            if( (yt > cutoff) || (yt < lcutoff) ) goto SKIP;
            if( (zt > cutoff) || (zt < lcutoff) ) goto SKIP;
            r = xt*xt+yt*yt+zt*zt;
            if( r < 1.) r = 1.;
            r0 = sqrt(r); r = r*r*r ;
            /* debugging
            *	printf(" %d %d %f %f %f \n", a1->serial,a2->serial,a1->q,a2->q,
            *	332.17752*a1->q*a2->q/r0);
            */
            *V += 332.17752*a1->q*a2->q/r0;
            *V -= a1->a*a2->a/r;
            *V += a1->b*a2->b/r/r;
SKIP:
            if( a2 == a1) break;
            a2 = a2->next;
        }
    }
    return 1;

}
/* f_periodic()
*
* f_periodic increments the forces in the atom structures by the force
* due to the periodic components.  NOTE THE WORD increment.
* the forces should first be zero'd.
* if not then this code will be invalid.  THIS IS DELIBERATE.
* on bigger (and better?) machines the different potential terms
* may be updated at random or in parrellel, if we assume that this routine
* will initialize the forces then we can't do this.
*/
int f_periodic(lambda)
float lambda;
/*  returns 0 if error, 1 if OK */
{
    float ux,uy,uz;
    float k,r,r0,xt,yt,zt,temp;
    float lcutoff,cutoff,get_f_variable();
    float acell,bcell,ccell;
    float am,bm,cm;
    int inbond,inangle,i,ibonded;
    ATOM *a1,*a2,*bonded[10],*angled[10];
    ATOM *a_next( ); /* returns first ATOM when called with -1 */

    /* periodicded potentials
    * do a double loop starting from the first atom to the 
    * last 
    * then from the second to the last 
    * etc
    *
    * also check to avoid bonded and 1-3 bonded atoms
    */
    cutoff = get_f_variable("cutoff");
    acell = get_f_variable("acell");
    bcell = get_f_variable("bcell");
    ccell = get_f_variable("ccell");
    if( (acell<1.) || (bcell<1.) || (ccell<1.))
    {
        aaerror("periodic - unit cell must be defined ");
        aaerror("use setf acell <value(>1.)> ");
        aaerror("use setf bcell <value(>1.)> ");
        aaerror("use setf ccell <value(>1.)> ");
        return 0;
    }
    am = acell/2;bm = bcell/2; cm = ccell/2;
    if( cutoff < 1.) cutoff = 1.e10;
    lcutoff = -cutoff;
    a1 = a_next(-1);
    while( (a1->next != a1) && (a1->next != NULL))
    {
        a1 = a1->next;
        get_bond( a1,bonded,10,&inbond);
        get_angle( a1,angled,10,&inangle);
        a2 = a_next(-1);
        while( (a2->next != a2) && (a2->next != NULL))
        {
            ibonded = 0;
            if( a1 == a2) ibonded = 1;
            for(i = 0; i< inbond; i++)
                if( a2 == bonded[i]) ibonded++;
            for(i = 0; i< inangle; i++)
                if( a2 == angled[i]) ibonded++;
            /* when bonded the closest image atom should be used  */
            if( ibonded > 0)
            {
                ux = (a1->x - a2->x +lambda*(a1->dx - a2->dx)+acell);
                temp = (a1->x - a2->x +lambda*(a1->dx - a2->dx)-acell);
                if( abs(temp) < abs(ux) ) ux = temp;
                uy = (a1->y - a2->y +lambda*(a1->dy - a2->dy)+bcell);
                temp = (a1->y - a2->y +lambda*(a1->dy - a2->dy)-bcell);
                if( abs(temp) < abs(uy) ) uy = temp;
                uz = (a1->z - a2->z +lambda*(a1->dz - a2->dz)+ccell);
                temp = (a1->z - a2->z +lambda*(a1->dz - a2->dz)-ccell);
                if( abs(temp) < abs(uz) ) uz = temp;
            }else
                /* otherwise just use the closest atom or image atom  */
            {
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
CHECK:
                if( ux > acell) { ux = ux-acell; goto CHECK;}
                if( uy > bcell) { uy = uy-bcell; goto CHECK;}
                if( uz > ccell) { uz = uz-ccell; goto CHECK;}
                if( ux < 0.) { ux = ux+acell; goto CHECK;}
                if( uy < 0.) { uy = uy+bcell; goto CHECK;}
                if( uz < 0.) { uz = uz+ccell; goto CHECK;}
                if( ux > am) {ux = - acell +ux; }
                if( uy > bm) {uy = - bcell +uy; }
                if( uz > cm) {uz = - ccell +uz; }
            }
            if( (ux > cutoff) || (ux < lcutoff) ) goto SKIP;
            if( (uy > cutoff) || (uy < lcutoff) ) goto SKIP;
            if( (uz > cutoff) || (uz < lcutoff) ) goto SKIP;
            r = ux*ux + uy*uy + uz*uz;
            /* watch for FP errors*/
            if( r <= 1.)
            { r = 1.; }
            r0 = sqrt(r); ux = ux/r0; uy = uy/r0; uz = uz/r0;
            k = -332.17752*a1->q*a2->q/r;
            r = r*r*r;
            /*	k += a1->a*a2->a/r/r0*6;
            	k -= a1->b*a2->b/r/r/r0*12;

            */
            k += ( a1->a*a2->a*6. - a1->b*a2->b/r*12.)/r/r0;
            a1->fx += ux*k;
            a1->fy += uy*k;
            a1->fz += uz*k;
            a2->fx -= ux*k;
            a2->fy -= uy*k;
            a2->fz -= uz*k;
SKIP:

            if( a2 == a1) break;
            a2 = a2->next;
        }
    }
    return 1;

}
