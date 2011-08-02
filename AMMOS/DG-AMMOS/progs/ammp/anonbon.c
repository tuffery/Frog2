/* nonbon.c
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
* but otherwise is self-contained. Note the hooks for Non-nonboned potentials
*/

/* a_nonbon()
* this function sums up the potentials
* for the atoms defined in the nonbon data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int a_nonbon( V, lambda,ilow,ihigh,op )
float *V,lambda;
int ilow,ihigh;
FILE *op;
{
    float r,r0,xt,yt,zt;
    float lcutoff,cutoff,get_f_variable();
    int inbond,inangle,i,ii;
    ATOM *a1,*a2;
    ATOM *a_next( ); /* returns first ATOM when called with -1 */
    float dielectric,ve,va,vh;
    float vel,val,vhl;
    float vtint,vtout, vtt;
    ATOM *a_m_serial();

    /* nonbonded potentials
    * do a double loop starting from the first atom to the 
    * last 
    * then from the second to the last 
    * etc
    *
    * also check to avoid bonded and 1-3 bonded atoms
    */
    dielectric = get_f_variable("dielec");
    if( dielectric < 1.) dielectric = 1.;
    dielectric = 332.17752/dielectric;
    cutoff = get_f_variable("cutoff");
    if( cutoff < 1.) cutoff = 1.e10;
    lcutoff = -cutoff;
    vtint = 0.; vtout = 0.; vtt = 0.;
    for( ii=ilow; ii<=ihigh; ii++)
    {
        a1 = a_m_serial(ii);
        if( a1 == NULL ) goto NOTANATOM;
        ve = 0.; va = 0.; vh = 0.;
        vel = 0.; val = 0.; vhl = 0.;
        a2 = a_next(-1);
        /*
        *	for(i = 0; i< a1->dontuse; i++)
        *	printf("%d ",a1->excluded[i]->serial);
        *	printf("\n");
        */
        /*
        	while(  (a2->next != a2) && (a2->next != NULL))
        	*/
        while(   (a2->next != NULL))
        {
            /* goto SKIP is used because this is one case where it makes sense */
            /*	if( a2 == a1) break;  */
            /*	if( a2 == a1) goto SKIP; */
            for(i = 0; i< a1->dontuse; i++)
                if( a2 == a1->excluded[i]) goto SKIP;
            /* non - bonded are only used when the atoms arent bonded */

            if( lambda == 0.)
            {
                xt = (a1->x - a2->x);
                if( (xt > cutoff) || (xt < lcutoff) ) goto SKIP;
                yt =  (a1->y - a2->y);
                if( (yt > cutoff) || (yt < lcutoff) ) goto SKIP;
                zt =  (a1->z - a2->z);
                if( (zt > cutoff) || (zt < lcutoff) ) goto SKIP;
            } else
            {
                xt = (a1->x - a2->x +lambda*(a1->dx - a2->dx));
                if( (xt > cutoff) || (xt < lcutoff) ) goto SKIP;
                yt = (a1->y - a2->y +lambda*(a1->dy - a2->dy));
                if( (yt > cutoff) || (yt < lcutoff) ) goto SKIP;
                zt = (a1->z - a2->z +lambda*(a1->dz - a2->dz));
                if( (zt > cutoff) || (zt < lcutoff) ) goto SKIP;
            }
            r = xt*xt+yt*yt+zt*zt;
            /*	if( r < 1.) r = 1.;  */

            r0 = sqrt(r); r = r*r*r ;
            /* debugging
            *	printf(" %d %d %f %f %f \n", a1->serial,a2->serial,a1->q,a2->q,
            *	332.17752*a1->q*a2->q/r0);
            */
            ve += dielectric*a1->q*a2->q/r0;
            va -= a1->a*a2->a/r;
            vh += a1->b*a2->b/r/r;
            if( a2->serial < ilow || a2->serial > ihigh)
            {
                vel += dielectric*a1->q*a2->q/r0;
                val -= a1->a*a2->a/r;
                vhl += a1->b*a2->b/r/r;
            }

SKIP:
            /*	if( a2->next == a1) break; */
            if( a2->next == a2) break;
            a2 = a2->next;
        }
        fprintf(op,"Vnonbon internal %s %d Eq %f E6 %f E12 %f\n",
                a1->name,ii,ve-vel,va-val,vh-vhl);
        fprintf(op,"Vnonbon external %s %d Eq %f E6 %f E12 %f\n",a1->name
                ,ii,vel,val,vhl);
        fprintf(op,"Vnonbon total    %s %d Eq %f E6 %f E12 %f\n",a1->name
                ,ii,ve,va,vh);
        *V += ve + va + vh;
        vtint += ve -vel+ va -val + vh -vhl;
        vtout += vel + val + vhl;
        vtt  += ve + va + vh;
NOTANATOM:
        i = i;
    }
    fprintf(op," Vnonbon total internal %f \n",vtint);
    fprintf(op," Vnonbon total external %f \n",vtout);
    fprintf(op," Vnonbon total          %f \n",vtt);
    return 1;

}

