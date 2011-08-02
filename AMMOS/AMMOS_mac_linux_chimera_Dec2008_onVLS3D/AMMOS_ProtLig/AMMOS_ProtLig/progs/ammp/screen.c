/* screen.c
*
* collection of routines to service screened non-locaL nonbonded potentials
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
*/
/*
*  copyright 1992, 1993, 1994 Robert W. Harrison
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

#ifdef GRACELESS
#define exp(x) ((x) > -100.? exp(x): 0.)
#endif

/* f_screen()
*
* f_screen increments the forces in the atom structures by the force
* due to the screened-nonbonded components.  NOTE THE WORD increment.
* the forces should first be zero'd.
* if not then this code will be invalid.  THIS IS DELIBERATE.
* on bigger (and better?) machines the different potential terms
* may be updated at random or in parrellel, if we assume that this routine
* will initialize the forces then we can't do this.
*/
int f_screen(lambda)
float lambda;
/*  returns 0 if error, 1 if OK */
{
    float ux,uy,uz;
    float k,r,r0,xt,yt,zt;
    float lcutoff,cutoff,get_f_variable();
    int inbond,inangle,i,test;
    ATOM *a1,*a2,*bonded[10],*angled[10];
    ATOM *a_next( ); /* returns first ATOM when called with -1 */
    int a_number(),inbuffer,imax;
    float (*buffer)[];
    int invector,atomsused,ii,jj;
    float (*vector)[];
    ATOM *(*atms)[],*(*atomall)[];
    float dielectric;
    float fx,fy,fz;
    float xt2,xt3,xt4;
    float yt2,yt3,yt4;
    float zt2,zt3,zt4;
    float alpha,ra;
    alpha = get_f_variable("alpha");
    if( alpha < 1.e-7) alpha = 1.;

    /* first update the lists
    *  this routine checks if any atom has
    *   broken the mxdq barrier and updates the
    * forces, potentials and expansions thereof */
    fv_update_nonbon( lambda);

    dielectric = get_f_variable("dielec");
    if( dielectric < 1.) dielectric = 1.;
    dielectric = 332.17752/dielectric;

    /*  get the number of atoms and allocate the memory for the array space */
    i = a_number();
    atomall = malloc( i*sizeof(ATOM *) );
    if( atomall == NULL)
    {aaerror("cannot allocate memory in f_nonbon"); return 0;}

    imax = a_number();
    for( i=0; i< imax; i++)
    {
        (*atomall)[i] = a_next(i);
    }
    for( i= 0; i< imax; i++)
    {
        fx = 0.; fy = 0.; fz = 0.;
        a1 = (*atomall)[i];
        xt = a1->dx*lambda +a1->x - a1->px;
        yt = a1->dy*lambda +a1->y - a1->py;
        zt = a1->dz*lambda +a1->z - a1->pz;


        fx = (a1->qxx*xt + a1->qxy*yt
              + a1->qxz*zt) ;
        fy = (a1->qxy*xt + a1->qyy*yt
              + a1->qyz*zt) ;
        fz = (a1->qxz*xt + a1->qyz*yt
              + a1->qzz*zt) ;
#ifdef CUBIC
        xt2 = xt*xt; yt2 = yt*yt; zt2 = zt*zt;
        fx += a1->qxxx*xt2/2. + a1->qxxy*xt*yt + a1->qxxz*xt*zt
              + a1->qxyy*yt/2. + a1->qxyz*yt*zt + a1->qxzz*zt2/2.;
        fy += a1->qxxy*xt2/2. + a1->qxyy*xt*yt + a1->qxyz*xt*zt
              + a1->qyyy*yt/2. + a1->qyyz*yt*zt + a1->qyzz*zt2/2.;
        fz += a1->qxxz*xt2/2. + a1->qxyz*xt*yt + a1->qxzz*xt*zt
              + a1->qyyz*yt/2. + a1->qyzz*yt*zt + a1->qzzz*zt2/2.;
#endif
#ifdef QUARTIC
        xt3 = xt*xt2; yt3 = yt*yt2; zt3 = zt*zt2;
        fx +=  a1->qxxxx*xt3/6. + a1->qxxxy*xt2*yt/2. + a1->qxxxz*xt2*zt/2.
               + a1->qxxyy*xt*yt/2. + a1->qxxyz*xt*yt*zt + a1->qxxzz*xt*zt2/2.
               + a1->qxyyy*yt3/6. + a1->qxyyz*yt2*zt/2. + a1->qxyzz*yt*zt2/2.
               + a1->qxzzz*zt3/6.;
        fy +=  a1->qxxxy*xt3/6. + a1->qxxyy*xt2*yt/2. + a1->qxxyz*xt2*zt/2.
               + a1->qxyyy*xt*yt/2. + a1->qxyyz*xt*yt*zt + a1->qxyzz*xt*zt2/2.
               + a1->qyyyy*yt3/6. + a1->qyyyz*yt2*zt/2. + a1->qyyzz*yt*zt2/2.
               + a1->qyzzz*zt3/6.;
        fz +=  a1->qxxxz*xt3/6. + a1->qxxyz*xt2*yt/2. + a1->qxxzz*xt2*zt/2.
               + a1->qxyyz*xt*yt/2. + a1->qxyzz*xt*yt*zt + a1->qxzzz*xt*zt2/2.
               + a1->qyyyz*yt3/6. + a1->qyyzz*yt2*zt/2. + a1->qyzzz*yt*zt2/2.
               + a1->qzzzz*zt3/6.;
#endif
#ifdef QUINTIC
        xt4 = xt*xt3; yt4 = yt*yt3; zt4 = zt*zt3;
        fx += a1->qxxxxx*xt4/24. + a1->qxxxxy*xt3*yt/6. + a1->qxxxxz*xt3*zt/6.
              + a1->qxxxyy*xt2*yt2/4. + a1->qxxxyz*xt2*yt*zt/2. + a1->qxxxzz*xt2*zt2/4.
              + a1->qxxyyy*xt*yt3/6. + a1->qxxyyz*xt*yt2*zt/2. + a1->qxxyzz*xt*yt*zt2/2.
              + a1->qxxzzz*xt*zt3/6. + a1->qxyyyy*yt4/24. + a1->qxyyyz*yt3*zt/6.
              + a1->qxyyzz*yt2*zt2/4. + a1->qxyzzz*yt*zt3/6. + a1->qxzzzz*zt4/24.;
        fy += a1->qxxxxy*xt4/24. + a1->qxxxyy*xt3*yt/6. + a1->qxxxyz*xt3*zt/6.
              + a1->qxxyyy*xt2*yt2/4. + a1->qxxyyz*xt2*yt*zt/2. + a1->qxxyzz*xt2*zt2/4.
              + a1->qxyyyy*xt*yt3/6. + a1->qxyyyz*xt*yt2*zt/2. + a1->qxyyzz*xt*yt*zt2/2.
              + a1->qxyzzz*xt*zt3/6. + a1->qyyyyy*yt4/24. + a1->qyyyyz*yt3*zt/6.
              + a1->qyyyzz*yt2*zt2/4. + a1->qyyzzz*yt*zt3/6. + a1->qyzzzz*zt4/24.;
        fz += a1->qxxxxz*xt4/24. + a1->qxxxyz*xt3*yt/6. + a1->qxxxzz*xt3*zt/6.
              + a1->qxxyyz*xt2*yt2/4. + a1->qxxyzz*xt2*yt*zt/2. + a1->qxxzzz*xt2*zt2/4.
              + a1->qxyyyz*xt*yt3/6. + a1->qxyyzz*xt*yt2*zt/2. + a1->qxyzzz*xt*yt*zt2/2.
              + a1->qxzzzz*xt*zt3/6. + a1->qyyyyz*yt4/24. + a1->qyyyzz*yt3*zt/6.
              + a1->qyyzzz*yt2*zt2/4. + a1->qyzzzz*yt*zt3/6. + a1->qzzzzz*zt4/24.;
#endif
        a1->fx += fx  + a1->dpx;
        a1->fy += fy  + a1->dpy;
        a1->fz += fz  + a1->dpz;
        /* do the close atoms */
        for( jj=0; jj< NCLOSE; jj++)
        { if( a1->close[jj] == NULL) break; }
        for( ii=0; ii< jj;ii++)
        {
            a2 = a1->close[ii];
            /* note ux is backwards from below */
            ux = (a2->dx -a1->dx)*lambda +(a2->x -a1->x);
            uy = (a2->dy -a1->dy)*lambda +(a2->y -a1->y);
            uz = (a2->dz -a1->dz)*lambda +(a2->z -a1->z);
            r = ux*ux + uy*uy + uz*uz; r0 = sqrt(r);
            ux = ux/r0; uy = uy/r0; uz = uz/r0;
            /* standard point atoms
                    k = -dielectric*a1->q*a2->q/r;
            */
            /* 1s 1s  screened atoms */
            k = dielectric*a1->q*a2->q;
            ra = r0*alpha ;
            /*
            * the full term
            	k *= (1.-exp(-ra)*(1.+11/16*ra + 3/16*ra*ra + ra*ra*ra/48))/r0;
            *  the econimized and derivatived version
            */
            k *= -(1.-exp(-ra)*( (((ra+9.)*ra + 33.)*ra)/48 +1.))/r
                 + (alpha*exp(-ra)*(
                        ((ra+6.)*ra + 15.)*ra/48 +1.  -   11./16
                    ) )/r0;
            /*
            ((ra+9.)*ra + 33.)*ra/48 +1.  -  ((ra+6.)*ra+ 11.)/16

            + (alpha*exp(-ra)*( (((ra+9.)*ra + 33.)*ra/48 +1.))/r0
            - (alpha*exp(-ra)*( ((ra+6.)*ra+ 11.)/16)))/r0;
            */
            r = r*r*r;
            k += a1->a*a2->a/r/r0*6;
            k -= a1->b*a2->b/r/r/r0*12;
            a1->fx += ux*k;
            a1->fy += uy*k;
            a1->fz += uz*k;
            a2->fx -= ux*k;
            a2->fy -= uy*k;
            a2->fz -= uz*k;
        }
    }

    a_inactive_f_zero();
    free( atomall); return 1;

}
/* v_screen()
* this function sums up the potentials
* for the atoms defined in the nonbon data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int v_screen( V, lambda )
float *V,lambda;
{
    float r,r0,xt,yt,zt;
    float get_f_variable();
    float k;
    int inbond,inangle,i;
    ATOM *a1,*a2,*bonded[10],*angled[10];
    ATOM *a_next( ); /* returns first ATOM when called with -1 */
    int a_number(),inbuffer;
    int invector,atomsused,ii,jj,imax;
    float (*vector)[];
    float vx;
    float k2;
    ATOM *(*atomall)[];
    float dielectric;
    float xt2,xt3,xt4,xt5;
    float yt2,yt3,yt4,yt5;
    float zt2,zt3,zt4,zt5;
    float alpha;

    alpha = get_f_variable("alpha");
    if( alpha < 1.e-7) alpha = 1.;

    fv_update_nonbon( lambda);

    dielectric = get_f_variable("dielec");
    if( dielectric < 1.) dielectric = 1.;
    dielectric = 332.17752/dielectric;

    /*  get the number of atoms and allocate the memory for the array space */
    i = a_number();
    atomall = malloc( i*sizeof(ATOM *) );
    if( atomall == NULL)
    {aaerror("cannot allocate memory in v_nonbon"); return 0;}

    imax = a_number();
    for( i=0; i< imax; i++)
    {
        (*atomall)[i] = a_next(i);
    }
    for( i= 0; i< imax; i++)
    {
        a1 = (*atomall)[i];
        vx = a1->VP;
        xt = a1->dx*lambda +a1->x - a1->px;
        yt = a1->dy*lambda +a1->y - a1->py;
        zt = a1->dz*lambda +a1->z - a1->pz;
        vx -= (a1->dpx*xt + a1->dpy*yt
               + a1->dpz*zt) ;
        vx -= ( (xt*(.5*a1->qxx*xt + a1->qxy*yt + a1->qxz*zt)
                 + yt*(.5*a1->qyy*yt + a1->qyz*zt) + .5*zt*a1->qzz*zt));
#ifdef CUBIC
        xt2 = xt*xt; yt2 = yt*yt;  zt2 = zt*zt;
        xt3 = xt2*xt; yt3 = yt2*yt; zt3 = zt2*zt;

        vx -= a1->qxxx*xt3/6. + a1->qxxy*xt2*yt/2 + a1->qxxz*xt2*zt/2
              + a1->qxyy*xt*yt2/2 + a1->qxyz*xt*yt*zt + a1->qxzz*xt*zt2/2
              + a1->qyyy*yt3/6 + a1->qyyz*yt2*zt/2 + a1->qyzz*yt*zt2/2
              + a1->qzzz*zt3/6.;
#endif
#ifdef QUARTIC
        xt4 = xt3*xt; yt4 = yt3*yt; zt4 = zt3*zt;
        vx -= a1->qxxxx*xt4/24. + a1->qxxxy*xt3*yt/6. + a1->qxxxz*xt3*yt/6. + a1->qxxyy*xt2*yt2/4.
              + a1->qxxyz*xt2*yt*zt/2. + a1->qxxzz*xt2*zt2/4. + a1->qxyyy*xt*yt3/6.
              + a1->qxyyz*xt*yt2*zt/2. + a1->qxyzz*xt*yt*zt2/2. + a1->qxzzz*xt*zt3/6.
              + a1->qyyyy*yt4/24. + a1->qyyyz*yt3*zt/6. + a1->qyyzz*yt2*zt2/4. + a1->qyzzz*yt*zt3/6.
              + a1->qzzzz*zt4/24.;
#endif
#ifdef QUINTIC
        xt5 = xt4*xt; yt5 = yt4*yt; zt5 = zt4*zt;
        vx -= a1->qxxxxx*xt5/120. + a1->qxxxxy*xt4*yt/24. + a1->qxxxxz*xt4*zt/24.
              + a1->qxxxyy*xt3*yt2/12. + a1->qxxxyz*xt3*yt*zt/6. + a1->qxxxzz*xt3*zt2/12.
              + a1->qxxyyy*xt2*yt3/12. + a1->qxxyyz*xt2*yt2*zt/4. + a1->qxxyzz*xt2*yt*zt2/4.
              + a1->qxxzzz*xt2*zt3/12. + a1->qxyyyy*xt*yt4/24.  + a1->qxyyyz*xt*yt3*zt/6.
              + a1->qxyyzz*xt*yt2*zt2/4. + a1->qxyzzz*xt*yt*zt3/6. + a1->qxzzzz*xt*zt4/24.
              + a1->qyyyyy*yt5/120. + a1->qyyyyz*yt4*zt/24 + a1->qyyyzz*yt3*zt2/12.
              + a1->qyyzzz*yt2*zt3/12. + a1->qyzzzz*yt*zt4/24. + a1->qzzzzz*zt5/120.;

#endif
        /* do the close atoms */
        for( jj=0; jj< NCLOSE; jj++)
        { if( a1->close[jj] == NULL) break; }
        for( ii=0; ii< jj;ii++)
        {
            a2 = a1->close[ii];
            xt = (a2->dx -a1->dx)*lambda +(a2->x -a1->x);
            yt = (a2->dy -a1->dy)*lambda +(a2->y -a1->y);
            zt = (a2->dz -a1->dz)*lambda +(a2->z -a1->z);
            r = xt*xt + yt*yt + zt*zt; r0 = sqrt(r);
            /*      xt = xt/r0; yt = yt/r0; zt = zt/r0;
            */
            /* standard point atom
                    k = dielectric*a1->q*a2->q/r0;
            */
            /* 1s 1s  screened atoms */
            k = dielectric*a1->q*a2->q/r0;
            r0 *= alpha ;
            /*
            * the full term
            	k *= (1.-exp(-r0)*(1.+11/16*r0 + 3/16*r0*r0 + r0*r0*r0/48));
            *  the econimized version
            */
            k *= 1.-exp(-r0)*( (((r0+9.)*r0 + 33.)*r0)/48 +1.);
            r = r*r*r;
            k -= a1->a*a2->a/r;
            k += a1->b*a2->b/r/r;
            vx += k;
        }
        *V += vx;
    }
    a_inactive_f_zero();
    free( atomall); return 1;
}

/*zone_screen()
* this function sums up the potentials
* for the atoms defined in the nonbon data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int zone_screen( V, lambda ,alist, inalist )
float *V,lambda;
ATOM *( *alist)[] ;
int inalist;
{
    float r,r0,xt,yt,zt;
    float lcutoff,cutoff,get_f_variable();
    int inbond,inangle,i,ii;
    ATOM *a1,*a2;
    ATOM *a_next( ); /* returns first ATOM when called with -1 */
    float dielectric,ve,va,vh;
    ATOM *a_m_serial();

    /* nonbonded potentials
    * do a double loop starting from the first atom to the 
    * last 
    * then from the second to the last 
    * etc
    *
    * also check to avoid bonded and 1-3 bonded atoms
    */
    if( inalist <= 0 ) return 1;
    dielectric = get_f_variable("dielec");
    if( dielectric < 1.) dielectric = 1.;
    dielectric = 332.17752/dielectric;
    cutoff = get_f_variable("cutoff");
    if( cutoff < 1.) cutoff = 1.e10;
    lcutoff = -cutoff;
    for( ii=0; ii< inalist; ii++)
    {
        a1 = (*alist)[ii];
        if( a1 == NULL ) goto NOTANATOM;
        ve = 0.; va = 0.; vh = 0.;
        a2 = a_next(-1);
        /*
        *	for(i = 0; i< a1->dontuse; i++)
        *	printf("%d ",a1->excluded[i]->serial);
        *	printf("\n");
        */
        /*
        	while(  (a2->next != a2) && (a2->next != NULL))
        	*/
        while(  (a2 != NULL) && (a2->next != NULL) && a2->next != a2)
        {
            /* goto SKIP is used because this is one case where it makes sense */
            /*	if( a2 == a1) break;  */
            if( a2 == a1) goto SKIP;
            for(i = 0; i< a1->dontuse; i++)
                if( a2 == a1->excluded[i]) goto SKIP;
            /* non - bonded are only used when the atoms arent bonded */

            xt = (a1->x -a2->x) + lambda*(a1->dx - a2->dx);
            if( (xt > cutoff) || (xt < lcutoff) ) goto SKIP;
            yt = (a1->y -a2->y) + lambda*(a1->dy - a2->dy);
            if( (yt > cutoff) || (yt < lcutoff) ) goto SKIP;
            zt = (a1->z -a2->z) + lambda*(a1->dz - a2->dz);
            if( (zt > cutoff) || (zt < lcutoff) ) goto SKIP;

            r = xt*xt+yt*yt+zt*zt;
            if( r < 1.) r = 1.;

            r0 = sqrt(r); r = r*r*r ;
            /* debugging
            	printf(" %d %d %f %f %f \n", a1->serial,a2->serial,a1->q,a2->q,
            	332.17752*a1->q*a2->q/r0);
            */
            ve += dielectric*a1->q*a2->q/r0;
            va -= a1->a*a2->a/r;
            vh += a1->b*a2->b/r/r;

SKIP:
            /*	if( a2->next == a1) break; */
            if( a2->next == a2) break;
            a2 = a2->next;
        }
        *V += ve + va + vh;
NOTANATOM:
        i = i;
    }
    return 1;

}


/* a_screen()
* this function sums up the potentials
* for the atoms defined in the nonbon data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int a_screen( V, lambda,ilow,ihigh,op )
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
    float ra,k,alpha;

    alpha = get_f_variable("alpha");
    if( alpha < 1.e-7) alpha = 1.;

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
            /*	if( a2 == a1) goto SKIP;  */
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
            if( r < .01) r = .01;

            r0 = sqrt(r); r = r*r*r ;
            /* debugging
            *	printf(" %d %d %f %f %f \n", a1->serial,a2->serial,a1->q,a2->q,
            *	332.17752*a1->q*a2->q/r0);
            */
            /*
            *	ve += dielectric*a1->q*a2->q/r0; 
            */
            /* 1s 1s  screened atoms */
            k = dielectric*a1->q*a2->q/r0;
            ra = r0* alpha ;
            /*
            * the full term
            *	k *= (1.-exp(-ra)*(1.+11/16*ra + 3/16*ra*ra + ra*ra*ra/48));
            *  the econimized version
            */
            k *= 1.-exp(-ra)*( (((ra+9.)*ra + 33.)*ra)/48 +1.);
            ve += k;
            va -= a1->a*a2->a/r;
            vh += a1->b*a2->b/r/r;
            if( a2->serial < ilow || a2->serial > ihigh)
            {
                /*
                	vel += dielectric*a1->q*a2->q/r0; 
                */
                vel += k;
                val -= a1->a*a2->a/r;
                vhl += a1->b*a2->b/r/r;
            }

SKIP:
            /*	if( a2->next == a1) break; */
            if( a2->next == a2) break;
            a2 = a2->next;
        }
        fprintf(op,"Vnonbon internal %d Eq %f E6 %f E12 %f\n",
                ii,ve-vel,va-val,vh-vhl);
        fprintf(op,"Vnonbon external %d Eq %f E6 %f E12 %f\n",ii,vel,val,vhl);
        fprintf(op,"Vnonbon total    %d Eq %f E6 %f E12 %f\n",ii,ve,va,vh);
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

