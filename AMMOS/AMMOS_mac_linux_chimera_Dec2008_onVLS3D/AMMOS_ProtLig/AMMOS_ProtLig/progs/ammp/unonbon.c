/* unonbon.c
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
* unonbon uses a 'use list' for the interactions
*
*  this is physically incorrect, but often done.
*  the use list is coded as
*  ap ... bp's for interaction ... ap
* a single array is malloc'd, it is 20 a_number() long
*  the global variable (in variable storage) nbdeep will allowriding
*  every nsteps ( 10 default ) will recalculate the list
*  again this may be over ridden with  nbstep
*  and if cutoff is not set (== 0) these routines silently call the
*  regular routines which don't care about cutoff
*  will always redo the list if a_number() changes
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

/* misc includes - ANSI and some are just to be safe */
#define ANSI 1
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

/* v_nonbon()
* this function sums up the potentials
* for the atoms defined in the nonbon data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int u_v_nonbon( V, lambda )
float *V,lambda;
{
    float r,r0,xt,yt,zt;
    float lcutoff,cutoff,get_f_variable();
    float rdebye;
    int inbond,inangle,i;
    ATOM *a1,*a2,*bonded[10],*angled[10];
    ATOM *a_next( ); /* returns first ATOM when called with -1 */
    ATOM *(*use)[];
    int uselist(),nuse,used;
    ATOM *cp,*bp;
    int a_number(),inbuffer;
    float (*buffer)[],xx,yy,zz;
    int invector,atomsused,ii,jj,imax;
    float (*vector)[];
    ATOM *(*atms)[];
    float dielectric;

    /* nonbonded potentials
    * do a double loop starting from the first atom to the 
    * last 
    * then from the second to the last 
    * etc
    *
    * also check to avoid bonded and 1-3 bonded atoms
    */
    cutoff = get_f_variable("cutoff");
    if( cutoff < 1.)
    {
        v_nonbon( V,lambda);
        return 1;
    }
    rdebye = cutoff/2.;
    dielectric = get_f_variable("dielec");
    if( dielectric < 1.) dielectric = 1.;
    dielectric = 332.17752/dielectric;
    if( !uselist( &use,&nuse,cutoff) )return 0;
    /*  get the number of atoms and allocate the memory for the array space */
    i = a_number();
    buffer = malloc( 3*i*sizeof(float) );
    if( buffer == NULL)
    { aaerror("cannot allocate memory in u_v_nonbon\n"); return 0;}
    vector = malloc( i*sizeof(float) );
    if( vector == NULL)
    { aaerror("cannot allocate memory in u_v_nonbon\n"); return 0;}
    atms = malloc( i*sizeof(ATOM *) );
    if( atms == NULL)
    { aaerror("cannot allocate memory in u_v_nonbon\n"); return 0;}


    a1 = a_next(-1);
    a1 = a1->next;
    imax = a_number();
    used = 0;
    for( jj=1; jj<imax; jj++,a1=bp)
    {
        bp = a1->next;
        inbuffer = 0;
        if( (*use)[used] == a1)
        { used += 1;}
        else { aaerror("error in uselist - must abort"); return 0;}
        while( (*use)[used] != a1)
        {
            (*atms)[inbuffer++] = (*use)[used];
            used += 1;
        }
        used += 1;
        /* (*atms) now contains the list of atoms to be  done
        *  there are inbuffer of them
        *  of course inbuffer can be zero so we must check for that
        */
        if( inbuffer > 0)
        {
#pragma _CNX no_recurrence 
            for( i=0; i< inbuffer; i++)
            {
                (*buffer)[3*i  ] = (*atms)[i]->x;
                (*buffer)[3*i+1] = (*atms)[i]->y;
                (*buffer)[3*i+2] = (*atms)[i]->z;
            }
            if( lambda != 0.)
            {
#pragma _CNX no_recurrence 
                for( i=0; i< inbuffer; i++)
                {
                    (*buffer)[3*i  ] = (*atms)[i]->x +(*atms)[i]->dx*lambda;
                    (*buffer)[3*i+1] = (*atms)[i]->y +(*atms)[i]->dy*lambda;
                    (*buffer)[3*i+2] = (*atms)[i]->z +(*atms)[i]->dz*lambda;
                }
            }
            xx = a1->x + lambda*a1->dx;
            yy = a1->y + lambda*a1->dy;
            zz = a1->z + lambda*a1->dz;
            /* now for the work */
#pragma _CNX no_recurrence 
            for( i=0;i< inbuffer; i++)
            {
                xt = xx - (*buffer)[3*i];
                yt = yy - (*buffer)[3*i+1];
                zt = zz - (*buffer)[3*i+2];
                r = xt*xt+yt*yt+zt*zt;
                if( r < 2.) r = 2.;
                r0 = sqrt(r); r = r*r*r ;
                /* the standard which follows is recursive */
                /*	 *V += 332.17752*a1->q*a2->q/r0;
                	*V -= a1->a*a2->a/r;
                	*V += a1->b*a2->b/r/r; 
                */
                /* use debye screen e(-r0/rdebye) */
                (*vector)[i] = a1->q*(*atms)[i]->q/r0*dielectric*exp(-r0/rdebye)
                               - a1->a*(*atms)[i]->a/r
                               + a1->b*(*atms)[i]->b/r/r;
            }
            for(i=0; i< inbuffer; i++)
                *V += (*vector)[i];

        } /* end of the inbuffer if check many lines ago */
    }
    free( atms); free( buffer);
    free( vector);
    return 1;

}
/* u_f_nonbon()
*
* u_f_nonbon increments the forces in the atom structures by the force
* due to the nonbon components.  NOTE THE WORD increment.
* the forces should first be zero'd.
* if not then this code will be invalid.  THIS IS DELIBERATE.
* on bigger (and better?) machines the different potential terms
* may be updated at random or in parrellel, if we assume that this routine
* will initialize the forces then we can't do this.
*/
int u_f_nonbon(lambda)
float lambda;
/*  returns 0 if error, 1 if OK */
{
    float r,r0,xt,yt,zt;
    float lcutoff,cutoff,get_f_variable();
    float rdebye;
    int inbond,inangle,i;
    ATOM *a1,*a2,*bonded[10],*angled[10];
    ATOM *a_next( ); /* returns first ATOM when called with -1 */
    ATOM *(*use)[];
    int uselist(),nuse,used;
    ATOM *cp,*bp;
    int a_number(),inbuffer;
    float (*buffer)[],xx,yy,zz,k;
    int invector,atomsused,ii,jj,imax;
    float (*vector)[];
    ATOM *(*atms)[];
    float dielectric;

    /* nonbonded potentials
    * do a double loop starting from the first atom to the 
    * last 
    * then from the second to the last 
    * etc
    *
    * also check to avoid bonded and 1-3 bonded atoms
    */
    cutoff = get_f_variable("cutoff");
    if( cutoff < 1.)
    {
        f_nonbon( lambda);
        return 1;
    }
    rdebye = cutoff/2.;
    dielectric = get_f_variable("dielec");
    if( dielectric < 1.) dielectric = 1.;
    dielectric = 332.17752/dielectric;
    if( !uselist( &use,&nuse,cutoff) )return 0;
    /*  get the number of atoms and allocate the memory for the array space */
    i = a_number();
    buffer = malloc( 3*i*sizeof(float) );
    if( buffer == NULL)
    { aaerror("cannot allocate memory in u_v_nonbon\n"); return 0;}
    vector = malloc( 3*i*sizeof(float) );
    if( vector == NULL)
    { aaerror("cannot allocate memory in u_v_nonbon\n"); return 0;}
    atms = malloc( i*sizeof(ATOM *) );
    if( atms == NULL)
    { aaerror("cannot allocate memory in u_v_nonbon\n"); return 0;}


    a1 = a_next(-1);
    a1 = a1->next;
    imax = a_number();
    used = 0;
    for( jj=1; jj<imax; jj++,a1=bp)
    {
        bp = a1->next;
        inbuffer = 0;
        if( (*use)[used] == a1)
        { used += 1;}
        else { aaerror("error in uselist - must abort"); return 0;}
        while( (*use)[used] != a1)
        {
            (*atms)[inbuffer++] = (*use)[used];
            used += 1;
        }
        used += 1;
        /* (*atms) now contains the list of atoms to be  done
        *  there are inbuffer of them
        *  of course inbuffer can be zero so we must check for that
        */
        if( inbuffer > 0)
        {
#pragma _CNX no_recurrence 
            for( i=0; i< inbuffer; i++)
            {
                (*buffer)[3*i  ] = (*atms)[i]->x;
                (*buffer)[3*i+1] = (*atms)[i]->y;
                (*buffer)[3*i+2] = (*atms)[i]->z;
            }
            if( lambda != 0.)
            {
#pragma _CNX no_recurrence 
                for( i=0; i< inbuffer; i++)
                {
                    (*buffer)[3*i  ] = (*atms)[i]->x +(*atms)[i]->dx*lambda;
                    (*buffer)[3*i+1] = (*atms)[i]->y +(*atms)[i]->dy*lambda;
                    (*buffer)[3*i+2] = (*atms)[i]->z +(*atms)[i]->dz*lambda;
                }
            }
            xx = a1->x + lambda*a1->dx;
            yy = a1->y + lambda*a1->dy;
            zz = a1->z + lambda*a1->dz;
            /* now for the work */
            for( i=0;i< inbuffer; i++)
            {
                xt = xx - (*buffer)[3*i];
                yt = yy - (*buffer)[3*i+1];
                zt = zz - (*buffer)[3*i+2];
                r = xt*xt+yt*yt+zt*zt;
                /* watch for FP errors*/
                if( r <= 1.)
                { r = 1.; }
                r0 = sqrt(r); xt = xt/r0; yt = yt/r0; zt = zt/r0;
                /* use debye screen e(-r0/rdebye) */
                /* d/dx(e(-r0/rdebye)/r0  = e(-r0/rdebye)*(-1/rdebye)/r0 + e(-r0/rdebye)/r) */
                k = -a1->q*(*atms)[i]->q*dielectric*exp(-r0/rdebye)*
                    (1./(rdebye*r0) +1./r) ;
                r = r*r*r;
                k += a1->a*(*atms)[i]->a/r/r0*6;
                k -= a1->b*(*atms)[i]->b/r/r/r0*12;
                (*vector)[3*i  ] = xt*k;
                (*vector)[3*i+1] = yt*k;
                (*vector)[3*i+2] = zt*k;
                /*
                *	a1->fx += ux*k; 
                *	a1->fy += uy*k; 
                *	a1->fz += uz*k; 
                *	a2->fx -= ux*k; 
                *	a2->fy -= uy*k; 
                *	a2->fz -= uz*k;
                */
            }
            for(i=0; i< inbuffer; i++)
            {
                a1->fx -= (*vector)[3*i  ];
                a1->fy -= (*vector)[3*i+1];
                a1->fz -= (*vector)[3*i+2];
                (*atms)[i] ->fx += (*vector)[3*i  ];
                (*atms)[i] ->fy += (*vector)[3*i+1];
                (*atms)[i] ->fz += (*vector)[3*i+2];
            }

        } /* end of the inbuffer if check many lines ago */
    }
    free( atms); free( buffer);
    free( vector);
    return 1;
}
/* uselist()
*  returns a pointer to an array of ATOM structure pointers
*  these are encoded as
*  a,bcdsfg,a where a is the outer most atom in the N^2 -N tree.
*
* not over brilliant
*
*  checks for change in atom number
*  other wise redoes the list every (nbstep (default= 10)) steps
*  makes a list with total storage (nbdeep (default = 20)) time a_number
*/
int uselist(  thelist,thesize, cutoff )
float cutoff;
int *thesize;
ATOM *(**thelist)[];
{
    /* static stuff used to keep information about status */
    static int  oldatomnumber = 0;
    static int  since = 0,lsize;
    static ATOM *(*local)[];
    static float oldcutoff = -1;
    int a_number();
    ATOM *a_next(),*a1,*a2,*ap,*bp;
    int i,j,k,max;
    int get_i_variable(),set_i_variable();
    float lcutoff;
    float x,y,z,r,rcut;

    /* check on wether to redo it or not  */
    i = a_number();
    j = get_i_variable("nbstep");
    if( j <= 0) j= 10;
    if( (i == oldatomnumber) && (since < j) && (cutoff == oldcutoff) )
    {
        *thelist = local;
        *thesize = lsize;
        since += 1;
        return 1;
    }
    /* a free and malloc are used because nbdeep may change */
RESET:
    /* don't free if it hasn't been malloc'd */
    if( oldatomnumber > 0) free(local);
    oldcutoff = cutoff;
    lcutoff = -cutoff;
    since = 0;
    oldatomnumber = i;
    j = get_i_variable("nbdeep");
    if( j<= 0) j = 20;
    max = i*j;
    local =  malloc( max*sizeof(ATOM * ) );
    if( local == NULL )
    { aaerror("cannot allocate uselist memory"); exit(0);}
    /* now have the uselist allocated */
    *thelist = local;

    *thesize = 0;
    rcut = cutoff*cutoff;
    a1 = a_next(-1);
    a1 = a1->next;
    for( i=1; i< oldatomnumber; i++,a1=ap)
    {
        ap = a1->next;
        (*local)[*thesize] = a1;
        *thesize += 1;
        a2 = a_next(-1);
        for( j=0;j<i; j++,a2=bp)
        {
            for(k=0; k< a1->dontuse; k++)
            {
                if( a2 == a1->excluded[k]) goto SKIP;
            }
            if( (a2->x-a1->x) > cutoff) goto SKIP;
            if( (a2->x-a1->x) < lcutoff) goto SKIP;
            if( (a2->y-a1->y) > cutoff) goto SKIP;
            if( (a2->y-a1->y) < lcutoff) goto SKIP;
            if( (a2->z-a1->z) > cutoff) goto SKIP;
            if( (a2->z-a1->z) < lcutoff) goto SKIP;
            /* now calculate the radius */
            x = a2->x -a1->x;
            y = a2->y -a1->y;
            z = a2->z -a1->z;
            r = x*x + y*y + z*z;
            if( r > rcut) goto SKIP;

            (*local)[*thesize] = a2;
            *thesize += 1;
            if( *thesize >= max)
            {aaerror("please increase nbdeep (seti nbdeep (>20);)");
                i = a_number();
                j = get_i_variable("nbdeep");
                if( j== 0) j = 20;
                if( j == i+2)
                {aaerror("Terrible error in uselist, too many interactions");
                    exit( 0) ; }
                j = 2*j; if( j > i+2) j = i+2;
                set_i_variable("nbdeep",j);
                goto RESET;
            }
SKIP:
            bp = a_next(1);
        }
        (*local)[*thesize] = a1;
        *thesize += 1;
        lsize = *thesize;

    }
    /*printf(" uselist finished %d %d\n",*thesize,max); */
    return 1;
}




