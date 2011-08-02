/* vshadow.c
*
* collection of routines to service shadowded potentials
*
* shadowed potentials have x,y,z,w as coords
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
*  copyright 1992, 1993 Robert W. Harrison
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
* but otherwise is self-contained. Note the hooks for Non-shadowed potentials
*/

/* v_shadow()
* this function sums up the potentials
* for the atoms defined in the shadow data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int fv_update_shadow(  lambda )
float lambda;
{
    float r,r0,xt,yt,zt,wt;
    float k,k1,k2,k3,k4,k5;
    float ka3,ka2;
    float kb3,kb2;
    float get_f_variable();
    int inbond,inangle,i;
    ATOM *a1,*a2,*bonded[10],*angled[10];
    ATOM *a_next( ); /* returns first ATOM when called with -1 */
    int (*indexes)[],inindex,in;
    int a_number();
    int ii,j,jj,imax,inclose;
    float (*vector)[];
    /*	float (*vecold)[];
    */
    ATOM *close[NCLOSE],*(*atomall)[];
    float mxdq,dielectric,mxcut;
    static float dielecold = -1.;


    mxdq = get_f_variable("mxdq");
    /*	if( mxdq < 0.) mxdq = 0.;
    */
    mxcut = get_f_variable("mxcut");
    if( mxcut < 0.) mxcut= 5.;
    if( mxdq > 0.) mxdq = mxdq*mxdq;

    dielectric = get_f_variable("dielec");
    if( dielectric < 1.) dielectric = 1.;
    if( dielecold != dielectric)
    {
        dielecold = dielectric;
        mxdq = -1.;
    }
    dielectric = 332.17752/dielectric;

    /*  get the number of atoms and allocate the memory for the array space */
    i = a_number();
    vector = malloc( 5*i*sizeof(float) );
    if( vector == NULL)
    { aaerror("cannot allocate memory in v_shadow\n"); return 0;}

    atomall = malloc( i*sizeof(ATOM *) );
    if( atomall == NULL)
    {aaerror("cannot allocate memory in v_shadow\n"); return 0;}

    imax = a_number();
    for( i=0; i< imax; i++)
    {
        (*atomall)[i] = a_next(i);
    }
    /* first check if anyone's moved and update the lists */
    /* note that this must be a look-ahead rather than
    *  look back search because
    * we cannot update ->px until we've used that atom !!! */
    for( ii=0; ii< imax; ii++)
    {
        a1 = (*atomall)[ii];
        xt = a1->dx*lambda +a1->x - a1->px;
        yt = a1->dy*lambda +a1->y - a1->py;
        zt = a1->dz*lambda +a1->z - a1->pz;
        wt = a1->dw*lambda +a1->w - a1->pw;
        r = xt*xt + yt*yt + zt*zt + wt*wt;
        if( r > mxdq ) goto DOIT;
    }

    free( vector);
    /*	free( vecold);
    */
    free (atomall);
    return 1;
DOIT:
    /* no MMBOND with shadow because we dont know if its any good
    	xt = get_f_variable("mmbox");
    	if( xt > 0.)
    	{ free(vector); free(atomall); mm_fv_update_shadow(lambda); return 1;}
    */
    indexes = malloc( imax* sizeof(int) );
    if( indexes == NULL ){
        aaerror(" cannot allocate memory in fv_update\n");
        return 0;}
    for( ii=0; ii< imax; ii++)
    {
        a1 = (*atomall)[ii];
        a1 -> VP = 0.;
        a1 -> dpx = 0.;
        a1 -> dpy = 0.;
        a1 -> dpz = 0.;
        a1 -> dpw = 0.;
        a1 -> qxx = 0.;
        a1 -> qxy = 0.;
        a1 -> qxz = 0.;
        a1 -> qyy = 0.;
        a1 -> qyz = 0.;
        a1 -> qzz = 0.;
        a1 -> qxw = 0.;
        a1 -> qyw = 0.;
        a1 -> qzw = 0.;
        a1 -> qww = 0.;
#ifdef CUBIC
        a1 -> qxxx = 0.;
        a1 -> qxxy = 0.;
        a1 -> qxxz = 0.;
        a1 -> qxyy = 0.;
        a1 -> qxyz = 0.;
        a1 -> qxzz = 0.;
        a1 -> qyyy = 0.;
        a1 -> qyyz = 0.;
        a1 -> qyzz = 0.;
        a1 -> qzzz = 0.;
#endif
        for( j=0; j< NCLOSE; j++)
            a1->close[j] = NULL;

    }
    for( ii=0; ii<  imax; ii++)
    { /* if this is met we update the expansion for this atom */
        a1 = (*atomall)[ii];
        inclose = 0;
        if( lambda != 0.)
        {
#pragma _CNX no_recurrence
            for( i=ii+1; i< imax; i++)
            {
                a2 = (*atomall)[i];
                j = i*5;
                (*vector)[j  ] = a2->x - a1->x + lambda*(a2->dx -a1->dx);
                (*vector)[j+1] = a2->y - a1->y + lambda*(a2->dy -a1->dy);
                (*vector)[j+2] = a2->z - a1->z + lambda*(a2->dz -a1->dz);
                (*vector)[j+3] = a2->w - a1->w + lambda*(a2->dw -a1->dw);
            }
        }else {
#pragma _CNX no_recurrence
            for( i=ii+1; i< imax; i++)
            {
                a2 = (*atomall)[i];
                j = i*5;
                (*vector)[j  ] = a2->x - a1->x ;
                (*vector)[j+1] = a2->y - a1->y ;
                (*vector)[j+2] = a2->z - a1->z ;
                (*vector)[j+3] = a2->w - a1->w ;
            }
        } /* end of difference position into vector loops */
#pragma _CNX no_recurrence
        for( i=ii+1; i< imax; i++)
        {
            j = i*5;
            (*vector)[j+4] = sqrt((*vector)[j]*(*vector)[j] +
                                  (*vector)[j+1]*(*vector)[j+1] +
                                  (*vector)[j+2]*(*vector)[j+2] +
                                  (*vector)[j+3]*(*vector)[j+3]);
        }
        /* add the new components */
        /* first extract indexes */
        inindex = 0;
        for( i=ii+1; i< imax; i++)
        {
            a2 = (*atomall)[i];
            for( j=0; j< a1->dontuse; j++)
            { if( a2 == a1->excluded[j]) goto SKIPNEW;}
            j = i*5;
            if( (*vector)[j+4] > mxcut || inclose > NCLOSE)
            {
                (*indexes)[inindex++] = i;
            }else {
                a1->close[inclose++] = (*atomall)[i];
            }
            if( inclose == NCLOSE)
            {
                aaerror(
                    " fv_update_shadow> too many atoms increase NCLOSE or decrease mxcut");
            }
SKIPNEW:   i = i;
        }
        /* and then use them */
#pragma _CNX no_recurrence
        for( in=0; in< inindex; in++)
        {
            i = (*indexes)[in];
            a2 = (*atomall)[i];
            j = i*5;
            r0 = (*vector)[j+4];
            r = r0*r0;
            r = r*r*r; /* r0^6 */
            xt = a1->q*a2->q*dielectric/r0;
            yt = a1->a*a2->a/r;
            zt = a1->b*a2->b/r/r;
            k = xt - yt + zt;
            xt = xt/r0; yt = yt/r0; zt = zt/r0;
            k1 = xt - yt*6 + zt*12;
            xt = xt/r0; yt = yt/r0; zt = zt/r0;
            k2 = xt*3; ka2 = - yt*6*8; kb2 =  zt*12*14;
#ifdef CUBIC
            xt = xt/r0; yt = yt/r0; zt = zt/r0;
            k3 = -xt*5*3; ka3 =   yt*6*8*10 ; kb3 =  -zt*12*14*16;
#endif
            k1 = -k1;
            xt = (*vector)[j]/r0 ;
            yt = (*vector)[j+1]/r0 ;
            zt = (*vector)[j+2] /r0;
            wt = (*vector)[j+3]/r0;

            /*
            xt = (*vector)[j] ;
            yt = (*vector)[j+1] ;
            zt = (*vector)[j+2] ;
            */
            a1->VP += k;
            a2->dpx -= k1*xt;
            a1->dpx += k1*xt;
            a2->dpy -= k1*yt;
            a1->dpy += k1*yt;
            a2->dpz -= k1*zt;
            a1->dpz += k1*zt;
            a2->dpw -= k1*wt;
            a1->dpw += k1*wt;
            a2->qxx -= k2*(xt*xt - 1./3) + ka2*(xt*xt - 1./8) + kb2*(xt*xt-1./14) ;
            a1->qxx -= k2*(xt*xt - 1./3) + ka2*(xt*xt - 1./8) + kb2*(xt*xt-1./14) ;
            a2->qxy -= (k2+ka2+kb2)*yt*xt;
            a1->qxy -= (k2+ka2+kb2)*yt*xt;
            a2->qxz -= (k2+ka2+kb2)*zt*xt;
            a1->qxz -= (k2+ka2+kb2)*zt*xt;
            a2->qxw -= (k2+ka2+kb2)*wt*xt;
            a1->qxw -= (k2+ka2+kb2)*wt*xt;
            a2->qyy -= k2*(yt*yt - 1./3) + ka2*(yt*yt - 1./8) + kb2*(yt*yt-1./14) ;
            a1->qyy -= k2*(yt*yt - 1./3) + ka2*(yt*yt - 1./8) + kb2*(yt*yt-1./14) ;
            a2->qyz -= (k2+ka2+kb2)*yt*zt;
            a1->qyz -= (k2+ka2+kb2)*yt*zt;
            a2->qyw -= (k2+ka2+kb2)*yt*wt;
            a1->qyw -= (k2+ka2+kb2)*yt*wt;
            a2->qzz -= k2*(zt*zt - 1./3) + ka2*(zt*zt - 1./8) + kb2*(zt*zt-1./14) ;
            a1->qzz -= k2*(zt*zt - 1./3) + ka2*(zt*zt - 1./8) + kb2*(zt*zt-1./14) ;
            a2->qzw -= (k2+ka2+kb2)*wt*zt;
            a1->qzw -= (k2+ka2+kb2)*wt*zt;
            a2->qww -= k2*(wt*wt - 1./3) + ka2*(wt*wt - 1./8) + kb2*(wt*wt-1./14) ;
            a1->qww -= k2*(wt*wt - 1./3) + ka2*(wt*wt - 1./8) + kb2*(wt*wt-1./14) ;
#ifdef CUBIC
            a2->qxxx -= k3*(xt*xt*xt - xt*( 9./15 )) ;
            a2->qxxx -= ka3*(xt*xt*xt - xt*( 24./80 )) ;
            a2->qxxx -= kb3*(xt*xt*xt - xt*( 42./(14*16)));
            a1->qxxx += k3*(xt*xt*xt - xt*( 9./15 )) ;
            a1->qxxx += ka3*(xt*xt*xt - xt*( 24./80 )) ;
            a1->qxxx += kb3*(xt*xt*xt - xt*( 42./(14*16)));
            a2->qxxy -= k3*(yt*xt*xt - yt*( 6./ 15));
            a2->qxxy -= ka3*(yt*xt*xt - yt*( 11./ 80));
            a2->qxxy -= kb3*(yt*xt*xt - yt*( 17./ (14*16)));
            a1->qxxy += k3*(yt*xt*xt - yt*( 6./ 15));
            a1->qxxy += ka3*(yt*xt*xt - yt*( 11./ 80));
            a1->qxxy += kb3*(yt*xt*xt - yt*( 17./ (14*16)));
            a2->qxxz -= k3*(zt*xt*xt - zt*( 6./ 15));
            a2->qxxz -= ka3*(zt*xt*xt - zt*( 11./ 80));
            a2->qxxz -= kb3*(zt*xt*xt - zt*( 17./ (14*16)));
            a1->qxxz += k3*(zt*xt*xt - zt*( 6./ 15));
            a1->qxxz += ka3*(zt*xt*xt - zt*( 11./ 80));
            a1->qxxz += kb3*(zt*xt*xt - zt*( 17./ (14*16)));
            a2->qxyy -= k3*(yt*yt*xt - xt*( 6./ 15));
            a2->qxyy -= ka3*(yt*yt*xt - xt*( 11./ 80));
            a2->qxyy -= kb3*(yt*yt*xt - xt*( 17./ (14*16)));
            a1->qxyy += k3*(yt*yt*xt - xt*( 6./ 15));
            a1->qxyy += ka3*(yt*yt*xt - xt*( 11./ 80));
            a1->qxyy += kb3*(yt*yt*xt - xt*( 17./ (14*16)));
            a2->qxyz -= (k3+ka3+kb3)*yt*zt*xt;
            a1->qxyz += (k3+ka3+kb3)*yt*zt*xt;
            a2->qxzz -= k3*(zt*zt*xt - xt*( 6./ 15));
            a2->qxzz -= ka3*(zt*zt*xt - xt*( 11./ 80));
            a2->qxzz -= kb3*(zt*zt*xt - xt*( 17./ (14*16)));
            a1->qxzz += k3*(zt*zt*xt - xt*( 6./ 15));
            a1->qxzz += ka3*(zt*zt*xt - xt*( 11./ 80));
            a1->qxzz += kb3*(zt*zt*xt - xt*( 17./ (14*16)));
            a2->qyyy -= k3*(yt*yt*yt - yt*( 9./15 )) ;
            a2->qyyy -= ka3*(yt*yt*yt - yt*( 24./80 )) ;
            a2->qyyy -= kb3*(yt*yt*yt - yt*( 42./(14*16)));
            a1->qyyy += k3*(yt*yt*yt - yt*( 9./15 )) ;
            a1->qyyy += ka3*(yt*yt*yt - yt*( 24./80 )) ;
            a1->qyyy += kb3*(yt*yt*yt - yt*( 42./(14*16)));
            a2->qyyz -= k3*(yt*yt*zt - zt*( 6./ 15));
            a2->qyyz -= ka3*(yt*yt*zt - zt*( 11./ 80));
            a2->qyyz -= kb3*(yt*yt*zt - zt*( 17./ (14*16)));
            a1->qyyz += k3*(yt*yt*zt - zt*( 6./ 15));
            a1->qyyz += ka3*(yt*yt*zt - zt*( 11./ 80));
            a1->qyyz += kb3*(yt*yt*zt - zt*( 17./ (14*16)));
            a2->qyzz -= k3*(zt*zt*yt - yt*( 6./ 15));
            a2->qyzz -= ka3*(zt*zt*yt - yt*( 11./ 80));
            a2->qyzz -= kb3*(zt*zt*yt - yt*( 17./ (14*16)));
            a1->qyzz += k3*(zt*zt*yt - yt*( 6./ 15));
            a1->qyzz += ka3*(zt*zt*yt - yt*( 11./ 80));
            a1->qyzz += kb3*(zt*zt*yt - yt*( 17./ (14*16)));
            a2->qzzz -= k3*(zt*zt*zt - zt*( 9./15 )) ;
            a2->qzzz -= ka3*(zt*zt*zt - zt*( 24./80 )) ;
            a2->qzzz -= kb3*(zt*zt*zt - zt*( 42./(14*16)));
            a1->qzzz += k3*(zt*zt*zt - zt*( 9./15 )) ;
            a1->qzzz += ka3*(zt*zt*zt - zt*( 24./80 )) ;
            a1->qzzz += kb3*(zt*zt*zt - zt*( 42./(14*16)));
#endif

            /* debugging
            	j = i *4;
            	fprintf(stderr," mxcut %f %f inclose %d who %d \n",mxcut,(*vector)[j+3],inclose,(*atomall)[i]->serial);
            	fprintf(stderr," vector %f %f %f \n", (*vector)[j],(*vector)[j+1],(*vector)[j+2]);
            */
        }/* end of loop i */
        /* merge the non-bond mxcut lists */
        a1->close[inclose] == NULL;
        /* set the position */
        a1->px = a1->dx*lambda + a1->x;
        a1->py = a1->dy*lambda + a1->y;
        a1->pz = a1->dz*lambda + a1->z;
        a1->pw = a1->dw*lambda + a1->w;

    }  /* end of ii loop */

    a_inactive_f_zero();

    free( indexes);
    free( vector);
    /*	free( vecold);
    */
    free (atomall);
    return 1;

}


/* f_shadow()
*
* f_shadow increments the forces in the atom structures by the force
* due to the shadow components.  NOTE THE WORD increment.
* the forces should first be zero'd.
* if not then this code will be invalid.  THIS IS DELIBERATE.
* on bigger (and better?) machines the different potential terms
* may be updated at random or in parrellel, if we assume that this routine
* will initialize the forces then we can't do this.
*/
int f_shadow(lambda)
float lambda;
/*  returns 0 if error, 1 if OK */
{
    float ux,uy,uz,uw;
    float k,r,r0,xt,yt,zt,wt;
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
    float fx,fy,fz,fw;
    float xt2,xt3,xt4;
    float yt2,yt3,yt4;
    float zt2,zt3,zt4;

    /* first update the lists
    *  this routine checks if any atom has
    *   broken the mxdq barrier and updates the
    * forces, potentials and expansions thereof */
    fv_update_shadow( lambda);

    dielectric = get_f_variable("dielec");
    if( dielectric < 1.) dielectric = 1.;
    dielectric = 332.17752/dielectric;

    /*  get the number of atoms and allocate the memory for the array space */
    i = a_number();
    atomall = malloc( i*sizeof(ATOM *) );
    if( atomall == NULL)
    {aaerror("cannot allocate memory in f_shadow"); return 0;}

    imax = a_number();
    for( i=0; i< imax; i++)
    {
        (*atomall)[i] = a_next(i);
    }
    for( i= 0; i< imax; i++)
    {
        fx = 0.; fy = 0.; fz = 0.; fw = 0.;
        a1 = (*atomall)[i];
        xt = a1->dx*lambda +a1->x - a1->px;
        yt = a1->dy*lambda +a1->y - a1->py;
        zt = a1->dz*lambda +a1->z - a1->pz;
        wt = a1->dw*lambda +a1->w - a1->pw;


        fx = (a1->qxx*xt + a1->qxy*yt
              + a1->qxz*zt + a1->qxw*wt) ;
        fy = (a1->qxy*xt + a1->qyy*yt
              + a1->qyz*zt + a1->qyw*wt) ;
        fz = (a1->qxz*xt + a1->qyz*yt
              + a1->qzz*zt + a1->qzw*wt) ;
        fw = (a1->qxw*xt + a1->qyw*yt
              + a1->qzw*zt + a1->qww*wt) ;
#ifdef CUBIC
        xt2 = xt*xt; yt2 = yt*yt; zt2 = zt*zt;
        fx += a1->qxxx*xt2/2. + a1->qxxy*xt*yt + a1->qxxz*xt*zt
              + a1->qxyy*yt/2. + a1->qxyz*yt*zt + a1->qxzz*zt2/2.;
        fy += a1->qxxy*xt2/2. + a1->qxyy*xt*yt + a1->qxyz*xt*zt
              + a1->qyyy*yt/2. + a1->qyyz*yt*zt + a1->qyzz*zt2/2.;
        fz += a1->qxxz*xt2/2. + a1->qxyz*xt*yt + a1->qxzz*xt*zt
              + a1->qyyz*yt/2. + a1->qyzz*yt*zt + a1->qzzz*zt2/2.;
#endif
        a1->fx += fx  + a1->dpx;
        a1->fy += fy  + a1->dpy;
        a1->fz += fz  + a1->dpz;
        a1->fw += fw  + a1->dpw;

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
            uw = (a2->dw -a1->dw)*lambda +(a2->w -a1->w);
            r = ux*ux + uy*uy + uz*uz +uw*uw; r0 = sqrt(r);
            ux = ux/r0; uy = uy/r0; uz = uz/r0; uw = uw/r0;
            k = -dielectric*a1->q*a2->q/r;
            r = r*r*r;
            k += a1->a*a2->a/r/r0*6;
            k -= a1->b*a2->b/r/r/r0*12;
            a1->fx += ux*k;
            a1->fy += uy*k;
            a1->fz += uz*k;
            a1->fw += uw*k;
            a2->fx -= ux*k;
            a2->fy -= uy*k;
            a2->fz -= uz*k;
            a2->fw -= uw*k;
        }
    }

    a_inactive_f_zero();
    free( atomall); return 1;

}
/* v_shadow()
* this function sums up the potentials
* for the atoms defined in the shadow data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int v_shadow( V, lambda )
float *V,lambda;
{
    float r,r0,xt,yt,zt,wt;
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

    fv_update_shadow( lambda);

    dielectric = get_f_variable("dielec");
    if( dielectric < 1.) dielectric = 1.;
    dielectric = 332.17752/dielectric;

    /*  get the number of atoms and allocate the memory for the array space */
    i = a_number();
    atomall = malloc( i*sizeof(ATOM *) );
    if( atomall == NULL)
    {aaerror("cannot allocate memory in v_shadow"); return 0;}

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
        wt = a1->dw*lambda +a1->w - a1->pw;
        vx -= (a1->dpx*xt + a1->dpy*yt
               + a1->dpz*zt + a1->dpw*wt) ;
        vx -= (xt*(.5*a1->qxx*xt + a1->qxy*yt + a1->qxz*zt +a1->qxw*wt)
               + yt*(.5*a1->qyy*yt + a1->qyz*zt + a1->qyw*wt) +
               zt*(.5*a1->qzz*zt+a1->qzw*wt) + .5*wt*wt*a1->qww );
        /* do the close atoms */
        for( jj=0; jj< NCLOSE; jj++)
        { if( a1->close[jj] == NULL) break; }
        for( ii=0; ii< jj;ii++)
        {
            a2 = a1->close[ii];
            xt = (a2->dx -a1->dx)*lambda +(a2->x -a1->x);
            yt = (a2->dy -a1->dy)*lambda +(a2->y -a1->y);
            zt = (a2->dz -a1->dz)*lambda +(a2->z -a1->z);
            wt = (a2->dw -a1->dw)*lambda +(a2->w -a1->w);
            r = xt*xt + yt*yt + zt*zt + wt*wt; r0 = sqrt(r);
            /*      xt = xt/r0; yt = yt/r0; zt = zt/r0;
            */
            k = dielectric*a1->q*a2->q/r0;
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

/*zone_shadow()
* this function sums up the potentials
* for the atoms defined in the shadow data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int zone_shadow( V, lambda ,alist, inalist )
float *V,lambda;
ATOM *( *alist)[] ;
int inalist;
{
    float r,r0,xt,yt,zt,wt;
    float lcutoff,cutoff,get_f_variable();
    int inbond,inangle,i,ii;
    ATOM *a1,*a2;
    ATOM *a_next( ); /* returns first ATOM when called with -1 */
    float dielectric,ve,va,vh;
    ATOM *a_m_serial();

    /* shadowded potentials
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
            wt = (a1->w -a2->w) + lambda*(a1->dw - a2->dw);
            if( (wt > cutoff) || (wt < lcutoff) ) goto SKIP;

            r = xt*xt+yt*yt+zt*zt+wt*wt;
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
/* routine to initialize w
*/
int init_fourd(howmuch )
float howmuch;
{
    ATOM *ap,*a_next();
    int i,numatm,a_number();
    float randf();

    if( howmuch < 1.e-3) howmuch =1.;
    numatm = a_number();
    if( numatm <= 0) return;
    for(i =0; i< numatm; i++)
    {
        ap = a_next(i);
        if( ap->active)
            ap->w = howmuch*2*randf()-1.;
    }
}
/* routines to force w to zero
* the idea is that after generating a 4-d structure you make it 3-d
* by slowly restraining w to zero
*/

int f_fourd( lambda)
float lambda;
{
    float kfourd;
    float get_f_variable();
    ATOM *ap,*a_next();
    int i,numatm,a_number();

    kfourd = get_f_variable("kfourd");
    if( kfourd <= 0) return;

    numatm = a_number();
    if( numatm <= 0 ) return;

    for( i=0; i< numatm; i++)
    {
        ap = a_next(i);
        ap->fw -= kfourd*(ap->w+lambda*ap->dw);
    }
}

int v_fourd( V,lambda)
float lambda,*V;
{
    float kfourd,wt;
    float get_f_variable();
    ATOM *ap,*a_next();
    int i,numatm,a_number();

    kfourd = get_f_variable("kfourd");
    if( kfourd <= 0) return;

    numatm = a_number();
    if( numatm <= 0 ) return;

    for( i=0; i< numatm; i++)
    {
        ap = a_next(i);
        wt = ap->w + lambda*ap->dw;
        *V += .5*kfourd*wt*wt;
    }
}
