/* bfgs.c
*
*  bfgs optimizer for AMMP
*  this is a bit of a memory hog when running
*  maximum size depends on the memory available on the
*  machine
*
*
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
#include <string.h>
#ifdef ANSI
#include <stdlib.h>
#endif
#include "ammp.h"


int bfgs( vfs,ffs,nfs, nstep, toler )
int nfs,(*vfs[])(),(*ffs[])();
int nstep;
float toler;

{
    int a_number(),i,j,n,na;
    int istep,iatm,cent;
    float  (*H)[];  /* space to be malloced for the matrix */
    float  (*Hup)[];  /* space to be malloced for the matrix */
    ATOM *(*atms)[], *a_next(); /* a_next returns the FIRST atom when called with -1*/
    int cngdel(); /* conjugate gradients -used as a fallback */
    float linmin(); /* line minimizer in optimist.c */
    float dg,gh,hg,ghg,ddtconst,alpha ;
    float a_l2_f();
    float (*gamma)[], (*delta)[];
    int a_f_zero();
    ATOM *ap,*bp;

    na = a_number();
    n = na*3;

    atms = malloc( na*sizeof( ATOM *) );
    if( atms == NULL)
    { aaerror(" cannot allocate memory in bfgs\n");
        return 0; /* if cannot allocate for atms then not worth continuing */
    }
    gamma = malloc( n*sizeof(float));
    if( gamma == NULL )
    { aaerror(" cannont allocate memory in bfgs\n");
        aaerror(" using CNGDEL instead - sorry \n");
        cngdel( vfs,ffs,nfs,nstep,nstep,toler,1);
        return 1;
    }
    delta = malloc( n*sizeof(float));
    if( delta == NULL )
    { aaerror(" cannont allocate memory in bfgs\n");
        aaerror(" using CNGDEL instead - sorry \n");
        cngdel( vfs,ffs,nfs,nstep,nstep,toler,1);
        return 1;
    }
    H = malloc( n*n*sizeof(float));
    if( H == NULL )
    { aaerror(" cannont allocate memory in bfgs\n");
        aaerror(" using CNGDEL instead - sorry \n");
        cngdel( vfs,ffs,nfs,nstep,nstep,toler,1);
        return 1;
    }
    Hup = malloc( n*n*sizeof(float));
    if( Hup == NULL )
    { aaerror(" cannont allocate memory in bfgs\n");
        aaerror(" using CNGDEL instead - sorry \n");
        cngdel( vfs,ffs,nfs,nstep,nstep,toler,1);
        return 1;
    }
    /* initialize the atms array */
    j = -1;
    for( i=0; i< na; i++)
    { (*atms)[i] = a_next(j); j = 1;}

    /* initialize H to I (scaled by gradient squared) */
    for( i=0; i< n*n; i++)
        (*H)[i] = 0.;
    /* and set up the gradient */
    /* the forces are the negative gradient */
    a_f_zero();
    for( i=0; i< nfs; i++)
        (*ffs[i])(0.);
    for( i=0; i< na; i++)
    {
        (*atms)[i]->fx *= -1;
        (*atms)[i]->fy *= -1;
        (*atms)[i]->fz *= -1;
    }
    gh = a_l2_f();
    gh = gh/n/2+1;
    for( i=0; i<n; i++)
        (*H)[n*i +i ] = 1./gh;
    /*
    * and now to do the work  
    */
    for( istep = 0; istep< nstep; istep++)
    {
        dg = 0.; gh = 0.;
        a_d_zero();
        /* loop over all of the atoms */
        for(iatm = 0; iatm<na; iatm++)
        {
            ap = (*atms)[iatm];
            (*gamma)[iatm*3] = -ap->fx;
            (*gamma)[iatm*3+1] = -ap->fy;
            (*gamma)[iatm*3+2] = -ap->fz;
            for( j=0; j< na; j++)
            {
                cent = 3*iatm*n +3*j ;
                ap->dx -=
                    (*H)[cent]* (*atms)[j]->fx +
                    (*H)[cent+1]* (*atms)[j]->fy +
                    (*H)[cent+2]* (*atms)[j]->fz;
                ap->dy -=
                    (*H)[cent+n]* (*atms)[j]->fx +
                    (*H)[cent+n+1]* (*atms)[j]->fy +
                    (*H)[cent+n+2]* (*atms)[j]->fz;
                ap->dz -=
                    (*H)[cent+n+n]* (*atms)[j]->fx +
                    (*H)[cent+n+n+1]* (*atms)[j]->fy +
                    (*H)[cent+n+n+2]* (*atms)[j]->fz;
            }
        }
        /* do the line search */
        gh = a_l2_f()/na;
        alpha = linmin( vfs,nfs, 4.);
        printf(" bfgs %d %f %f\n",istep, gh,alpha);
        if( gh < toler*toler) goto CLEANUP;
        if( gh < 1.e-7) goto CLEANUP;
        if( alpha < 1.e-10) goto CLEANUP;
        a_inc_d( alpha );
        a_f_zero();
        for(i=0;i<nfs; i++)
            (*ffs[i])(0.);
        for( i=0; i< na; i++)
        {
            (*atms)[i]->fx *= -1;
            (*atms)[i]->fy *= -1;
            (*atms)[i]->fz *= -1;
        }

        /* make up delta grad */
        for( iatm = 0; iatm< na; iatm++)
        {
            ap = (*atms)[iatm];
            (*gamma)[iatm*3] += ap->fx;
            (*gamma)[iatm*3+1] += ap->fy;
            (*gamma)[iatm*3+2] += ap->fz;
            (*delta)[iatm*3] = alpha*ap->dx;
            (*delta)[iatm*3+1] = alpha*ap->dy;
            (*delta)[iatm*3+2] = alpha*ap->dz;
        }
        dg = 0.; gh = 0.; ghg = 0.; hg = 0;
        ddtconst = 0.;
        for( i=0; i< n*n; i++)
            (*Hup)[i] = 0;
        for(i=0;i<n;i++)
        {
            dg+= (*delta)[i]*(*gamma)[i];
            hg = 0;
            for(j=0;j<n;j++)
                hg+= (*H)[i*n+j]*(*gamma)[j];
            ghg += (*gamma)[i]*hg;
        }
        if( dg < 1.e-6*ghg) goto CLEANUP;
        ddtconst = (1.+ghg/dg)/dg;
        for(i=0; i<n; i++)
        {
            for(j=0;j<n;j++)
                (*Hup)[i*n+j] += (*delta)[i]*(*delta)[j]*ddtconst;
        }
        for( i=0; i<n; i++)
        {
            gh = 0.;
            for(j=0; j<n; j++)
                gh += (*gamma)[j]*(*H)[j*n+i];
            gh = -gh/dg;
            for( j=0; j<n; j++)
                (*Hup)[j*n+i] += (*delta)[j]*gh;
        }
        for( i=0; i<n; i++)
        {
            hg = 0.;
            for(j=0; j<n; j++)
                hg += (*gamma)[j]*(*H)[i*n+j];
            hg = -hg/dg;
            for( j=0; j<n; j++)
                (*Hup)[i*n+j] += (*delta)[j]*hg;
        }
    }
CLEANUP:
    free( H );
    free(Hup);
    free(delta); free(gamma); free(atms);

    return 1;
} /*end of routine */
