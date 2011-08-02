/* doubletime.c
*
*  numerical perturbation theory
*
*  per atom 
*  X(t) = x(0) +  xr cos( at + bsin(gt )) + xi sin(at + bsin(gt) )
*
*  pac estimates X better than V so use X(t)
*
*  perform a nonlinear fit then extrapolate
*
*
*/

/*
*  copyright 1993,1994,1995 Robert W. Harrison
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

#define NPREDICT 6

int doubletime( forces,nforces,nstep,dlong,dshort,temper)
int (*forces[])(),nforces;
int nstep;
float dlong,dshort,temper;
{
    float (*positions)[];
    float (*velos)[];
    float (*theforces)[];
    int i,numatm,a_number();
    int istep;
    int pac();
    int nlocal,ilocal;
    ATOM *ap,*a_next();
    float steps[2*NPREDICT],coef[NPREDICT];
    float dtime;


    numatm = a_number();
    if( numatm < 1) return 0;
    if( nstep < 1) return 0;
    if( dshort == 0. ) dshort = 0.00001;
    if( dlong == 0.) dlong = 0.001;
    dtime = dshort;
    nlocal = 2*NPREDICT;
    nlocal = NPREDICT;
    /* predict in cycles */
    if( dlong > 1.) dlong = NPREDICT*dlong*dshort;
    if( temper < 1. ) temper = 300.;


    positions = malloc( 3*nlocal*numatm*sizeof( float ));
    /*
    velos = malloc( 3*nlocal*numatm*sizeof( float ));
    theforces = malloc( 3*nlocal*numatm*sizeof( float ));
    */
    if( positions == NULL ){
        aaerror("cannot allocate memory in doubletime\n");
        return 1;
    }
    for( istep = 0; istep < nstep; istep++)
    {

        /* predict nlocal dshort steps */
        /* actually use current position and nlocal -1 updates */
        v_rescale (temper);
        for( ilocal=0; ilocal< nlocal ; ilocal++)
        {
            for( i=0; i< numatm; i++)
            {
                ap = a_next(i);
                (*positions)[3*(i*nlocal + ilocal)  ] = ap->x;
                (*positions)[3*(i*nlocal + ilocal)+1] = ap->y;
                (*positions)[3*(i*nlocal + ilocal)+2] = ap->z;
                /*
                (*velos)[3*(i*nlocal + ilocal)  ] = ap->vx;
                (*velos)[3*(i*nlocal + ilocal)+1] = ap->vy;
                (*velos)[3*(i*nlocal + ilocal)+2] = ap->vz;
                (*theforces)[3*(i*nlocal + ilocal)  ] = ap->fx;
                (*theforces)[3*(i*nlocal + ilocal)+1] = ap->fy;
                (*theforces)[3*(i*nlocal + ilocal)+2] = ap->fz;
                */
            }
            if( ilocal != nlocal-1)
                pac( forces,nforces,1,dshort);
        } /* ilocal */

        /* now fit and predict */
        for( i = 0; i< numatm; i++)
        {
            ap = a_next( i );
            /* x */
            for( ilocal = 0; ilocal < nlocal ; ilocal++)
            {
                steps[ilocal] = (*positions)[3*(i*nlocal + ilocal)];
            }/* ilocal */
            fittraj( steps,coef,nlocal,NPREDICT);

            projtraj( &ap->x,&ap->vx,steps[0],coef,NPREDICT,dlong,dshort);
            /* y */
            for( ilocal = 0; ilocal < nlocal ; ilocal++)
            {
                steps[ilocal] = (*positions)[3*(i*nlocal + ilocal)+1];
            }/* ilocal */
            fittraj( steps,coef,nlocal,NPREDICT);
            projtraj( &ap->y,&ap->vy,steps[0],coef,NPREDICT,dlong,dshort);
            /* z */
            for( ilocal = 0; ilocal < nlocal ; ilocal++)
            {
                steps[ilocal] = (*positions)[3*(i*nlocal + ilocal)+2];
            }/* ilocal */
            fittraj( steps,coef,nlocal,NPREDICT);
            projtraj( &ap->z,&ap->vz,steps[0],coef,NPREDICT,dlong,dshort);
        } /* i */

    }/* istep */
    free (positions);
    /*
    free (velos);
    free (theforces);
    */
    return 0;
}/* end of routine */

int projtraj( x,vx,x0,coef,ncoef,dlong,dshort)
float *x,*vx,x0,coef[],dlong,dshort;
int ncoef;
{
    float dt ;
    float dts;
    float xc,xs;
    dt = dlong/dshort;
    dts = dlong;
    xs = sin(coef[5]*dt);
    xc = cos(coef[1]*dt +coef[4]*xs );
    xs = sin(coef[1]*dt +coef[4]*xs );
    dts = cos( coef[5]*dt);
    *vx =  -coef[2]*( coef[1] + coef[4]*coef[5]*dts )*xs +
           coef[3]*( coef[1] + coef[4]*coef[5]*dts )* xc ;
    *vx = *vx / (dshort*1.414213562373095) ;
    /*
    *vx = *vx / dshort ;
    *vx = *vx / dshort *.5;
    *x = coef[0] + (coef[2]*xc + coef[3]*xs)*.5;
    */
    *x = coef[0] + (coef[2]*xc + coef[3]*xs)/1.414213562373095;
}

/*
*  fit with
*  a  + a( cos b t ) + a (sin c t);
*/
int fittraj( steps,coef,nstep,ncoef)
int nstep,ncoef;
float steps[],coef[];
{
    float deltaco[NPREDICT];
    float derivs[2*NPREDICT];
    float normal[NPREDICT][NPREDICT] ;
    float jacobean[2*NPREDICT][NPREDICT] ;
    float xold;
    float x,xc,xs,xxs;
    float beta,betad,thestep;
    int i,j,k,iter;
    /* initialize */

    coef[0] = steps[0];
    coef[1] = 0.;
    coef[1] = 0.3;
    coef[2] = .01;
    coef[2] = .0;
    coef[3] = 0.1;
    coef[4] = .0;
    coef[5] = .01;
    xold = 10e10;
    for( iter = 0; iter< 64; iter++)
    {

        /* this is calc'd every time because
        * in the production version it will change
        * every time 
        */
        x = 0.;
        for( i=0 ;i< nstep; i++)
        {
            xxs =  sin(coef[5]*i);
            xc = cos( coef[1]*i + coef[4]*xxs );
            xs = sin( coef[1]*i + coef[4]*xxs );
            jacobean[i][0] = 1;
            jacobean[i][1] = (-coef[2]*xs+ coef[3]*xc)*i;
            jacobean[i][2] = xc;
            jacobean[i][3] = xs;
            jacobean[i][4] = (-coef[2]*xs+ coef[3]*xc)*xxs;
            jacobean[i][5] = i*cos(coef[5]*i)*coef[4]*
                             (-coef[2]*xs + coef[3]*xc);
            /* the error vector */
            derivs[i] =  coef[0] + coef[2]*xc + coef[3]*xs ;
            derivs[i] =  steps[i] - derivs[i];
            x += derivs[i]*derivs[i];
        }
        /* now for the normal matrix and the delta vector */

        for( i=0; i< ncoef; i++)
        {
            deltaco[i] = 0.;
            for( j=0; j< ncoef; j++)
            {
                normal[i][j] = 0.;
            }
        }

        for(i=0; i< ncoef; i++)
        {
            for( j=0 ; j < ncoef ; j++)
            {
                for( k=0 ; k < nstep ; k++)
                {
                    normal[i][j] += jacobean[k][i]*jacobean[k][j];
                }
            }
            for( k= 0; k< nstep; k++)
            {
                deltaco[i] -= jacobean[k][i]*derivs[k];
            }
            normal[i][i] += .1;
        }


        mom_solve( &normal[0][0],&deltaco[0], NPREDICT,NPREDICT );

        coef[0] -= deltaco[0];
        coef[1] -= deltaco[1];
        coef[2] -= deltaco[2];
        coef[3] -= deltaco[3];
        coef[4] -= deltaco[4];
        coef[5] -= deltaco[5];
        if( x < 1.e-7)break;
        if( fabs(x-xold) < 1.e-9)break;
        xold = x;
    }/*iter*/
    /*
    printf("final error %f\n", x);
    */
}
