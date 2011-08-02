/*
* tncnewt.c
*
* truncated newtons method
* optimiser routine for AMMP
*
* theory is to use a finite difference appx to the
* hessian step product and find the optimal step vector
* that way.
*
* basically this is fixed point iteration for the
* difference in gradients being colinear and 
* descent direction along the step vector.
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
#include <stdio.h>
#include <math.h>
#ifdef ANSI
#include <stdlib.h>
#endif
#include "ammp.h"

/* function tncnewt()
* truncated newtons method optimiser
*
* passed a set of functions vfs for the potential
* passed a set of functions ffs for the force
*  how many functions  nfs
*  how many steps to try
*  when to stop
*/
int tncnewt( vfs,ffs,nfs, nstep,dtoler,toler )
int nfs,(*vfs[])(),(*ffs[])();
int nstep;
float toler,dtoler;
{
    int i,j,ifs,istep;
    float vb,vt,vto,vtk,alpha,atop;
    float dfx,dfy,dfz;
    float linmin();
    float a_max_f(),a_max_d(),lam, dstep;
    ATOM *ap,*a_next();
    /*  do at most nstep steps */
    printf(" entering tncnewt \n");
    if( dtoler < 1.e-3) dtoler = .001;
    for( istep=0; istep< nstep ; istep++)
    {
        vb = 0.;
        for( ifs = 0; ifs < nfs; ifs++)
            (*vfs[ifs])( &vb,  0.);
        a_f_zero();
        for( ifs = 0; ifs < nfs; ifs++)
            (*ffs[ifs])( 0.);
        lam = a_max_f();
        printf(" vb, maxf %f %f \n",vb,(lam));
        if( lam <= toler) return 1;
        dstep = dtoler/lam;
        a_ftodx(1.,0.);
        a_ftovx(1.,0.);
        a_ftogx(1.,0.);
        /* initial iteration to find the step vector */
        /* using a small number of trial searches
        * rather than the largish number in tnpack
        *
        * but not checking for convergence 
        */
        for( i=0; i< 5; i++)
        {
            a_f_zero();
            for( ifs = 0; ifs < nfs; ifs++)
                (*ffs[ifs])( dstep);

            ap = a_next(-1);
            /* form the Hd product and alpha */
            alpha = 0.; atop = 0.;
            while( 1)
            {
                if( ap == NULL) break;
                /* this is the product
                * well use it on the fly 
                	ap->fx = (ap->fx-ap->dx)/dtoler;
                	ap->fy = (ap->fy-ap->dy)/dtoler;
                	ap->fz = (ap->fz-ap->dz)/dtoler;
                	remember the force is -dV/dx
                	this is steepest descents :o<
                */
                dfx = (ap->fx-ap->gx) -ap->gx*dstep  ;
                dfy = (ap->fy-ap->gy) -ap->gy*dstep;
                dfz = (ap->fz-ap->gz)  -ap->gz*dstep;
                /*	alpha += ap->dx*ap->dx;
                	alpha += ap->dy*ap->dy;
                	alpha += ap->dz*ap->dz;
                */
                alpha += dfx*dfx + dfy*dfy + dfz*dfz;
                /*	atop += (dfx )*ap->dx;
                	atop += (dfy )*ap->dy;
                	atop += (dfz )*ap->dz;
                *//*
                	atop += dfx*(dfx + ap->gx*dstep);
                	atop += dfy*(dfy + ap->gy*dstep);
                	atop += dfz*(dfz + ap->gz*dstep);
                	*/
                dfx += ap->gx*dstep;
                dfy += ap->gy*dstep;
                dfz += ap->gz*dstep;
                atop += dfx*dfx + dfy*dfy + dfz*dfz;
                if( ap == ap->next) break;
                ap = ap->next;
            }
            lam = sqrt(a_max_f());if( lam < 1.) lam = 1.;

            ap = a_next(-1);

            if( fabs(alpha) > 1.e-5){
                while( 1)
                {
                    if( ap == NULL) break;
                    /* need optimum step size  0.1 is a guess */
                    ap->vx += atop/alpha * ((ap->fx-ap->gx) - ap->gx*dstep ); ap->dx = ap->vx;
                    ap->vy += atop/alpha * ((ap->fy-ap->gy) - ap->gy*dstep ); ap->dy = ap->vy;
                    ap->vz += atop/alpha * ((ap->fz-ap->gz) - ap->gz*dstep ); ap->dz = ap->vz;

                    if( ap == ap->next) break;
                    ap = ap->next;
                }
            }
            lam = sqrt( a_max_d()); if( lam < 1.) lam = 1.;
            lam = 1./lam;
            ap = a_next(-1);
            while(1)
            {
                if( ap == NULL ) break;
                ap->dx *= lam;
                ap->dy *= lam;
                ap->dz *= lam;
                if( ap == ap->next) break;
                ap = ap->next;
            }
            dstep = dtoler;
        }
        ap = a_next(-1);
        while( 1)
        {
            if( ap == NULL) break;
            ap->dx  =  ap->vx;
            ap->dy  =  ap->vy;
            ap->dz  =  ap->vz;
            if( ap == ap->next) break;
            ap = ap->next;
        }

        lam = sqrt( a_max_d()); if( lam < 1.) lam = 1.;
        lam = 1./lam;
        ap = a_next(-1);
        while(1)
        {
            if( ap == NULL ) break;
            ap->dx *= lam;
            ap->dy *= lam;
            ap->dz *= lam;
            if( ap == ap->next) break;
            ap = ap->next;
        }


        lam = linmin( vfs,nfs, 1.);
        if( lam < 1.e-6)
        {
            a_f_zero();
            for( ifs = 0; ifs < nfs; ifs++)
                (*ffs[ifs])( 0.);
            lam = sqrt(a_max_f()); if( lam < 1.) lam = 1.;
            a_ftodx(1./lam, 0.);
            lam = linmin(vfs,nfs,.01);
            if( lam < 1.e-8) return 1;
            vt = 0.;
            for( ifs = 0; ifs < nfs; ifs++)
                (*vfs[ifs])( &vt,lam);
            printf(" steep sub-iter  %f %f\n",vt,vb);
            if( vt > vb) return 0;
        }
        a_inc_d( lam );
    }
    return 0;
}
