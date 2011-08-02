/*
* optimist.c
*
* optimiser routines for AMMP
*
*  steepest descents
*  conjugate gradients
*  line search    routines
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

/* function steep()
* steepest descents optimiser
*
* passed a set of functions vfs for the potential
* passed a set of functions ffs for the force
*  how many functions  nfs
*  how many steps to try
*  when to stop
*/
int steep( vfs,ffs,nfs, nstep,toler )
int nfs,(*vfs[])(),(*ffs[])();
int nstep;
float toler;
{
    int i,j,ifs;
    float vb,vt,vto;
    float a_max_f(),lam;
    /*  do at most nstep steps */
    for( i=0; i< nstep ; i++)
    {
        vb = 0.;
        for( ifs = 0; ifs < nfs; ifs++)
            (*vfs[ifs])( &vb,  0.);
        a_f_zero();
        for( ifs = 0; ifs < nfs; ifs++)
            (*ffs[ifs])( 0.);
        lam = a_max_f();
        if( lam <= toler) return 1;
        if( lam <= 1. )
        { a_ftodx(1.,0.); } else { a_ftodx( 1./sqrt(lam),0. ); }
        vto = vb;
        lam = 0.;
        /* this is a fibonacci line search
        * rather odd, the steps keep getting
        * bigger until
        * either the end
        * or we go up
        */
        for( j=0; j< 200; j++)
        {
            vt = 0.;
            lam = lam + .01*j;
            for( ifs = 0; ifs < nfs; ifs++)
                (*vfs[ifs])(&vt, lam);
            if( vt > vto && j == 0 )
            { a_ftodx(0., .1);   }
            else if( vt > vto )
            { lam = lam -.01*j; break; }
            vto = vt;
        }
        a_inc_d( lam );
    }
    return 0;
}

/* function cngdel()
* conjugate gradients optimiser
*
* passed a set of functions vfs for the potential
* passed a set of functions ffs for the force
*  how many functions  nfs
*  how many steps to try
*  when to stop
*/
int cngdel( vfs,ffs,nfs, nstep,nreset,toler,echo )
int nfs,(*vfs[])(),(*ffs[])();
int nstep,nreset;
float toler;
{
    int i,j,ifs;
    float vb,vt,vto;
    float a_max_f(),lam;
    float linmin(),a_max_d();
    float a_pr_beta();
    float beta,betad,a_l2_f(),a_l2_g();
    int insane();/* check status */
    int crazy_force();
    /*  do at most nstep steps */
    if( nreset < 1) nreset = nstep;
    a_g_zero();
    a_d_zero();
    for( i=0; i< nstep ; i++)
    {
        if( insane() ) return 0;
        vb = 0.;
        for( ifs = 0; ifs < nfs; ifs++)
            (*vfs[ifs])( &vb,  0.);
        a_f_zero();
        for( ifs = 0; ifs < nfs; ifs++)
            (*ffs[ifs])( 0.);
        if( crazy_force() ) return 0;
        lam = a_max_f();
        /* make up the beta *
        /* use just the simple one */
        /*	beta = a_l2_f();
        	betad = a_l2_g();
        	if( betad < 1.e-4) {betad = 1.;beta = 0.;}
        	beta = -beta/betad;
        */
        beta = a_pr_beta();

        if( (i%nreset) == 0) beta = 0.;
        if( echo) printf(" vb, maxf %f %f %f \n",vb,lam,beta);
        if( lam <= toler) return 1;
        /*	printf(" max f max d  %f %f \n",a_max_f(),a_max_d() ); */
        /* and make up the descent direction */
        lam = a_max_f();
        a_ftodx(1.,beta);
        a_ftogx( 1.,0.);
        lam = linmin( vfs,nfs, sqrt(a_max_d()) );
        if( lam < 1.e-6)
        {
            a_f_zero();
            for( ifs = 0; ifs < nfs; ifs++)
                (*ffs[ifs])(0. );
            lam = sqrt(a_max_f()); if( lam < 1.) lam = 1.;
            a_ftodx( 1./(lam), 0.);
            lam = linmin(vfs,nfs,1. );
            if( lam < 1.e-6) return 0;
        }
        a_inc_d( lam );
    }
    return 0;
}



/*line minimization routine linmin( x, search, n, f )
*
* returns step size for approximate line minimizer
*  recursive line minimizer
*/
float linmin(ffs,nfs ,damp )
int (*ffs[])();
int nfs;
float damp;

{
    int i,iter,jter,imin;
    float alpha[401],fval[401],dstep,step,stpmin,fvt,fmin,fold;
    int nostep,get_i_variable();
    float mxdq, get_f_variable();

    nostep = get_i_variable("nostep");
    if( nostep < 1) nostep = 8;
    mxdq = get_f_variable("mxdq");
    set_f_variable("mxdq" , 100.);

    dstep = 1.;
    step = 0;
    stpmin = 0.;
    imin = 0;
    alpha[0] = 0;
    fval[0] = 0;
    /*	 if( damp < 1.) {dstep = .25;}else{ dstep = .01/sqrt(damp);} */
if( damp < 1.) {dstep = .25;}else{ dstep = 1./damp;}

    for( i=0; i< nfs ; i++)
        (*ffs[i])(&fval[0],0.);
    fmin = fval[0]; fold = fmin;
    /* recheck is to find a descent stp length */
    imin++;
recheck:
    alpha[imin] = dstep;
    fval[imin] = 0.;
    for( i=0; i< nfs ; i++)
        (*ffs[i])(&fval[imin],dstep);
    if( fval[imin] > fval[0])
    {
        dstep = dstep*.25;
        if( dstep > 1.e-8) goto recheck;
        goto cleanup;
    }
    /* if here we have dstep small enough to use */
    /*	for( iter=0; iter< 8; iter++)
    */
    for( iter=0; iter< nostep; iter++)
    {
        for(jter=1; jter<100; jter++)
        {
            step += dstep;
            /* get the function value */
            for( i=0; i< imin; i++)
            {
                if( alpha[i] == step)
                {
                    fvt = fval[i]; goto FOUND ;
                }
            }
            fvt = 0.;
            for( i=0; i< nfs ; i++)
                (*ffs[i])(&fvt,step);
            alpha[imin] = step;  fval[imin++] = fvt;
            if( imin > 400) imin = 400;
FOUND:
            /*		printf("\n %f %f %f %f", step,fvt,fmin,fold);  */
            if( fvt< fmin)
            {
                fmin = fvt; stpmin = step;
            }
            if( fvt > fold)
            {
                dstep = -dstep/2; break;
            }
            fold = fvt;
        }

    }
cleanup:
    /* do not allow 0 steps */
    /*	if( stpmin == 0.) stpmin = dstep/2; */
    set_f_variable("mxdq" , mxdq);
    return(stpmin);
}
