/*
* trust.c
*
* trust based
* optimiser routine for AMMP
*
* theory is to move in a steepest descent direction until
* the disagreement between observed and calculated energies 
* is too high  ( use a .1 v_init tolerance zB) 
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
/* function trust()
* trust region based optimiser
*
* passed a set of functions vfs for the potential
* passed a set of functions ffs for the force
*  how many functions  nfs
*  how many steps to try
*  when to stop
*/
int trust( vfs,ffs,nfs, nstep,dtoler,toler )
int nfs,(*vfs[])(),(*ffs[])();
int nstep;
float toler,dtoler;
{
    int i,j,ifs;
    float vb,vt,vto,vtk;
    float a_max_f(),lam,vtoler,v_predict();
    /*  do at most nstep steps */
    printf(" entering trust \n");
    a_v_zero();
    if( dtoler < 1.e-3) dtoler = .1;
    for( i=0; i< nstep ; i++)
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
        if( lam <= 1. )
        { a_ftodx(1.,0.); } else { a_ftodx( .1/sqrt(lam),0. ); }
        vto = vb;
        vtk = vb;
        lam = 0.;
        /* this is a fibonacci line search
        * rather odd, the steps keep getting
        * bigger until
        * either the end
        * or we go up
        */
        vtoler = abs(vb*dtoler);
        if( vtoler < .1) vtoler = .1;
        for( j=0; j< 200; j++)
        {
            vt = 0.;
            lam = lam + .001*j;
            for( ifs = 0; ifs < nfs; ifs++)
                (*vfs[ifs])(&vt, lam);
            vto = v_predict( lam ,vb);
            if( (vt -vto) > vtoler )
            {  lam = lam -.001*j; break; }
            if( vt > vtk && j != 0 )
            { lam = lam -.001*j; break; }
            vtk = vt;

        }
        a_inc_d( lam );
    }
    vb = 0.;
    for( ifs = 0; ifs < nfs; ifs++)
        (*vfs[ifs])( &vb,  0.);
    printf(" vb %f \n",vb);
    return 0;
}
/* function v_predict( lam,vinit )
*
* predicts the potential based on the initial force
* vector
*
* dv =  -F dx 
*/
float v_predict( lambda,vinit )
float lambda,vinit;
{
    ATOM *ap,*a_next();
    float vp;
    ap = a_next(-1);
    vp = vinit;
    if( ap == NULL) return vp;
    while(1)
    {
        if( ap->next == NULL) return vp;
        vp -= ap->fx*ap->dx*lambda;
        vp -= ap->fy*ap->dy*lambda;
        vp -= ap->fz*ap->dz*lambda;
        if( ap == ap->next) return vp;
        ap = ap->next;
    }
}
