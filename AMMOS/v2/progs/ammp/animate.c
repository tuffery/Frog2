/*  animate.c
*
* routines for performing molecular dynamics
*
* v_maxwell( float Temperature, driftx,drifty,driftz)  initialize velocities
*                           with maxwell distribution, assuming
*                           simple kinetic theory relating T and v
*                           driftx,drifty,driftz allow the use of a constant 
*                           drift velocity.
*
* int v_rescale( float Temperature)
*  rescale velocities so that ke is 3RT/2M
*
* int verlet(forces,nforces, nstep,dtime)
*           perform verlet (forward Euler) dynamics
*                           
* int pac( forces,nforces, nstep,dtime)
*            predict and correct dynamics
*
* int pacpac( forces,nforces,nstep,dtime)
*             iterated pac dynamics
*
*
* int tpac( forces,nforces, nstep,dtime, T)
*  nose constrained dynamics
* int ppac( forces,nforces, nstep,dtime, P)
*   pressure only constrained
* int ptpac( forces,nforces, nstep,dtime,P, T)
*   pressure and temperature constrained
* int hpac( forces,nforces, nstep,dtime, H)
*  total energy  constrained dynamics
*/
/*  7/99 add a sanity check and return if
*  we have a nan */
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
* but otherwise is self-contained. Note the hooks for Non-bonded potentials
*/
/* routine insane checks for nan and returns on error conditions */
int insane()
{
    ATOM *ap, *a_next();
    int na,a_number();
    int sane;
    int i;

    na = a_number();
    if( na == 0) return (1==0);
    sane = (1==0);
    for( i=0; i< na; i++)
    {
        ap = a_next(i);
        if( finite(ap->x) == 0){ sane= (1==1);}
        if( finite(ap->y) == 0){ sane= (1==1);}
        if( finite(ap->z) == 0){ sane= (1==1);}

    }
    return sane;
}
int crazy_force()
{
    ATOM *ap, *a_next();
    int na,a_number();
    int sane;
    int i;
    na = a_number();
    if( na == 0) return (1==0);
    sane = (1==0);
    for( i=0; i< na; i++)
    {
        ap = a_next(i);
        if( finite(ap->fx) == 0){ sane= (1==1);}
        if( finite(ap->fy) == 0){ sane= (1==1);}
        if( finite(ap->fz) == 0){ sane= (1==1);}

    }
    return sane;
}
int v_maxwell( T,dx,dy,dz)
float T,dx,dy,dz;
{
    float randg();
    void rand3();
    ATOM *ap,*a_next(),*bonded[10];
    int iflag,inbond;
    float vmag;
    float R;
    R = 1.987 ; /* kcal/mol/K */

    iflag = -1;
    while( (ap= a_next(iflag++)) != NULL)
    {
        iflag = 1;
        if( ap->mass > 0.)
        {
            /* convert from kcal to mks */
            /* 4.184 to joules */
            /* 1000 grams to kg */
            vmag = sqrt( 3.*R*T/ap->mass*4.184*1000.)* randg();
            rand3( &ap->vx,&ap->vy,&ap->vz);
	ap->vw = 0.;
	ap->dw = 0.;
            if( ap->active ){
                ap->vx = ap->vx*vmag + dx;
                ap->vy = ap->vy*vmag + dy;
                ap->vz = ap->vz*vmag + dz;
                ap->vw =0.; 
                ap->dw =0.; 
            }else{
                ap->vx = 0.;
                ap->vy = 0.;
                ap->vz = 0.;
                ap->vw =0.; 
                ap->dw =0.; 
            }
        }
    }
    /* now check those who are zero mass */
    /* and give them the velocity of the first bonded atom */
    iflag = -1;
    while( (ap= a_next(iflag)) != NULL)
    {
        iflag = 1;
        if( ap->mass <= 0.)
        {
            get_bond(ap,bonded,10,&inbond);
            if( inbond >= 0)
            {
                ap->vx = bonded[0]->vx;
                ap->vy = bonded[0]->vy;
                ap->vz = bonded[0]->vz;
                ap->vw =0.; 
                ap->dw =0.; 
            }
        }
    }
    return 1;
}
/* v_rescale(T)
*  rescale the velocities for constant KE  == Temperature 
*/
int v_rescale( T )
float T;
{
    ATOM *ap,*a_next();
    int iflag,a_number();
    float vmag,KE,target;
    float R;
    R = 1.987 ; /* kcal/mol/K */

    target = 0.;
    target += .5*(3.*R*T)*4.184*1000*a_number();
    KE = 0.;
    iflag = -1;
    while( (ap= a_next(iflag++)) != NULL)
    {
	ap->vw = 0.;
	ap->dw = 0.;
        iflag = 1;
        if( ap->mass > 0.)
        {
            vmag = ap->vx*ap->vx+ap->vy*ap->vy+ap->vz*ap->vz;
            KE += ap->mass*vmag;
        }}
    KE = KE *.5;
    if( KE == 0.)
    { aaerror(" Cannot rescale a zero velocity field -use v_maxwell");
        return 0;
    }
    vmag = sqrt(target/KE);

    iflag = -1;
    while( (ap= a_next(iflag++)) != NULL)
    {
        iflag = 1;
        ap->vx = ap->vx*vmag;
        ap->vy = ap->vy*vmag;
        ap->vz = ap->vz*vmag;
                ap->vw =0.; 
                ap->dw =0.; 
    }
    return 1;
}
/* routine verlet( nstep,dtime)
*int verlet(forces,nforces, nstep,dtime)
*
* perform nstep leapfrogging dynamics with dtime
*/
int verlet(forces,nforces, nstep,dtime)
int (*forces[])(),nforces;
int nstep;
float dtime;
{
    ATOM *bp,*ap,*a_next(),*bonded[10];
    int inbond,iflag;
    int a_f_zero(),a_inc_v();
    int istep,iforces;
    int i,imax,a_number();
    for( istep = 0.; istep< nstep; istep++)
    {
        if( insane() ) return 0;

        /*  find the force at the midpoint */
        a_f_zero();
        for( iforces=0;iforces<nforces; iforces++)
            (*forces[iforces])( 0.);
        /* update velocities */
        imax = a_number();
        ap = a_next(-1);
        bp = ap;
        for( i=0; i< imax; i++,ap = bp)
        {
            bp = a_next(1);
            if( ap->mass > 0.)
            {
                /* the magic number takes kcal/molA to mks */
                ap->vx += ap->fx/ap->mass*dtime*4.184e6;
                ap->vy += ap->fy/ap->mass*dtime*4.184e6;
                ap->vz += ap->fz/ap->mass*dtime*4.184e6;
            }
        }
        imax = a_number();
        ap = a_next(-1);
        bp = ap;
        for( i=0; i< imax; i++,ap = bp)
        {
            bp = a_next(1);
            if( ap->mass <= 0.)
            {
                get_bond(ap,bonded,10,&inbond);
                if( inbond >= 0)
                {
                    ap->vx = bonded[0]->vx;
                    ap->vy = bonded[0]->vy;
                    ap->vz = bonded[0]->vz;
                }
            }
        }
        /* update positions */
        a_inc_v(dtime);
    }/* end of istep loop */
    return 1;
}
/* routine pac( nstep,dtime)
*int pac(forces,nforces, nstep,dtime)
*
* perform nstep pac dynamics with dtime
*
* predict the path given current velocity
* integrate the force (simpson's rule)
*  predict the final velocity
*  update the position using trapezoidal correction
*  
*  ideally several cycles are good
*/
int pac(forces,nforces, nstep,dtime)
int (*forces[])(),nforces;
int nstep;
float dtime;
{
    ATOM *ap,*bp,*a_next(),*bonded[10];
    int inbond,iflag;
    int a_f_zero(),a_inc_v();
    int istep,iforces;
    int i,imax,a_number();

    for( istep = 0.; istep< nstep; istep++)
    {
        if( insane() ) {printf(" insane returns true\n"); return 0;}
        /*  move the velocity vector into the displacment slot */
        imax = a_number();
        ap = a_next(-1);
        bp = ap;
        for( i=0; i< imax; i++,ap = bp)
        {
            bp = a_next(1);
            ap->dx = ap->vx;
            ap->dy = ap->vy;
            ap->dz = ap->vz;
        }

        /*  find the force at the midpoint */
        a_f_zero();
        for( iforces=0;iforces<nforces; iforces++)
            iflag = (*forces[iforces])( dtime/2.);
        /* update velocities */
        imax = a_number();
        ap = a_next(-1);
        bp = ap;
#pragma _CNX no_recurrence
        for( i=0; i< imax; i++,ap = bp)
        {
            bp = a_next(1);
            if( ap->mass > 0.)
            {
                /* the magic number takes kcal/molA to mks */
                /*		ap->vx += ap->fx/ap->mass*dtime*4.184e6/6.;
                *		ap->vy += ap->fy/ap->mass*dtime*4.184e6/6.;
                *		ap->vz += ap->fz/ap->mass*dtime*4.184e6/6.;
                */
                ap->vx =  ap->dx +  ap->fx/ap->mass*dtime*4.184e6;
                ap->vy = ap->dy  + ap->fy/ap->mass*dtime*4.184e6;
                ap->vz = ap->dz  + ap->fz/ap->mass*dtime*4.184e6;
            }
        }
        imax = a_number();
        ap = a_next(-1);
        bp = ap;
#pragma _CNX no_recurrence
        for( i=0; i< imax; i++,ap = bp)
        {
            bp = a_next(1);
            if( ap->mass <= 0.)
            {
                get_bond(ap,bonded,10,&inbond);
                if( inbond >= 0)
                {
                    ap->vx = bonded[0]->vx;
                    ap->vy = bonded[0]->vy;
                    ap->vz = bonded[0]->vz;
                }
            }
        }
        /* update positions */
        imax = a_number();
        ap = a_next(-1);
        bp = ap;
#pragma _CNX no_recurrence
        for( i=0; i< imax; i++,ap = bp)
        {
            bp = a_next(1);
            iflag = 1;
            if( ap->active){
                ap->x += .5*(ap->vx + ap->dx)*dtime;
                ap->y += .5*(ap->vy + ap->dy)*dtime;
                ap->z += .5*(ap->vz + ap->dz)*dtime;
            }
        }

    }/* end of istep loop */
    return 1;
}
/* routine tpac( nstep,dtime)
*int tpac(forces,nforces, nstep,dtime,T)
*
* perform nstep pac dynamics with dtime
* kinetic energy constraint to (3*natom-1) kT/2
*
* predict the path given current velocity
* integrate the force (simpson's rule)
*  predict the final velocity
*  update the position using trapezoidal correction
*  
*  ideally several cycles are good
*
* adaptive steps (6/19/96)
*  if the rescale is too large (i.e. > 2) do two half steps
*
*/
int tpac(forces,nforces, nstep,dtime_real,T)
int (*forces[])(),nforces;
int nstep;
float dtime_real,T;
{
    ATOM *ap,*bp,*a_next(),*bonded[10];
    float ke,Tke,R;
    float alpha;
    float dtime;
    int inbond,iflag;
    int a_f_zero(),a_inc_v();
    int istep,iforces;
    int i,imax,a_number();
    R = 1.987; /* kcal/mol/K */
    for( istep = 0.; istep< nstep; istep++)
    {

        /*	if( insane() ) return 0; */
        /*  move the velocity vector into the displacment slot */
        ke = 0.;
        imax = a_number();
        ap = a_next(-1);
        bp = ap;
        for( i=0; i< imax; i++,ap = bp)
        {
            bp = a_next(1);
            ke += ap->mass*(
                      ap->vx*ap->vx + ap->vy*ap->vy + ap->vz*ap->vz);
            ap->dx = ap->vx;
            ap->dy = ap->vy;
            ap->dz = ap->vz;
        }
        Tke = 3*(imax)*R*4.184*1000;  /* converted into MKS */
        Tke = ke/Tke;  /* Tke is now the current temperature */
        /* scale the current velocities */
        dtime = dtime_real;
        if( Tke > 1.e-6)
        {
            ke = sqrt(T/Tke); /* ke is the scaled shift value */
            dtime = dtime_real/ke;
            /* 0.00002 is 2fs, this is near the limit so don't use it */
            if( dtime > 0.000015 ){
                tpac(forces,nforces,1,dtime_real*0.5,T);
                tpac(forces,nforces,1,dtime_real*0.5,T);
                goto SKIP;
            }
            ap = a_next(-1);
            bp =  ap;
            for( i=0; i< imax;  i++, ap = bp)
            {
                bp = a_next(1);
                ap->dx *= ke;
                ap->dy *= ke;
                ap->dz *= ke;
            }
        }

        /*  find the force at the midpoint */
        a_f_zero();
        for( iforces=0;iforces<nforces; iforces++)
            iflag = (*forces[iforces])( dtime/2.);
        /* update velocities */
        imax = a_number();
        ap = a_next(-1);
        bp = ap;
#pragma _CNX no_recurrence
        for( i=0; i< imax; i++,ap = bp)
        {
            bp = a_next(1);
            if( ap->mass > 0.)
            {
                ap->vx = ap->dx  + ap->fx/ap->mass*dtime*4.184e6;
                ap->vy = ap->dy  + ap->fy/ap->mass*dtime*4.184e6;
                ap->vz = ap->dz  + ap->fz/ap->mass*dtime*4.184e6;
            }
        }
        imax = a_number();
        ap = a_next(-1);
        bp = ap;
#pragma _CNX no_recurrence
        for( i=0; i< imax; i++,ap = bp)
        {
            bp = a_next(1);
            if( ap->mass <= 0.)
            {
                get_bond(ap,bonded,10,&inbond);
                if( inbond >= 0)
                {
                    ap->vx = bonded[0]->vx;
                    ap->vy = bonded[0]->vy;
                    ap->vz = bonded[0]->vz;
                }
            }
        }
        /* update positions */
        imax = a_number();
        ap = a_next(-1);
        bp = ap;
#pragma _CNX no_recurrence
        for( i=0; i< imax; i++,ap = bp)
        {
            bp = a_next(1);
            iflag = 1;
            ap->x += .5*(ap->vx + ap->dx)*dtime;
            ap->y += .5*(ap->vy + ap->dy)*dtime;
            ap->z += .5*(ap->vz + ap->dz)*dtime;
        }
SKIP: ; /* if we are here from goto we have done two half steps (or more)*/

    }/* end of istep loop */
    return 1;
}
/* routine pacpac( nstep,dtime)
*int pacpac(forces,nforces, nstep,dtime)
*
* perform nstep pac dynamics with dtime
*
* predict the path given current velocity
* integrate the force (simpson's rule)
*  predict the final velocity
*  update the position using trapezoidal correction
*  
*  ideally several cycles are good
*/
int pacpac(forces,nforces, nstep,dtime)
int (*forces[])(),nforces;
int nstep;
float dtime;
{
    ATOM *ap,*a_next(),*bp,*bonded[10];
    int inbond,iflag;
    int a_f_zero(),a_inc_v();
    int istep,iforces,icorrect;
    int i,imax,a_number();
    for( istep = 0.; istep< nstep; istep++)
    {
        if( insane() ) return 0;

        /*  move the velocity vector into the displacment slot */
        iflag = -1;
        while( (ap=a_next(iflag)) != NULL)
        {
            iflag = 1;
            ap->dx = ap->vx;
            ap->dy = ap->vy;
            ap->dz = ap->vz;
        }

        /*  find the force at the midpoint */
        a_f_zero();
        for( iforces=0;iforces<nforces; iforces++)
            iflag = (*forces[iforces])( dtime/2.);
        /* update velocities */
        iflag = -1;
        while( (ap=a_next(iflag)) != NULL)
        {
            iflag = 1;
            if( ap->mass > 0.)
            {
                /* the magic number takes kcal/molA to mks */
                ap->gx = ap->vx;
                ap->gy = ap->vy;
                ap->gz = ap->vz;
                ap->vx += ap->fx/ap->mass*dtime*4.184e6;
                ap->vy += ap->fy/ap->mass*dtime*4.184e6;
                ap->vz += ap->fz/ap->mass*dtime*4.184e6;
            }
        }
        iflag = -1;
        while( (ap=a_next(iflag)) != NULL)
        {
            iflag = 1;
            if( ap->mass <= 0.)
            {
                ap->gx = ap->vx;
                ap->gy = ap->vy;
                ap->gz = ap->vz;
                get_bond(ap,bonded,10,&inbond);
                if( inbond >= 0)
                {
                    ap->vx = bonded[0]->vx;
                    ap->vy = bonded[0]->vy;
                    ap->vz = bonded[0]->vz;
                }
            } /* end of mass check */
        }
        /* make up the new prediction direction */
        imax = a_number();
        ap = a_next(-1);
        bp = ap;
#pragma _CNX no_recurrence
        for( i=0; i< imax; i++,ap = bp)
        {
            bp = a_next(1);
            ap->dx = ap->vx + ap->gx;
            ap->dy = ap->vy + ap->gy;
            ap->dz = ap->vz + ap->gz;
        }
        for( icorrect = 0;icorrect < 2; icorrect ++)
        {
            a_f_zero();
            for( iforces=0;iforces<nforces; iforces++)
                iflag = (*forces[iforces])( dtime/4.);
            /* update velocities */
            imax = a_number();
            ap = a_next(-1);
            bp = ap;
#pragma _CNX no_recurrence
            for( i=0; i< imax; i++,ap = bp)
            {
                bp = a_next(1);
                if( ap->mass > 0.)
                {
                    /* the magic number takes kcal/molA to mks */
                    ap->vx = ap->gx + ap->fx/ap->mass*dtime*4.184e6;
                    ap->vy = ap->gy + ap->fy/ap->mass*dtime*4.184e6;
                    ap->vz = ap->gz + ap->fz/ap->mass*dtime*4.184e6;
                }
            }
            imax = a_number();
            ap = a_next(-1);
            bp = ap;
#pragma _CNX no_recurrence
            for( i=0; i< imax; i++,ap = bp)
            {
                bp = a_next(1);
                if( ap->mass <= 0.)
                {
                    get_bond(ap,bonded,10,&inbond);
                    if( inbond >= 0)
                    {
                        ap->vx = bonded[0]->vx;
                        ap->vy = bonded[0]->vy;
                        ap->vz = bonded[0]->vz;
                    }
                } /* end of mass check */
            }
            /* make up the new prediction direction */
            imax = a_number();
            ap = a_next(-1);
            bp = ap;
#pragma _CNX no_recurrence
            for( i=0; i< imax; i++,ap = bp)
            {
                bp = a_next(1);
                ap->dx = ap->vx + ap->gx;
                ap->dy = ap->vy + ap->gy;
                ap->dz = ap->vz + ap->gz;
            }
        }
        /* update positions */
        imax = a_number();
        ap = a_next(-1);
        bp = ap;
#pragma _CNX no_recurrence
        for( i=0; i< imax; i++,ap = bp)
        {
            bp = a_next(1);
            ap->x += .5*(ap->vx + ap->gx)*dtime;
            ap->y += .5*(ap->vy + ap->gy)*dtime;
            ap->z += .5*(ap->vz + ap->gz)*dtime;
        }

    }/* end of istep loop */
    return 1;
}
/* routine hpac( nstep,dtime)
*int hpac(forces,nforces, nstep,dtime,H)
*
* perform nstep pac dynamics with dtime
* kinetic energy adusted for constant H
*
* predict the path given current velocity
* integrate the force (simpson's rule)
*  predict the final velocity
*  update the position using trapezoidal correction
*  
*  ideally several cycles are good
*/
int hpac(forces,poten,nforces,nstep,dtime_real,H)
int (*forces[])(),(*poten[])(),nforces;
int nstep;
float dtime_real,H;
{
    ATOM *ap,*bp,*a_next(),*bonded[10];
    float ke,Tke;
    float alpha;
    float dtime;
    int inbond,iflag;
    int a_f_zero(),a_inc_v();
    int istep,iforces;
    int i,imax,a_number();


    for( istep = 0.; istep< nstep; istep++)
    {
        if( insane() ) return 0;

        /*  move the velocity vector into the displacment slot */
        ke = 0.;
        imax = a_number();
        ap = a_next(-1);
        bp = ap;
        for( i=0; i< imax; i++,ap = bp)
        {
            bp = a_next(1);
            ke += ap->mass*(
                      ap->vx*ap->vx + ap->vy*ap->vy + ap->vz*ap->vz);
            ap->dx = ap->vx;
            ap->dy = ap->vy;
            ap->dz = ap->vz;
        }
        ke = ke*.5/4.184/1000/1000;  /* ke in kcal/mol */
        /* get the current potential */
        Tke = 0.;
        for(i=0; i< nforces; i++)
            (*poten[i])(&Tke,0.);
        /* scale the current velocities */
        dtime = dtime_real;
        if( Tke < H )
        {
            ke = sqrt((H-Tke)/ke); /* ke is the scaled shift value */
            dtime = dtime_real/ke;
            /* 0.00002 is 2fs, this is near the limit so don't use it */
            if( dtime > 0.000020 ){
                hpac(forces,poten,nforces,1,dtime_real*0.5,H);
                hpac(forces,poten,nforces,1,dtime_real*0.5,H);
                goto SKIP;
            }
            ap = a_next(-1);
            bp =  ap;
            for( i=0; i< imax;  i++, ap = bp)
            {
                bp = a_next(1);
                ap->dx *= ke;
                ap->dy *= ke;
                ap->dz *= ke;
            }
        } else {
            aaerror("Warning in Hpac, Potential energy higher than target\n");
            a_v_zero();
            a_d_zero();
        }

        /*  find the force at the midpoint */
        a_f_zero();
        for( iforces=0;iforces<nforces; iforces++)
            iflag = (*forces[iforces])( dtime/2.);
        /* update velocities */
        imax = a_number();
        ap = a_next(-1);
        bp = ap;
#pragma _CNX no_recurrence
        for( i=0; i< imax; i++,ap = bp)
        {
            bp = a_next(1);
            if( ap->mass > 0.)
            {
                ap->vx = ap->dx  + ap->fx/ap->mass*dtime*4.184e6;
                ap->vy = ap->dy  + ap->fy/ap->mass*dtime*4.184e6;
                ap->vz = ap->dz  + ap->fz/ap->mass*dtime*4.184e6;
            }
        }
        imax = a_number();
        ap = a_next(-1);
        bp = ap;
#pragma _CNX no_recurrence
        for( i=0; i< imax; i++,ap = bp)
        {
            bp = a_next(1);
            if( ap->mass <= 0.)
            {
                get_bond(ap,bonded,10,&inbond);
                if( inbond >= 0)
                {
                    ap->vx = bonded[0]->vx;
                    ap->vy = bonded[0]->vy;
                    ap->vz = bonded[0]->vz;
                }
            }
        }
        /* update positions */
        imax = a_number();
        ap = a_next(-1);
        bp = ap;
#pragma _CNX no_recurrence
        for( i=0; i< imax; i++,ap = bp)
        {
            bp = a_next(1);
            iflag = 1;
            ap->x += .5*(ap->vx + ap->dx)*dtime;
            ap->y += .5*(ap->vy + ap->dy)*dtime;
            ap->z += .5*(ap->vz + ap->dz)*dtime;
        }
SKIP: ; /* if goto here we've had too large a step and used half steps */

    }/* end of istep loop */
    return 1;
}
/* routine ppac( nstep,dtime)
*int ppac(forces,nforces, nstep,dtime,P)
*
* force the pressure to be constant
* use P = integral ( f . r )dV as 	
* the basis for a diffeomorphism
*   P => kP or Integral( kf.r)dV
*         to enforce pressure
*   r => r/k to enforce physical reality
*   may need to damp this.
*  
*
* perform nstep pac dynamics with dtime
*
* predict the path given current velocity
* integrate the force (simpson's rule)
*  predict the final velocity
*  update the position using trapezoidal correction
*  
*  ideally several cycles are good
*/
int ppac(forces,nforces, nstep,dtime_real,P)
int (*forces[])(),nforces;
int nstep;
float dtime_real,P;
{
    ATOM *ap,*bp,*a_next(),*bonded[10];
    float p,Tp,R;
    float dtime,cx,cy,cz;
    float alpha;
    int inbond,iflag;
    int a_f_zero(),a_inc_v();
    int istep,iforces;
    int i,imax,a_number();
    R = 1.987; /* kcal/mol/K */

    imax = a_number();
    if( imax <= 0 )return 0;
    for( istep = 0.; istep< nstep; istep++)
    {
        if( insane() ) return 0;

        cx = 0.; cy = 0.; cz = 0.;
        /*  move the velocity vector into the displacment slot */
        ap = a_next(-1);
        bp = ap;
        for( i=0; i< imax; i++,ap = bp)
        {
            bp = a_next(1);
            ap->dx = ap->vx;
            ap->dy = ap->vy;
            ap->dz = ap->vz;
            cx += ap->x;
            cy += ap->y;
            cz += ap->z;
        }
        cx /= imax;
        cy /= imax;
        cz /= imax;

        /* calculate the pressure */

        p = 0.;
        Tp = 0.;
        ap = a_next(-1);
        bp = ap;
        for( i=0; i< imax; i++,ap = bp)
        {
            bp = a_next(1);
            p += ap->vx*ap->vx*ap->mass;
            p += ap->vy*ap->vy*ap->mass;
            p += ap->vz*ap->vz*ap->mass;
            Tp += (ap->x-cx)*(ap->x-cx);
            Tp += (ap->y-cy)*(ap->y-cy);
            Tp += (ap->z-cz)*(ap->z-cz);
        }
        Tp = sqrt(Tp/imax);
        Tp = 4*PI/3*Tp*Tp*Tp;
        p = p/imax/Tp*.5; /* now mks molar */
        printf("P %f p %f Tp %f\n",P,p,Tp);
        /* moment shift
        	p = sqrt( P/p);
        	dtime = dtime_real/p;
        */
        dtime = dtime_real;
        /* this is about the steepest volume correction which works !!
          1. + .2/1.2 and 1 + .5/1.5 fail
        */
        p = (1.+.1*pow( p/P, 1./3.))/1.1;

        /* temporary kludge to understand problem */
        ap = a_next(-1);
        bp = ap;
        for( i=0; i< imax; i++,ap = bp)
        {
            bp = a_next(1);
            /*
            		ap->vx *= p;	
            		ap->vy *= p;	
            		ap->vz *= p;	
            		ap->dx *= p;	
            		ap->dy *= p;	
            		ap->dz *= p;	
            */

            ap->x *= p;
            ap->y *= p;
            ap->z *= p;
        }
        /*  find the force at the midpoint */
        a_f_zero();
        for( iforces=0;iforces<nforces; iforces++)
            iflag = (*forces[iforces])( dtime/2.);

        /* update velocities */
        imax = a_number();
        ap = a_next(-1);
        bp = ap;
#pragma _CNX no_recurrence
        for( i=0; i< imax; i++,ap = bp)
        {
            bp = a_next(1);
            if( ap->mass > 0.)
            {
                ap->vx = ap->dx  + ap->fx/ap->mass*dtime*4.184e6;
                ap->vy = ap->dy  + ap->fy/ap->mass*dtime*4.184e6;
                ap->vz = ap->dz  + ap->fz/ap->mass*dtime*4.184e6;
            }
        }
        imax = a_number();
        ap = a_next(-1);
        bp = ap;
#pragma _CNX no_recurrence
        for( i=0; i< imax; i++,ap = bp)
        {
            bp = a_next(1);
            if( ap->mass <= 0.)
            {
                get_bond(ap,bonded,10,&inbond);
                if( inbond >= 0)
                {
                    ap->vx = bonded[0]->vx;
                    ap->vy = bonded[0]->vy;
                    ap->vz = bonded[0]->vz;
                }
            }
        }
        /* update positions */
        imax = a_number();
        ap = a_next(-1);
        bp = ap;
#pragma _CNX no_recurrence
        for( i=0; i< imax; i++,ap = bp)
        {
            bp = a_next(1);
            iflag = 1;
            ap->x += .5*(ap->vx + ap->dx)*dtime;
            ap->y += .5*(ap->vy + ap->dy)*dtime;
            ap->z += .5*(ap->vz + ap->dz)*dtime;
        }

    }/* end of istep loop */
    return 1;
}
/* routine ptpac( nstep,dtime)
*int ptpac(forces,nforces, nstep,dtime,P,T)
*
* force the pressure to be constant
*  
*
* perform nstep pac dynamics with dtime
*
* predict the path given current velocity
* integrate the force (simpson's rule)
*  predict the final velocity
*  update the position using trapezoidal correction
*  
*  ideally several cycles are good
*/
int ptpac(forces,nforces, nstep,dtime_real,P,T)
int (*forces[])(),nforces;
int nstep;
float dtime_real,P,T;
{
    ATOM *ap,*bp,*a_next(),*bonded[10];
    float p,Tp,R;
    float Tk;
    float dtime,cx,cy,cz;
    float alpha;
    int inbond,iflag;
    int a_f_zero(),a_inc_v();
    int istep,iforces;
    int i,imax,a_number();
    R = 1.987; /* kcal/mol/K */

    imax = a_number();
    if( imax <= 0 )return 0;
    for( istep = 0.; istep< nstep; istep++)
    {
        if( insane() ) return 0;

        cx = 0.; cy = 0.; cz = 0.;
        /*  move the velocity vector into the displacment slot */
        ap = a_next(-1);
        bp = ap;
        for( i=0; i< imax; i++,ap = bp)
        {
            bp = a_next(1);
            ap->dx = ap->vx;
            ap->dy = ap->vy;
            ap->dz = ap->vz;
            cx += ap->x;
            cy += ap->y;
            cz += ap->z;
        }
        cx /= imax;
        cy /= imax;
        cz /= imax;

        /* calculate the pressure */

        p = 0.;
        Tp = 0.;
        ap = a_next(-1);
        bp = ap;
        for( i=0; i< imax; i++,ap = bp)
        {
            bp = a_next(1);
            p += ap->vx*ap->vx*ap->mass;
            p += ap->vy*ap->vy*ap->mass;
            p += ap->vz*ap->vz*ap->mass;
            Tp += (ap->x-cx)*(ap->x-cx);
            Tp += (ap->y-cy)*(ap->y-cy);
            Tp += (ap->z-cz)*(ap->z-cz);
        }
        Tp = sqrt(Tp/imax);
        Tp = 4*PI/3*Tp*Tp*Tp;
        Tk = 3*imax*R*4.184*1000;
        Tk = p/Tk;  /* Tk is now the temperature */
        if( Tk < 1.e-5) Tk = 1.;
        p = p/imax/Tp*.5; /* now mks molar  ( kilopascal's because of grams)*/
        printf("P %f p %f Tp %f\n",P,p,Tp);
        /* momentum shift */
        Tk = sqrt(T/Tk);
        dtime = dtime_real/Tk;
        /* 0.00002 is 2fs, this is near the limit so don't use it */
        if( dtime > 0.000020 ){
            ptpac(forces,nforces,1,dtime_real*0.5,P,T);
            ptpac(forces,nforces,1,dtime_real*0.5,P,T);
            goto SKIP;
        }
        /* this is about the steepest volume correction which works !!
          1. + .2/1.2 and 1 + .5/1.5 fail
        also checked that the current 'pressure' is the best to use
        for stable running  
        */
        p = (1.+.1*pow( p/P, 1./3.))/1.1;

        /* temporary kludge to understand problem */
        ap = a_next(-1);
        bp = ap;
        for( i=0; i< imax; i++,ap = bp)
        {
            bp = a_next(1);
            ap->vx *= Tk;
            ap->vy *= Tk;
            ap->vz *= Tk;
            ap->dx *= Tk;
            ap->dy *= Tk;
            ap->dz *= Tk;

            ap->x *= p;
            ap->y *= p;
            ap->z *= p;
        }
        /*  find the force at the midpoint */
        a_f_zero();
        for( iforces=0;iforces<nforces; iforces++)
            iflag = (*forces[iforces])( dtime/2.);

        /* update velocities */
        imax = a_number();
        ap = a_next(-1);
        bp = ap;
#pragma _CNX no_recurrence
        for( i=0; i< imax; i++,ap = bp)
        {
            bp = a_next(1);
            if( ap->mass > 0.)
            {
                ap->vx = ap->dx  + ap->fx/ap->mass*dtime*4.184e6;
                ap->vy = ap->dy  + ap->fy/ap->mass*dtime*4.184e6;
                ap->vz = ap->dz  + ap->fz/ap->mass*dtime*4.184e6;
            }
        }
        imax = a_number();
        ap = a_next(-1);
        bp = ap;
#pragma _CNX no_recurrence
        for( i=0; i< imax; i++,ap = bp)
        {
            bp = a_next(1);
            if( ap->mass <= 0.)
            {
                get_bond(ap,bonded,10,&inbond);
                if( inbond >= 0)
                {
                    ap->vx = bonded[0]->vx;
                    ap->vy = bonded[0]->vy;
                    ap->vz = bonded[0]->vz;
                }
            }
        }
        /* update positions */
        imax = a_number();
        ap = a_next(-1);
        bp = ap;
#pragma _CNX no_recurrence
        for( i=0; i< imax; i++,ap = bp)
        {
            bp = a_next(1);
            iflag = 1;
            ap->x += .5*(ap->vx + ap->dx)*dtime;
            ap->y += .5*(ap->vy + ap->dy)*dtime;
            ap->z += .5*(ap->vz + ap->dz)*dtime;
        }

SKIP: ; /* if goto here we've had too large a step and used half steps */
    }/* end of istep loop */
    return 1;
}
