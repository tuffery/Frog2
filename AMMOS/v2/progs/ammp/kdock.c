/* kdock.c
*  modified version of kohonen.c for docking
*
* the idea is to use the algorithm in a spatially limited 
*  region on an active zone.
*  use mxdq as a limit on step size and skip if the closest
*  atom is inactive.
*
*/
/*
*  kohonen.c
*
*  implement a kohonen neural net to find a 
* 3-d space filling curve corresponding to the structure;
*
*  net links correspond to chemical bonds, angle distances, noel distances
*  eventually the distance update will also include vdw interactions
*  and chirality.
*
*/
/*
*  copyright 1993-1997 Robert W. Harrison
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

#define noel_target
#define MAX_LOCAL 200

void kdock( vfs,ffs,nfs,niter, radius, initialize ,rx,ry,rz)
int (*vfs[])(),(*ffs[])(),nfs;
int niter;
float radius;
float rx,ry,rz;
int initialize;
{
    ATOM *ap,*bp[MAX_LOCAL],*cp,*a_next();
    ATOM *actp[600]; /* pointers which are frozen */
    int inactp;
    int numatm , a_number();
    int iter,local_iter, i,j;
    float x,y,z,r,rc,randf();
    float cx,cy,cz, dx,dy,dz;
    int v_bond(),v_angle(),u_v_nonbon(),v_noel();
    int v_mmbond(),v_mmangle(),v_c_angle(),v_nonbon();
    int v_hard(),v_abc();
    int v_hybrid();
    int v_step();
    int steep();
    int use_noel,use_hybrid,use_step;
#ifdef WINDOWS
    void force_screen_update(); /* in animate.c */
#endif
    void get_bond(); /* getbond(ap,atom*bp[],max,&howmany) */
    void get_angle(); /* same call */
    void get_noel(); /* same call */
    void get_bond_and_length(); /* get...(ap,atom* bp[],float r[],max &howmany) */
    void get_noel_and_length(); /* same call */
    void get_noel_and_bounds(); /* get...(ap atom* bp,float r[],rup[],max,&howmany) */
    void get_step_and_bounds(); /* similar to above but get r,rmid,rup */
    float upper[MAX_LOCAL],middle[MAX_LOCAL];
    float target[MAX_LOCAL];
    void rand3();
    float k;
    float mxdq, get_f_variable();

    numatm = a_number();
    if( numatm < 2) return;
    if(radius < 1.) radius = sqrt((float)numatm);
    mxdq = get_f_variable("mxdq");
    if( mxdq < 0.1) mxdq = 0.1;

    cx = rx;
    cy = ry;
    cz = rz;

    if( initialize > 0 )
    {
        for( i=0; i< numatm; i++)
        {ap = a_next(i);
REDO_INITIAL_VALUE: ;
            x = 2*randf()-1.;
            y = 2*randf()-1.;
            z = 2*randf()-1.;
            if( (x*x + y*y + z*z)> 1.) goto REDO_INITIAL_VALUE;
            if( ap->active){
                ap->x = x*radius+cx;
                ap->y = y*radius+cy;
                ap->z = z*radius+cz;}
        }
    }
    use_noel = (1==0);
    use_hybrid = use_noel;
    use_step = use_noel;
    for( i=0; i< nfs; i++)
    {
        if( vfs[i] == v_noel){ use_noel = 1==1;}
        if( vfs[i] == v_hybrid) { use_hybrid = 1==1;}
        if( vfs[i] == v_step) { use_step = 1==1;}
    }


    for( iter = 0; iter< niter; iter++)
    {
#ifdef WINDOWS
        if( iter%10) force_screen_update();
#endif



#ifdef WINDOWs
        BE_NICE();
#endif
        k = 0.2;
RE_RANDOM:;
        x = 2.*randf()-1.;
        y = 2.*randf()-1.;
        z = 2.*randf()-1.;
        rc = x*x + y*y + z*z;
        if(rc > 1.) goto RE_RANDOM;
        x = x*radius +cx; y = y*radius +cy ; z = z*radius +cz;


        /* find the closest atom */
        inactp = 0;
        cp = a_next(-1);
        rc = 10.e10;
        for(i=0;i < numatm; i++)
        {ap = a_next(i);
            r  = (ap->x-x)*(ap->x-x);
            r += (ap->y-y)*(ap->y-y);
            r += (ap->z-z)*(ap->z-z);
            if( r < rc){ cp = ap; rc = r;}
        }

        if( !cp->active && rc < 9.) goto RE_RANDOM;
        /* now find the closest active atom */
        cp = a_next(-1);
        rc = 10.e10;
        for( i=0; i< numatm; i++)
        {ap = a_next(i);
            if( ap->active){
                r  = (ap->x-x)*(ap->x-x);
                r += (ap->y-y)*(ap->y-y);
                r += (ap->z-z)*(ap->z-z);
                if( r < rc){ cp = ap; rc = r;}
            }
        }
        /* truncate the step to no more than mxdq */
        if( rc > mxdq*mxdq)
        {
            rc = mxdq/sqrt(rc);
            x = (x-cp->x)*rc+cp->x;
            y = (y-cp->y)*rc+cp->y;
            z = (z-cp->z)*rc+cp->z;
        }
        /* update it and its neighbors */

        for( i=0; i< numatm; i++)
        {
            ap = a_next(i);
            if( ap->active && ap != cp){
                r  = (ap->x-cp->x)*(ap->x-cp->x);
                r += (ap->y-cp->y)*(ap->y-cp->y);
                r += (ap->z-cp->z)*(ap->z-cp->z);
                if( r < 25.)/* 25 (r==5) works well with good noe */
                {
                    ap->x += (x-cp->x)*k;
                    ap->y += (y-cp->y)*k;
                    ap->z += (z-cp->z)*k;
                    actp[inactp++] = ap;
                }
            }
        }


        if( cp->active){
            dx = x-cp->x;
            dy = y-cp->y;
            dz = z-cp->z;
            cp->x += (dx)*k;
            cp->y += (dy)*k;
            cp->z += (dz)*k;

            actp[inactp++] = cp;
        } else { /* put the center at the atom */
            x = cp->x; y = cp->y; z = cp->z;
            dx = 0.; dy = 0.; dz = 0.;
        }

        get_bond_and_length(cp,bp,target,20,&j);
        for( i=0; i< j ; i++)
        {ap = bp[i];
            if( ap->active){
                ap->x += (x-ap->x)*k;
                ap->y += (y-ap->y)*k;
                ap->z += (z-ap->z)*k;
                actp[inactp++] = ap;
                r = (x - ap->x)*(x -ap->x);
                r += (y - ap->y)*(y -ap->y);
                r += (z - ap->z)*(z -ap->z);
                if( r > 0.){
                    r = target[i]/sqrt(r);
                    ap->x = x + (ap->x-x)*r;
                    ap->y = y + (ap->y-y)*r;
                    ap->z = z + (ap->z-z)*r;
                }
            }
        }

        if( use_step){
            get_step_and_bounds(cp,bp,target,middle,upper,MAX_LOCAL,&j);
            for( i=0; i<j; i++)
            {	ap = bp[i];
                if( ap->active){
                    /* define either step_metric for metrical usage
                    or step_neighborhood for neighborhood useage
                    */
#define STEP_NEIGHBORHOOD
#ifdef STEP_METRIC
                    r = (x - ap->x)*(x -ap->x);
                    r += (y - ap->y)*(y -ap->y);
                    r += (z - ap->z)*(z -ap->z);
                    if( r <= upper[i]*upper[i] && r > 0.){
                        if( r > middle[i]*middle[i] ) r = middle[i]/sqrt(r);
                        else if( r < target[i]*target[i] ) r = target[i]/sqrt(r);
                        else r = 1.;
                        ap->x = x + (ap->x -x)*r;
                        ap->y = y + (ap->y -y)*r;
                        ap->z = z + (ap->z -z)*r;
#endif
#ifdef STEP_NEIGHBORHOOD
                        r = (cp->x -ap->x)*(cp->x -ap->x);
                        r += (cp->y -ap->y)*(cp->y - ap->y);
                        r += (cp->z - ap->z)*(cp->z - ap->z);
                        if( r < upper[i]*upper[i] ){

                            ap->x += dx*k;
                            ap->y += dy*k;
                            ap->z += dz*k;
#endif
                            actp[inactp++] = ap;
                        }/* close enough to work on */
                    }/* active */
                }/* for i */
            }/* use_step*/
            if( use_noel){
#ifdef noel_target
                get_noel_and_length(cp,bp,target,MAX_LOCAL,&j);
#else
get_noel_and_bounds(cp,bp,target,upper,MAX_LOCAL,&j);
#endif
                for( i=0; i< j ; i++)
                {ap = bp[i];
                    if( ap->active){

                        ap->x += dx*k;
                        ap->y += dy*k;
                        ap->z += dz*k;
                        actp[inactp++] = ap;

                        r = (x-ap->x)*(x-ap->x);
                        r += (y-ap->y)*(y-ap->y);
                        r += (z-ap->z)*(z-ap->z);

                        if( r > 0.){
#ifdef noel_target
                            r = target[i]/sqrt(r);

#else
                            r = sqrt(r);
                            if( r < target[i])
                            {
                                r = target[i]/r;
                                ap->x = x + (ap->x -x)*r;
                                ap->y = y + (ap->y -y)*r;
                                ap->z = z + (ap->z -z)*r;
                            }else if( r > upper[i])
                            {
                                r = upper[i]/r;
                                ap->x = x + (ap->x -x)*r;
                                ap->y = y + (ap->y -y)*r;
                                ap->z = z + (ap->z -z)*r;

                            }
#endif
                        }

                    }
                }
            }

            if( inactp > 0 ){
                if( use_hybrid){
                    for( i=0; i< inactp; i++)
                        gsdg_hybrid(actp[i]);}
                for( i=0; i< inactp; i++)
                    actp[i]->active = (1==0);
                kohonen_minimizer(vfs,ffs,nfs,1);
                for( i=0; i< inactp; i++)
                    actp[i]->active = (1==1);
            }


        }/* end of iter loop */
    }/* end of routine kdock*/
