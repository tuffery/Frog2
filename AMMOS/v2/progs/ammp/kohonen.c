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
#define MAX_LOCAL 700

#ifdef GTK_VERSION
void AMMP_whiz_bang(void);
#endif

void kohonen( vfs,ffs,nfs,niter, radius, initialize ,rx,ry,rz)
int (*vfs[])(),(*ffs[])(),nfs;
int niter;
float radius;
float rx,ry,rz;
int initialize;
{
    ATOM *ap,*bp[MAX_LOCAL],*cp,*a_next();
    ATOM *actp[100+MAX_LOCAL]; /* pointers which are frozen */
    int inactp;
    int numatm , a_number();
    int iter,local_iter, i,j;
    float x,y,z,randf();
    critical_precision r,rc;
    float cx,cy,cz;
    float dx,dy,dz;
    int v_bond(),v_angle(),u_v_nonbon(),v_noel(),v_noe_lin();
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
    void get_step_and_bounds();
    float upper[MAX_LOCAL];
    float middle[MAX_LOCAL];
    float target[MAX_LOCAL];
    void rand3();
    float k;
    int prior_type;

    numatm = a_number();
    if( numatm < 2) return;
    if(radius < 1.) radius = sqrt((float)numatm);
    cx = 0.;
    cy = 0.;
    cz = 0.;
    j = 0;
    for( i=0; i< numatm; i++)
    {
        ap = a_next(i);
        if( !ap->active){ j++;
            cx += ap->x; cy += ap->y; cz += ap->z;}
    }

    if( j > 0){
        r = 1./(float)j;
        cx *= r; cy *= r; cz *= r;
    }

    if( initialize > 0 )
    {
        for( i=0; i< numatm; i++)
        {ap = a_next(i);
            x = 2*randf()-1.;
            y = 2*randf()-1.;
            z = 2*randf()-1.;
            if( ap->active){
                ap->x = x*radius+cx;
                ap->y = y*radius+cy;
                ap->z = z*radius+cz;}
        }
    }
    /* figure out what kind of prior to use
    *  rx >= 0 ry == 0 rz == 0 spherical 
    *  rx > 0 ry > 0 rz == 0  cylindrical
    *  rx > 0 ry > 0 rz > 0   elipsoidal
    *
    *  cylindrical and elipsoidal ignore the radius 
    *  spherical uses the radius value
    */
    prior_type = 0;  /* default to sphere */
    if( rx >= 0. && ry < 1.e-7 && rz < 1.e-7) prior_type = 0;
    if( rx > 0. && ry > 0. && rz < 1.e-7) prior_type = 1;
    if( rx > 0. && ry > 0. && rz > 0.) prior_type = 2;
    use_noel = (1==0);
    use_hybrid = use_noel;
    use_step = use_noel;
    for( i=0; i< nfs; i++)
    {
        if( vfs[i] == v_noel){ use_noel = 1==1;}
        if( vfs[i] == v_noe_lin){ use_noel = 1==1;}
        if( vfs[i] == v_hybrid) { use_hybrid = 1==1;}
        if( vfs[i] == v_step) { use_step = 1==1;}
    }
    for( iter = 0; iter< niter; iter++)
    {


#ifdef WINDOWS
        force_screen_update();
#endif
#ifdef GTK_VERSION
        AMMP_whiz_bang();
#endif

        /* k = (float)(niter-iter)/(float)niter; */
        k = 0.01 + (float)(niter-iter)/(float)niter*0.2;
        for( local_iter=0; local_iter< numatm; local_iter++)
        {
#ifdef WINDOWS 
            BE_NICE();
#endif	
            k = 0.2;
            /* pick a random point in the sphere */
            /* trial of rigid chiral enforcement */
            /*	if( use_hybrid)
            	{
            		for( i=0; i< numatm; i++)
            		{
            			ap = a_next(i);
            			if(ap->active)
            				gsdg_hybrid(ap);

            		}
            	}
            */
            if( prior_type == 0){
                r = radius;
RE_RANDOM:;
                x = 2.*randf()-1.;
                y = 2.*randf()-1.;
                z = 2.*randf()-1.;
                rc = x*x + y*y + z*z;
                if(rc > 1.) goto RE_RANDOM;
                x = x*r +cx; y = y*r +cy ; z = z*r +cz;
            } else if( prior_type == 1){
RE_RANDOM_2D: ;
                x = 2.*randf()-1.;
                y = 2.*randf()-1.;
                rc = x*x + y*y;
                if( rc > 1.) goto RE_RANDOM_2D;
                z = 2.*randf()-1.;
                x = x*rx + cx; y = y*rx + cy; z = z*ry;
            } else if( prior_type == 2){
RE_RANDOM_ELIPSE: ;
                x = 2.*randf()-1.;
                y = 2.*randf()-1.;
                z = 2.*randf()-1.;
                rc = x*x + y*y + z*z;
                if(rc > 1.) goto RE_RANDOM_ELIPSE;
                x = x*rx +cx; y = y*ry +cy ; z = z*rz +cz;
            }
            /* find the closest atom */
            inactp = 0;
            /*
            */
            rc = 10.e10;
            for(i=0;i < numatm; i++)
            {ap = a_next(i);
                if( ap->active ){
                    r  = (ap->x-x)*(ap->x-x);
                    r += (ap->y-y)*(ap->y-y);
                    r += (ap->z-z)*(ap->z-z);
                    if( r < rc){ cp = ap; rc = r;}
                }
            }
            /* update it and its neighbors */
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
                        /*			r = target[i]/sqrt(r)-1.;
                        			ap->x -= (x-ap->x)*r;
                        			ap->y -= (y-ap->y)*r;
                        			ap->z -= (z-ap->z)*r;
                        */			
                        r = target[i]/sqrt(r);
                        ap->x = x + (ap->x-x)*r;
                        ap->y = y + (ap->y-y)*r;
                        ap->z = z + (ap->z-z)*r;
                    }
                }
            }

            /*
            		get_angle(cp,bp,20,&j);
            		for( i=0; i< j ; i++)
            		{ap = bp[i];
            		  ap->x += (x-ap->x)*k*.1;
            		  ap->y += (y-ap->y)*k*.1;
            		  ap->z += (z-ap->z)*k*.1;
            		}
            */
            if( use_step){
                get_step_and_bounds(cp,bp,target,middle,upper,MAX_LOCAL,&j);
                for( i=0; i<j; i++)
                {
                    ap = bp[i];
                    if( ap->active){
#define STEP_NEIGHBORHOOD
#ifdef STEP_METRIC
                        r = (x-ap->x)*(x-ap->x);
                        r += (y-ap->y)*(y-ap->y);
                        r += (z-ap->z)*(z-ap->z);
                        if( r <= upper[i]*upper[i] && r > 0.){
                            if( r > middle[i]*middle[i] ) r = middle[i]/sqrt(r);
                            else if( r < target[i]*target[i] ) r = target[i]/sqrt(r);
                            else r = 1.;
                            ap->x = x+ (ap->x -x)*r;
                            ap->y = y+ (ap->y -y)*r;
                            ap->z = z+ (ap->z -z)*r;
#endif
#ifdef STEP_NEIGHBORHOOD
                            r = (cp->x-ap->x)*(cp->x-ap->x);
                            r += (cp->y-ap->y)*(cp->y-ap->y);
                            r += (cp->z-ap->z)*(cp->z-ap->z);
                            if( r <= upper[i]*upper[i] ){
                                /*
                                			ap->x +=  (x - cp->x)*k;
                                			ap->y +=  (y - cp->y)*k;
                                			ap->z +=  (z - cp->z)*k;
                                */
                                ap->x += (dx)*k;
                                ap->y += (dy)*k;
                                ap->z += (dz)*k;
#endif
                                actp[inactp++] = ap;
                            } /* r in bounds */

                        }/* ap is active */
                    }/* for i */
                }/* end of use_step */
                if( use_noel){
#ifdef noel_target
                    get_noel_and_length(cp,bp,target,MAX_LOCAL,&j);
#else
                    get_noel_and_bounds(cp,bp,target,upper,MAX_LOCAL,&j);
#endif
                    for( i=0; i< j ; i++)
                    {ap = bp[i];
                        if( ap->active){
                            ap->x += (dx)*k;
                            ap->y += (dy)*k;
                            ap->z += (dz)*k;
                            actp[inactp++] = ap;
                            r = (x - ap->x)*(x -ap->x);
                            r += (y - ap->y)*(y -ap->y);
                            r += (z - ap->z)*(z -ap->z);
                            if( r > 0.){
#ifdef noel_target
                                /*
                                			r = target[i]/sqrt(r);
                                			ap->x = x + (ap->x-x)*r;
                                			ap->y = y + (ap->y-y)*r;
                                			ap->z = z + (ap->z-z)*r;
                                */
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

            }/* end of local_iter */
        }/* end of iter loop */
    }/* end of routine */
    /* copy of cngdel ,buggerized for one step */
    int kohonen_minimizer( vfs,ffs,nfs, nstep )
    int nfs,(*vfs[])(),(*ffs[])();
    int nstep;
    {
        int i,ifs;
        float a_max_f(),lam,vb;
        float linmin(),a_max_d();
        int a_ftodx();
#ifdef WINDOWS
        void force_screen_update(); /* in animate.c */
#endif
        /*  do at most nstep steps */

        a_g_zero();
        a_d_zero();
        for( i=0; i< nstep ; i++)
        {

            vb = 0.;
            for( ifs = 0; ifs < nfs; ifs++)
            {
                (*vfs[ifs])( &vb,  0.);
#ifdef WINDOWS 
                BE_NICE
#endif;
            }
            a_f_zero();
            for( ifs = 0; ifs < nfs; ifs++)
            {
                (*ffs[ifs])( 0.);
#ifdef WINDOWS 
                BE_NICE
#endif;
            }
            lam = a_max_f();
            a_ftodx(1.,0.);

            lam = linmin( vfs,nfs, sqrt(a_max_d()) );
            a_inc_d( lam );
        }
#ifdef WINDOWS
        force_screen_update();
#endif
#ifdef GTK_VERSION
        AMMP_whiz_bang();
#endif

        return 0;
    }


