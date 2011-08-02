/* gsdg.c
*
*  Gauss-Siedel Distance Geometry
*  
*  iteratively solve distance geometry equations
*  one atom at a time, but update the calculated distance estimates
*  each time.  (we know from the PMDG experiments that this is not
*  too expensive)
*/
/*
*  copyright 1993,1994 Robert W. Harrison
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

void gsdg( vfs,nfs,niter,low_serial,high_serial )
int (*vfs[])(),nfs;
int niter;
int low_serial,high_serial;
{
    ATOM *ap,*bp,*a_next();
    int numatm , a_number();
    int iter, i,j,k;
    float svec[3],rvec[3]; /* storage for search vectors */
    float x,y,z;
    float gsdg_line_search();
    int v_bond(),v_angle(),u_v_nonbon(),v_noel();
    int v_mmbond(),v_mmangle(),v_c_angle(),v_nonbon();
    int v_hard(),v_abc(),v_hybrid();

    numatm = a_number();

    if( low_serial > high_serial)
    { i = low_serial; low_serial = high_serial; high_serial = i;}

    if( high_serial <= 0 )
        { for( i=0; i< numatm; i++)
        { ap = a_next(i);
            if( high_serial < ap->serial ) high_serial = ap->serial;
        }
    }

    for( iter = 0; iter< niter; iter++)
    {
        printf(" iter %d starting ",iter);
        fflush(stdout);
        ap = a_next(-1);/* cannot use a_next in this loop */
        for( i=0; i< numatm; i++)
        {
            if( ap == NULL ) break;
            if( ap->serial >= low_serial && ap->serial <= high_serial )
            {
                if( ap->active){

                    for( j=0;j< numatm; j++)
                    { bp = a_next(j);
                        bp->vx = 16.; /* default to 4A separation */
                        bp->vy = 0.; /* but only as a lower bound */
                    }

                    for( j=0; j< nfs; j++)
                    {
                        if( vfs[j] == v_nonbon || vfs[j] == u_v_nonbon)
                            { for( k=0; k< numatm; k++)
                            {bp = a_next(k);
                                bp->vy = -10; }
                            break;
                        }
                    }
                    for( j=0; j< nfs; j++)
                    {
                        if( vfs[j] == v_bond || vfs[j] == v_mmbond || vfs[j] == v_abc) gsdg_bond(ap);
                        if( vfs[j] == v_angle || vfs[j] == v_mmangle || vfs[j] == v_c_angle) gsdg_angle(ap);
                        if( vfs[j] == v_noel) gsdg_noel(ap);
                        if( vfs[j] == v_hybrid) gsdg_hybrid(ap);
                        /*
                        	if( vfs[j] == v_hard) gsdg_nonbon(ap)
                        	if( vfs[j] == v_nonbon || vfs[j] == u_v_nonbon) gsdg_nonbon(ap)
                        */
                    }

                    rvec[0] = 0;
                    rvec[1] = 0;
                    rvec[2] = 0;
                    rand3( &svec[0],&svec[1],&svec[2]);
                    x = gsdg_line_search( svec, &y,ap);
                    rvec[0] += y*svec[0];
                    rvec[1] += y*svec[1];
                    rvec[2] += y*svec[2];
                    rand3( &svec[0],&svec[1],&svec[2]);
                    x = gsdg_line_search( svec, &y,ap);
                    rvec[0] += y*svec[0];
                    rvec[1] += y*svec[1];
                    rvec[2] += y*svec[2];
                    rand3( &svec[0],&svec[1],&svec[2]);
                    x = gsdg_line_search( svec, &y,ap);
                    rvec[0] += y*svec[0];
                    rvec[1] += y*svec[1];
                    rvec[2] += y*svec[2];

                    x = gsdg_line_search( rvec,&y,ap);

                    ap->x += y*rvec[0];
                    ap->y += y*rvec[1];
                    ap->z += y*rvec[2];
                }/* end of active if */
            }/* end of serial number bounds if */

            if( ap == ap->next ) break;
            ap = ap->next;
        }
        printf(" done \n");
    }/* end of iter loop */
}/* end of routine */

float gsdg_line_search( vect, step,who )
float vect[3],*step;
ATOM *who;
{
    float val;
    float vt,lam;
    int i,j;
    float dstep;

    float gsdg_dgeom();

    val = gsdg_dgeom(vect,0.,who);
    lam = 0;
    *step = 0;
    dstep = -.5;
    for( i=0; i< 3; i++)
    {
        dstep *= -.5;
        for( j = 0; j< 200 ; j++)
        {
            lam += dstep;
            vt =  gsdg_dgeom(vect,lam,who);
            if( vt < val ){ *step = lam; val = vt;} else {break;}
        }
        if( j == 200) dstep *= -2;
    }
    return val;
}/*end of routine */

float gsdg_dgeom( vect,lam,who)
ATOM *who;
float vect[3],lam;
{
    int numatm,a_number();
    int i;
    float x,y,z;
    ATOM *ap,*a_next();
    float dt;
    float dsum;

    numatm = a_number();
    x = who->x + vect[0]*lam;
    y = who->y + vect[1]*lam;
    z = who->z + vect[2]*lam;

    dsum = 0.;
    for( i=0; i< numatm; i++)
    {
        ap = a_next(i);
        if( ap != who )
        {
            dt = (x -ap->x)*(x-ap->x);
            dt += (y -ap->y)*(y-ap->y);
            dt += (z -ap->z)*(z-ap->z);

            if( ap->vy > 0 )
            {
                dsum += ap->vy*(ap->vx -dt)*(ap->vx -dt);
            } else {
                if( ap->vx > dt)
                    dsum -= ap->vy*(ap->vx -dt)*(ap->vx -dt);
            }

        }
    }
    return dsum;
}/* end of routine */
/* dgeom trace routines */
int v_trace(V ,lambda )
float *V,lambda;
{
    int numatm,i,a_number();
    ATOM *ap,*a_next();
    float xc,yc,zc;
    float xt,yt,zt;
    float l_trace,get_f_variable();

    numatm = a_number();
    if( numatm < 2 ) return ;
    l_trace = get_f_variable("trace");
    if( l_trace == 0.) l_trace = 1./numatm;

    xc = 0.; yc = 0.; zc = 0.;

    for( i=0; i< numatm; i++)
    {
        ap = a_next(i);
        xc += ap->x + lambda*ap->dx;
        yc += ap->y + lambda*ap->dy;
        zc += ap->z + lambda*ap->dz;
    }
    xc /= numatm; yc /= numatm; zc /= numatm;
    for( i=0; i< numatm; i++)
    {
        ap = a_next(i);
        xt = ap->x + lambda*ap->dx - xc;
        yt = ap->y + lambda*ap->dy - yc;
        zt = ap->z + lambda*ap->dz - zc;
        *V -= l_trace*( xt*xt + yt*yt + zt*zt);
    }

}
int f_trace(lambda )
float lambda;
{
    int numatm,i,a_number();
    ATOM *ap,*a_next();
    float xc,yc,zc;
    float xt,yt,zt;
    float l_trace,get_f_variable();

    numatm = a_number();
    if( numatm < 2 ) return ;
    l_trace = get_f_variable("trace");
    if( l_trace == 0.) l_trace = 1./numatm;

    xc = 0.; yc = 0.; zc = 0.;

    for( i=0; i< numatm; i++)
    {
        ap = a_next(i);
        xc += ap->x + lambda*ap->dx;
        yc += ap->y + lambda*ap->dy;
        zc += ap->z + lambda*ap->dz;
    }
    xc /= numatm; yc /= numatm; zc /= numatm;
    l_trace = 2*l_trace*(1.-1./numatm);
    for( i=0; i< numatm; i++)
    {
        ap = a_next(i);
        xt = ap->x + lambda*ap->dx - xc;
        yt = ap->y + lambda*ap->dy - yc;
        zt = ap->z + lambda*ap->dz - zc;
        /*
        		*V -= l_trace*( xt*xt + yt*yt + zt*zt);
        		ap->fx += 2*l_trace*xt*(1.-1./numatm);
        		ap->fy += 2*l_trace*yt*(1.-1./numatm);
        		ap->fz += 2*l_trace*zt*(1.-1./numatm);
        */
        ap->fx += l_trace*xt;
        ap->fy += l_trace*yt;
        ap->fz += l_trace*zt;
    }

}
