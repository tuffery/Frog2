/* gsdx.c
* just find a feasible solution. 
*
*  instead of the error being minimized as a target
*  (i.e. points on spheres)
*  minimize the number of constraints being used
*  in the L1 sense
*   The error is the violations of the inequalities
*    from the constraint graph
*    ( 1,..,-1,...) <= v 
*    where the graph matix is one for the atom and -1 for its
*    constraint target.
*  x,y,z are searched at once and the error is evaluated in a 
*   3-d spherical sense.
*  This then becomes the Bellman-Ford graph algorithm.
*/
/* bellman.c
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

void bellman( vfs,nfs,niter,low_serial,high_serial )
int (*vfs[])(),nfs;
int niter;
int low_serial,high_serial;
{
    ATOM *ap,*bp,*a_next();
    int numatm , a_number();
    int iter, i,j,k;
    float svec[3],rvec[3]; /* storage for search vectors */
    float x,y,z;
    float bellman_line_search();
    int v_bond(),v_angle(),u_v_nonbon(),v_noel();
    int v_mmbond(),v_mmangle(),v_c_angle(),v_nonbon();
    int v_hard(),v_abc();
    int v_hybrid();
    int local_iter;

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
                /*if( ap->active)*/{
                    for( local_iter=0; local_iter<1; local_iter++)
                    {

                        for( j=0;j< numatm; j++)
                        { bp = a_next(j);
                            bp->vx = 16.; /* default to 4A separation */
                            bp->vx = 9.; /* try 3 */
                            bp->vy = 0.; /* but only as a lower bound */
                        }

                        for( j=0; j< nfs; j++)
                        {
                            if( vfs[j] == v_nonbon || vfs[j] == u_v_nonbon)
                                { for( k=0; k< numatm; k++)
                                {bp = a_next(k);
                                    bp->vy = -10;}
                                break;
                            }
                        }


                        for( j=0; j< nfs; j++)
                        {
                            if( vfs[j] == v_bond || vfs[j] == v_mmbond || vfs[j] == v_abc) gsdg_bond(ap);
                            if( vfs[j] == v_angle || vfs[j] == v_mmangle || vfs[j] == v_c_angle) gsdg_angle(ap);
                            if( vfs[j] == v_noel) gsdg_noel(ap);
                            if( vfs[j] == v_hybrid) gsdg_hybrid(ap);
                        }
                        ap->vy = 0.;

                        bellman_relax(ap);

                    }/* end of local_iter */
                }/* end of active if */
            }/* end of serial number bounds if */

            if( ap == ap->next ) break;
            ap = ap->next;
        }
        printf(" done \n");
    }/* end of iter loop */
}/* end of routine */

int bellman_relax(who)
ATOM *who;
{
    int na,a_mumber();
    int i;
    ATOM *ap,*a_next();
    float r,dx,dy,dz,target;
    /*	float delta ; */

    na = a_number();
    if( na < 2) return 0;

    for(i=0; i< na; i++)
    {ap = a_next(i);
        if( ap->active ){
            dx = ap->x - who->x;
            dy = ap->y - who->y;
            dz = ap->z - who->z;
            r = dx*dx + dy*dy + dz*dz;

            if( r < 1.e-5) {
                rand3( &dx,&dy,&dz);
                r = dx*dx + dy*dy + dz*dz;
            }

            if( ap->vy > 0. && r > ap->vx ){
                target = sqrt(ap->vx/r);
                if( target < 0.75) target = 0.75;
                /*
                //		who->x = 0.5*(who->x + ap->x - dx*target);
                //		who->y = 0.5*(who->y + ap->y - dy*target);
                //		who->z = 0.5*(who->z + ap->z - dz*target);
                */
                /*
                		ap->x = 0.25*( 3*ap->x + who->x + dx*target);
                		ap->y = 0.25*( 3*ap->y + who->y + dy*target);
                		ap->z = 0.25*( 3*ap->z + who->z + dz*target);
                */
                ap->x = 0.5*( ap->x + who->x + dx*target);
                ap->y = 0.5*( ap->y + who->y + dy*target);
                ap->z = 0.5*( ap->z + who->z + dz*target);

            }

            if( ap->vy < 0. && r < ap->vx){
                target = sqrt(ap->vx/r);

                /*
                //		who->x = 0.5*(who->x + ap->x - dx*target);
                //		who->y = 0.5*(who->y + ap->y - dy*target);
                //		who->z = 0.5*(who->z + ap->z - dz*target);
                */

                ap->x = 0.25*(3*ap->x + who->x + dx*target);
                ap->y = 0.25*(3*ap->y + who->y + dy*target);
                ap->z = 0.25*(3*ap->z + who->z + dz*target);

            }
        }/* ap->active ?*/
    }/* i */

    return 0;
}/* end of bellman_relax */

float bellman_line_search( vect, step,who )
float vect[3],*step;
ATOM *who;
{
    float val;
    float vt,lam;
    int i,j;
    float dstep;

    float bellman_dgeom();

    val = bellman_dgeom(vect,0.,who);
    lam = 0;
    *step = 0;
    dstep = -.5;
    dstep = -1.;
    for( i=0; i< 3; i++)
    {
        dstep *= -.5;
        for( j = 0; j< 200 ; j++)
        {
            lam += dstep;
            vt =  bellman_dgeom(vect,lam,who);
            if( vt < val ){ *step = lam; val = vt;} else {break;}
        }
        if( j == 200) dstep *= -2;
    }
    return val;
}/*end of routine */

float bellman_dgeom( vect,lam,who)
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



            if( ap->vy > 0 )
            {
                /*
                dsum += ap->vy*(ap->vx -dt)*(ap->vx -dt);
                */
                dt = (x -ap->x)*(x-ap->x);
                dt += (y -ap->y)*(y-ap->y);
                dt += (z -ap->z)*(z-ap->z);

                if( dt > ap->vx) dsum += fabs(dt - ap->vx);
            } else {
                /*
                if( ap->vx > dt)
                dsum -= ap->vy*(ap->vx -dt)*(ap->vx -dt);
                */
                dt = (x -ap->x)*(x-ap->x);
                dt += (y -ap->y)*(y-ap->y);
                dt += (z -ap->z)*(z-ap->z);

                if( dt < ap->vx) dsum += fabs(ap->vx - dt);

            }

        }
    }
    return dsum;
}/* end of routine */
