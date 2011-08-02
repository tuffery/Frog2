/* angle.c
*
* collection of routines to service bond angle  potentials
*
* POOP (Poor-mans Object Oriented Programming) using scope rules
*
* these routines hold a data base (in terms of array indeces)
* of angle, with the associated length and force constant
*
* (this could be table driven but what the hell memories cheap)
*
* the routines for potential value, force and (eventually) second
* derivatives are here also
*
* force and 2nd derivative routines assume zero'd arrays for output
* this allows for parralellization if needed (on a PC?)
*
* forces are angle wise symmetric - so we don't have to fuck around with
* s matrices and the like.
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
#ifdef ANSI
#include <stdlib.h>
#endif

#include "ammp.h"
/* ATOM structure contains a serial number for indexing into
* arrays and the like (a Hessian)
* but otherwise is self-contained. Note the hooks for Non-angleed potentials
*/
typedef struct{
    ATOM *atom1,*atom2,*atom3;
    float target,k;
    void *next;
}  ANGLE;
#define ALONG sizeof(ANGLE)

ANGLE *angle_first = NULL;
ANGLE *angle_last = NULL;
/* function angle adds a angle to the angle list
* returns 1 if ok
* returns 0 if not
*  is passed the array pointers, length and constant
* allocates the new memory, initializes it and
* returns
*/
int angle( p1,p2,p3,fk,bl)
int p1,p2,p3;
float bl,fk ;
{
    ANGLE *new;
    ATOM *ap1,*ap2,*ap3,*a_m_serial();
    char line[BUFSIZ];
    int i;
    /* get the atom pointers for the two serial numbers */
    ap1 = a_m_serial( p1 );
    ap2 = a_m_serial( p2 );
    ap3 = a_m_serial( p3 );
    if( (ap1 == NULL) || (ap2 == NULL) || (ap3==NULL) )
    {
        sprintf( line,"undefined atom in angle %d %d %d \0",p1,p2,p3);
        aaerror( line );
        return 0;
    }

    if( ( new = malloc( ALONG ) ) == NULL)
    {
        return 0;
    }
    /* initialize the pointers */
    if( angle_first == NULL) angle_first = new;
    if( angle_last == NULL) angle_last = new;
    new -> atom1 = ap1;
    new -> atom2 = ap2;
    new -> atom3 = ap3;
    new -> target = bl;
    new -> k = fk;
    new -> next = new;
    /* update the exclude list in the atoms structure */
    if( ap1->dontuse < NEXCLUDE)
    {
        for( i=0; i< ap1->dontuse; i++)
        {
            if( ap1->excluded[i] == ap3) goto excluded1;
        }
        ap1->excluded[ap1->dontuse] = ap3; (ap1->dontuse)++;
    }else{
        aaerror(" too many bonds to an atom increase NEXCLUDE in ammp.h");
        exit(0);
    }
excluded1:
    if( ap3->dontuse < NEXCLUDE)
    {
        for( i=0; i< ap3->dontuse; i++)
        {
            if( ap3->excluded[i] == ap1) goto excluded3;
        }
        ap3->excluded[ap3->dontuse] = ap1; (ap3->dontuse)++;
    }else{
        aaerror(" too many bonds to an atom increase NEXCLUDE in ammp.h");
        exit(0);
    }
excluded3:
    angle_last -> next = new;
    angle_last = new;
    return 1;
}


/* v_angle()
* this function sums up the potentials
* for the atoms defined in the angle data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int v_angle( V, lambda )
float *V,lambda;
{
    ANGLE *bp;
    float r,x1,y1,z1,x2,y2,z2;
    float dp;
    ATOM *a1,*a2,*a3;


    bp = angle_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1; a2 = bp->atom2; a3 = bp->atom3;
        if( a1->active || a2->active || a3->active )
        {
            x1 = (a1->x -a2->x +lambda*(a1->dx-a2->dx));
            y1 = (a1->y -a2->y +lambda*(a1->dy-a2->dy));
            z1 = (a1->z -a2->z +lambda*(a1->dz-a2->dz));
            x2 = (a3->x -a2->x +lambda*(a3->dx-a2->dx));
            y2 = (a3->y -a2->y +lambda*(a3->dy-a2->dy));
            z2 = (a3->z -a2->z +lambda*(a3->dz-a2->dz));
            dp = x1*x2+y1*y2+z1*z2;
            r = (( x1*x1+y1*y1+z1*z1)*(x2*x2+y2*y2+z2*z2));
            if( r > 1.e-8){
                r = sqrt(r);
                dp = dp/r;  if( dp > 1.) dp = 1.; if( dp < -1.) dp = -1.;
                dp = acos(dp);
                *V += bp->k * (bp->target-dp)*(bp->target-dp);
            }
        }
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* f_angle()
*
* f_angle increments the forces in the atom structures by the force
* due to the angle components.  NOTE THE WORD increment.
* the forces should first be zero'd.
* if not then this code will be invalid.  THIS IS DELIBERATE.
* on bigger (and better?) machines the different potential terms
* may be updated at random or in parrellel, if we assume that this routine
* will initialize the forces then we can't do this.
*/
int f_angle(lambda)
float lambda;
/*  returns 0 if error, 1 if OK */
{
    ANGLE *bp;
    float r,k,ux1,uy1,uz1,ux2,uy2,uz2;
    ATOM *a1,*a2,*a3;
    float x1,y1,z1,x2,y2,z2;
    float r1,r2,dtheta,dp;
    float r11,r22,sdth;

    bp = angle_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        k = bp->k;
        a1 = bp->atom1; a2 = bp->atom2; a3 = bp->atom3;
        if( a1->active|| a2->active || a3->active )
        {
            x1 = (a1->x -a2->x +lambda*(a1->dx-a2->dx));
            y1 = (a1->y -a2->y +lambda*(a1->dy-a2->dy));
            z1 = (a1->z -a2->z +lambda*(a1->dz-a2->dz));
            x2 = (a3->x -a2->x +lambda*(a3->dx-a2->dx));
            y2 = (a3->y -a2->y +lambda*(a3->dy-a2->dy));
            z2 = (a3->z -a2->z +lambda*(a3->dz-a2->dz));
            dp = x1*x2+y1*y2+z1*z2;
            r1 = sqrt(x1*x1+y1*y1+z1*z1);
            r2 = sqrt(x2*x2+y2*y2+z2*z2);
            if( r1 < 1.e-5 || r2 < 1.e-5) goto SKIP;
            r = r1*r2;
            if( r > 1.e-8){

                dp = dp/r;  if( dp > 1.) dp = 1.; if( dp < -1.) dp = -1.;
                dtheta = acos(dp);
                sdth = sin(dtheta); if( sdth < 1.e-3) sdth = 1.e-3;
                r11 = r1*sdth; r22 = r2*sdth;
                ux1 = x2/r2 - dp*x1/r1;
                uy1 = y2/r2 - dp*y1/r1;
                uz1 = z2/r2 - dp*z1/r1;
                ux2 = x1/r1 - dp*x2/r2;
                uy2 = y1/r1 - dp*y2/r2;
                uz2 = z1/r1 - dp*z2/r2;
                dtheta = -2.*k*(bp->target - dtheta);
                ux1 = ux1*dtheta/r11;
                uy1 = uy1*dtheta/r11;
                uz1 = uz1*dtheta/r11;
                ux2 = ux2*dtheta/r22;
                uy2 = uy2*dtheta/r22;
                uz2 = uz2*dtheta/r22;
                if( a1->active)
                {
                    a1->fx += ux1;
                    a1->fy += uy1;
                    a1->fz += uz1;
                }

                if( a2->active)
                {
                    a2->fx += -ux1 - ux2;
                    a2->fy += -uy1 - uy2;
                    a2->fz += -uz1 - uz2;
                }

                if( a3->active)
                {
                    a3->fx += ux2;
                    a3->fy += uy2;
                    a3->fz += uz2;
                }

            }
        }
SKIP:
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* function get_angle( a1,bonded,10,inbond);
* check the ANGLE list for atoms 1-3 ed to a1
*/
void get_angle( a1,bonded,mbond,inbond)
ATOM *a1, *bonded[];
int mbond,*inbond ;
{
    ANGLE *mine;
    mine = angle_first;
    *inbond = 0;
    while(1)
    {
        if( (mine == NULL) )
        {
            return;
        }
        if( mine->atom1 == a1)
        {
            bonded[(*inbond)++] = mine->atom3;
        }
        if( mine->atom3 == a1)
        {
            bonded[(*inbond)++] = mine->atom1;
        }
        if( mine == mine->next) return;
        mine = mine->next;
        if( *inbond == mbond ) return;
    }
}
/* routine dump_angles
* this function outputs the angle parameters
* and does it in a simple form
* angle ser1,ser2,ser3,k,theta (in degrees )
* the rest is just free format
*/
void dump_angles( where )
FILE *where;
{
    ANGLE *b;
    ATOM *a1,*a2,*a3;
    float rtodeg;
    b = angle_first;
    if( b == NULL ) return;
    rtodeg = 180./acos(-1.);
    while( (b->next != b) )
    {
        if( b->next == NULL) return;
        a1 = b->atom1; a2 = b->atom2;a3 = b->atom3;
        fprintf( where,"angle %d %d %d %f %f \;\n",a1->serial,a2->serial,
                 a3-> serial,b->k,b->target*rtodeg);
        b = b->next;
    }
    if( b->next == NULL) return;
    a1 = b->atom1; a2 = b->atom2;a3 = b->atom3;
    fprintf( where,"angle %d %d %d %f %f \;\n",a1->serial,a2->serial,
             a3-> serial,b->k,b->target*rtodeg);
}

/* MM3 angle function  */

/* v_mmangle()
* this function sums up the potentials
* for the atoms defined in the angle data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int v_mmangle( V, lambda )
float *V,lambda;
{
    ANGLE *bp;
    float r,x1,y1,z1,x2,y2,z2;
    float dp;
    ATOM *a1,*a2,*a3;


    bp = angle_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1; a2 = bp->atom2; a3 = bp->atom3;
        if( a1->active || a2->active || a3->active )
        {
            x1 = (a1->x -a2->x +lambda*(a1->dx-a2->dx));
            y1 = (a1->y -a2->y +lambda*(a1->dy-a2->dy));
            z1 = (a1->z -a2->z +lambda*(a1->dz-a2->dz));
            x2 = (a3->x -a2->x +lambda*(a3->dx-a2->dx));
            y2 = (a3->y -a2->y +lambda*(a3->dy-a2->dy));
            z2 = (a3->z -a2->z +lambda*(a3->dz-a2->dz));
            dp = x1*x2+y1*y2+z1*z2;
            r = (( x1*x1+y1*y1+z1*z1)*(x2*x2+y2*y2+z2*z2));
            if( r > 1.e-8){
                r = sqrt(r);
                dp = dp/r;  if( dp > 1.) dp = 1.; if( dp < -1.) dp = -1.;
                dp = acos(dp);
                /*	*V += bp->k * (bp->target-dp)*(bp->target-dp);
                */
                dp = dp - bp->target;
                *V += bp->k*dp*dp*(1.-.014*dp+5.6e-5*dp*dp
                                   -7.e-7*dp*dp*dp +9e-10*dp*dp*dp*dp);
            }
        }
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* f_mmangle()
*
* f_mmangle increments the forces in the atom structures by the force
* due to the angle components.  NOTE THE WORD increment.
* the forces should first be zero'd.
* if not then this code will be invalid.  THIS IS DELIBERATE.
* on bigger (and better?) machines the different potential terms
* may be updated at random or in parrellel, if we assume that this routine
* will initialize the forces then we can't do this.
*/
int f_mmangle(lambda)
float lambda;
/*  returns 0 if error, 1 if OK */
{
    ANGLE *bp;
    float r,k,ux1,uy1,uz1,ux2,uy2,uz2;
    ATOM *a1,*a2,*a3;
    float x1,y1,z1,x2,y2,z2;
    float r1,r2,dtheta,dp;
    float r11,r22,sdth;

    bp = angle_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        k = bp->k;
        a1 = bp->atom1; a2 = bp->atom2; a3 = bp->atom3;
        if( a1->active|| a2->active || a3->active )
        {
            x1 = (a1->x -a2->x +lambda*(a1->dx-a2->dx));
            y1 = (a1->y -a2->y +lambda*(a1->dy-a2->dy));
            z1 = (a1->z -a2->z +lambda*(a1->dz-a2->dz));
            x2 = (a3->x -a2->x +lambda*(a3->dx-a2->dx));
            y2 = (a3->y -a2->y +lambda*(a3->dy-a2->dy));
            z2 = (a3->z -a2->z +lambda*(a3->dz-a2->dz));
            dp = x1*x2+y1*y2+z1*z2;
            r1 = sqrt(x1*x1+y1*y1+z1*z1);
            r2 = sqrt(x2*x2+y2*y2+z2*z2);
            if( r1 < 1.e-5 || r2 < 1.e-5) goto SKIP;
            r = r1*r2;
            if( r > 1.e-8){

                dp = dp/r;  if( dp > 1.) dp = 1.; if( dp < -1.) dp = -1.;
                dtheta = acos(dp);
                sdth = sin(dtheta); if( sdth < 1.e-3) sdth = 1.e-3;
                r11 = r1*sdth; r22 = r2*sdth;
                ux1 = x2/r2 - dp*x1/r1;
                uy1 = y2/r2 - dp*y1/r1;
                uz1 = z2/r2 - dp*z1/r1;
                ux2 = x1/r1 - dp*x2/r2;
                uy2 = y1/r1 - dp*y2/r2;
                uz2 = z1/r1 - dp*z2/r2;
                /*	*V += bp->k*dp*dp*(1.-.014*dp+5.6e-5*dp*dp
                	-7.e-7*dp*dp*dp +9e-10*dp*dp*dp*dp);
                */
                dp = dtheta - bp->target;
                dtheta = k*dp*(2.-.014*3*dp+4*(5.6e-5)*dp*dp
                               -5*(7.e-7)*dp*dp*dp +6*(9.e-10)*dp*dp*dp*dp);
                ux1 = ux1*dtheta/r11;
                uy1 = uy1*dtheta/r11;
                uz1 = uz1*dtheta/r11;
                ux2 = ux2*dtheta/r22;
                uy2 = uy2*dtheta/r22;
                uz2 = uz2*dtheta/r22;
                if( a1->active){
                    a1->fx += ux1;
                    a1->fy += uy1;
                    a1->fz += uz1;
                }

                if( a2->active){
                    a2->fx += -ux1 - ux2;
                    a2->fy += -uy1 - uy2;
                    a2->fz += -uz1 - uz2;
                }

                if( a3->active){
                    a3->fx += ux2;
                    a3->fy += uy2;
                    a3->fz += uz2;
                }

            }
        }
SKIP:
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* a_angle()
* this function sums up the potentials
* for the atoms defined in the angle data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int a_angle( V, lambda, ilow,ihigh,op )
float *V,lambda;
int ilow,ihigh;
FILE *op;
{
    ANGLE *bp;
    float r,x1,y1,z1,x2,y2,z2;
    float dp;
    ATOM *a1,*a2,*a3;


    bp = angle_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1; a2 = bp->atom2; a3 = bp->atom3;
        if( (a1->serial >= ilow && a1->serial <= ihigh)
                ||  (a2->serial >= ilow && a2->serial <= ihigh)
                ||  (a3->serial >= ilow && a3->serial <= ihigh) )
        {
            x1 = (a1->x -a2->x +lambda*(a1->dx-a2->dx));
            y1 = (a1->y -a2->y +lambda*(a1->dy-a2->dy));
            z1 = (a1->z -a2->z +lambda*(a1->dz-a2->dz));
            x2 = (a3->x -a2->x +lambda*(a3->dx-a2->dx));
            y2 = (a3->y -a2->y +lambda*(a3->dy-a2->dy));
            z2 = (a3->z -a2->z +lambda*(a3->dz-a2->dz));
            dp = x1*x2+y1*y2+z1*z2;
            r = (( x1*x1+y1*y1+z1*z1)*(x2*x2+y2*y2+z2*z2));
            if( r > 1.e-8){
                r = sqrt(r);
                dp = dp/r;  if( dp > 1.) dp = 1.; if( dp < -1.) dp = -1.;
                dp = acos(dp);
                z2  = bp->k * (bp->target-dp)*(bp->target-dp);
                *V += z2;
                fprintf(op,"Angle %s %d %s %d %s %d E %f value %f error %f\n",
                        a1->name,a1->serial,a2->name,a2->serial,a3->name,a3->serial
                        ,z2,dp*180./3.14159265,
                        (dp-bp->target)*180./3.14159265);
            }
        }
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* MM3 angle function  */

/* a_mmangle()
* this function sums up the potentials
* for the atoms defined in the angle data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int a_mmangle( V, lambda ,ilow,ihigh,op)
float *V,lambda;
int ilow,ihigh;
FILE *op;
{
    ANGLE *bp;
    float r,x1,y1,z1,x2,y2,z2;
    float dp;
    ATOM *a1,*a2,*a3;


    bp = angle_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1; a2 = bp->atom2; a3 = bp->atom3;
        if( (a1->serial >= ilow && a1->serial <= ihigh)
                ||  (a2->serial >= ilow && a2->serial <= ihigh)
                ||  (a3->serial >= ilow && a3->serial <= ihigh) )
        {
            x1 = (a1->x -a2->x +lambda*(a1->dx-a2->dx));
            y1 = (a1->y -a2->y +lambda*(a1->dy-a2->dy));
            z1 = (a1->z -a2->z +lambda*(a1->dz-a2->dz));
            x2 = (a3->x -a2->x +lambda*(a3->dx-a2->dx));
            y2 = (a3->y -a2->y +lambda*(a3->dy-a2->dy));
            z2 = (a3->z -a2->z +lambda*(a3->dz-a2->dz));
            dp = x1*x2+y1*y2+z1*z2;
            r = (( x1*x1+y1*y1+z1*z1)*(x2*x2+y2*y2+z2*z2));
            if( r > 1.e-8){
                r = sqrt(r);
                dp = dp/r;  if( dp > 1.) dp = 1.; if( dp < -1.) dp = -1.;
                dp = acos(dp);
                /*	*V += bp->k * (bp->target-dp)*(bp->target-dp);
                */
                dp = dp - bp->target;
                z2 = bp->k*dp*dp*(1.-.014*dp+5.6e-5*dp*dp
                                  -7.e-7*dp*dp*dp +9e-10*dp*dp*dp*dp);
                *V += z2;
                dp = dp + bp->target;
                fprintf(op,"mmAngle %s %d %s %d %s %d E %f value %f error %f\n",
                        a1->name,a1->serial,a2->name,a2->serial,a3->name,a3->serial
                        ,z2,dp*180./3.14159265,
                        (dp-bp->target)*180./3.14159265);
                /*
                	fprintf(op,"mm Angle %d %d %d E %f value %f error %f\n",
                		a1->serial,a2->serial,a3->serial,z2,dp*180./3.14159265,
                		(dp-bp->target)*180./3.14159265);
                */
                /*	fprintf(op,"mm Angle %d %d %d  E %f error %f\n",
                		a1->serial,a2->serial,a3->serial,z2,dp*180./3.14159265);
                */
            }
        }
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* UFF cos angle types */
/* v_angle()
* this function sums up the potentials
* for the atoms defined in the angle data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int v_c_angle( V, lambda )
float *V,lambda;
{
    ANGLE *bp;
    float r,x1,y1,z1,x2,y2,z2;
    float dp;
    ATOM *a1,*a2,*a3;
    float C0,C1,C2;


    bp = angle_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1; a2 = bp->atom2; a3 = bp->atom3;
        if( a1->active || a2->active || a3->active )
        {
            x1 = (a1->x -a2->x +lambda*(a1->dx-a2->dx));
            y1 = (a1->y -a2->y +lambda*(a1->dy-a2->dy));
            z1 = (a1->z -a2->z +lambda*(a1->dz-a2->dz));
            x2 = (a3->x -a2->x +lambda*(a3->dx-a2->dx));
            y2 = (a3->y -a2->y +lambda*(a3->dy-a2->dy));
            z2 = (a3->z -a2->z +lambda*(a3->dz-a2->dz));
            dp = x1*x2+y1*y2+z1*z2;
            r = (( x1*x1+y1*y1+z1*z1)*(x2*x2+y2*y2+z2*z2));
            if( r > 1.e-8){
                r = sqrt(r);
                dp = dp/r;  if( dp > 1.) dp = 1.; if( dp < -1.) dp = -1.;
                r = dp;
                dp = acos(dp);
                /*	*V += bp->k * (bp->target-dp)*(bp->target-dp); */
                C0 = cos( bp->target);
                C2 = 1./(4. - 4*C0*C0);
                C1 = -4.*C2*C0;
                C0 = C2*(2*C0*C0 + 1);
                /* the 2* is extra because of unit differences between my quadratic and uff*/
                *V += 2*bp->k *( C0 + C1*r + C2*cos(dp*2));
            }
        }
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* f_angle()
*
* f_angle increments the forces in the atom structures by the force
* due to the angle components.  NOTE THE WORD increment.
* the forces should first be zero'd.
* if not then this code will be invalid.  THIS IS DELIBERATE.
* on bigger (and better?) machines the different potential terms
* may be updated at random or in parrellel, if we assume that this routine
* will initialize the forces then we can't do this.
*/
int f_c_angle(lambda)
float lambda;
/*  returns 0 if error, 1 if OK */
{
    ANGLE *bp;
    float r,k,ux1,uy1,uz1,ux2,uy2,uz2;
    ATOM *a1,*a2,*a3;
    float x1,y1,z1,x2,y2,z2;
    float r1,r2,dtheta,dp;
    float r11,r22,sdth;
    float C0,C1,C2;

    bp = angle_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        k = bp->k;
        a1 = bp->atom1; a2 = bp->atom2; a3 = bp->atom3;
        if( a1->active|| a2->active || a3->active )
        {
            x1 = (a1->x -a2->x +lambda*(a1->dx-a2->dx));
            y1 = (a1->y -a2->y +lambda*(a1->dy-a2->dy));
            z1 = (a1->z -a2->z +lambda*(a1->dz-a2->dz));
            x2 = (a3->x -a2->x +lambda*(a3->dx-a2->dx));
            y2 = (a3->y -a2->y +lambda*(a3->dy-a2->dy));
            z2 = (a3->z -a2->z +lambda*(a3->dz-a2->dz));
            dp = x1*x2+y1*y2+z1*z2;
            r1 = sqrt(x1*x1+y1*y1+z1*z1);
            r2 = sqrt(x2*x2+y2*y2+z2*z2);
            if( r1 < 1.e-5 || r2 < 1.e-5) goto SKIP;
            r = r1*r2;
            if( r > 1.e-8){

                dp = dp/r;  if( dp > 1.) dp = 1.; if( dp < -1.) dp = -1.;
                dtheta = acos(dp);
                sdth = sin(dtheta); if( sdth < 1.e-3) sdth = 1.e-3;
                r11 = r1*sdth; r22 = r2*sdth;
                ux1 = x2/r2 - dp*x1/r1;
                uy1 = y2/r2 - dp*y1/r1;
                uz1 = z2/r2 - dp*z1/r1;
                ux2 = x1/r1 - dp*x2/r2;
                uy2 = y1/r1 - dp*y2/r2;
                uz2 = z1/r1 - dp*z2/r2;
                C0 = cos( bp->target);
                C2 = 1./(4. - 4*C0*C0);
                C1 = -4.*C2*C0;
                C0 = C2*(2*C0*C0 + 1);
                /*	*V += bp->k *( C0 + C1*r + C2*cos(dp*2)); */
                dtheta = -2.*bp->k*(C1*sdth + 2*C2*sin(dtheta*2) );
                /*	dtheta = -2.*k*(bp->target - dtheta); */
                ux1 = ux1*dtheta/r11;
                uy1 = uy1*dtheta/r11;
                uz1 = uz1*dtheta/r11;
                ux2 = ux2*dtheta/r22;
                uy2 = uy2*dtheta/r22;
                uz2 = uz2*dtheta/r22;
                if( a1->active){
                    a1->fx += ux1;
                    a1->fy += uy1;
                    a1->fz += uz1;
                }

                if( a2->active){
                    a2->fx += -ux1 - ux2;
                    a2->fy += -uy1 - uy2;
                    a2->fz += -uz1 - uz2;
                }

                if( a3->active){
                    a3->fx += ux2;
                    a3->fy += uy2;
                    a3->fz += uz2;
                }

            }
        }
SKIP:
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
int a_c_angle( V, lambda, ilow,ihigh,op )
float *V,lambda;
int ilow,ihigh;
FILE *op;
{
    ANGLE *bp;
    float r,x1,y1,z1,x2,y2,z2;
    float dp;
    ATOM *a1,*a2,*a3;
    float C0,C1,C2;


    bp = angle_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1; a2 = bp->atom2; a3 = bp->atom3;
        if( (a1->serial >= ilow && a1->serial <= ihigh)
                ||  (a2->serial >= ilow && a2->serial <= ihigh)
                ||  (a3->serial >= ilow && a3->serial <= ihigh) )
        {
            x1 = (a1->x -a2->x +lambda*(a1->dx-a2->dx));
            y1 = (a1->y -a2->y +lambda*(a1->dy-a2->dy));
            z1 = (a1->z -a2->z +lambda*(a1->dz-a2->dz));
            x2 = (a3->x -a2->x +lambda*(a3->dx-a2->dx));
            y2 = (a3->y -a2->y +lambda*(a3->dy-a2->dy));
            z2 = (a3->z -a2->z +lambda*(a3->dz-a2->dz));
            dp = x1*x2+y1*y2+z1*z2;
            r = (( x1*x1+y1*y1+z1*z1)*(x2*x2+y2*y2+z2*z2));
            if( r > 1.e-8){
                r = sqrt(r);
                dp = dp/r;  if( dp > 1.) dp = 1.; if( dp < -1.) dp = -1.;
                r = dp;
                dp = acos(dp);
                /*	z2  = bp->k * (bp->target-dp)*(bp->target-dp);
                	*V += z2;
                */
                C0 = cos( bp->target);
                C2 = 1./(4. - 4*C0*C0);
                C1 = -4.*C2*C0;
                C0 = C2*(2*C0*C0 + 1);
                /* the 2* is extra because of unit differences between my quadratic and uff*/
                *V += 2*bp->k *( C0 + C1*r + C2*cos(dp*2));
                fprintf(op,"c Angle %s %d %s %d %s %d E %f value %f error %f\n",
                        a1->name,a1->serial,a2->name,a2->serial,a3->name,a3->serial
                        ,z2,dp*180./3.14159265,
                        (dp-bp->target)*180./3.14159265);
                /*
                	fprintf(op,"Angle %d %d %d E %f value %f error %f\n",
                		a1->serial,a2->serial,a3->serial,z2,dp*180./3.14159265,
                		(dp-bp->target)*180./3.14159265);
                */
            }
        }
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
gsdg_angle( who )
ATOM *who;
{
    ANGLE *bp;
    ATOM *ap,*ap1;
    float r,r1,theta;
    float bond_length();

    bp = angle_first;

    while( 1)
    { if( bp == NULL ) return ;
        if( bp->atom1 == who ){
            ap1 = bp->atom2;
            ap = bp->atom3;
            r = bond_length(who,ap1); r1 = bond_length(ap1,ap);
            theta = r*r + r1*r1 -2*cos( bp->target) *r*r1;
            ap->vx = theta;
            ap->vy = bp->k;
        }
        if( bp->atom3 == who ){
            ap1 = bp->atom2;
            ap = bp->atom1;
            r = bond_length(who,ap1); r1 = bond_length(ap1,ap);
            theta = r*r + r1*r1 -2*cos( bp->target) *r*r1;
            ap->vx = theta;
            ap->vy = bp->k;
        }
        if( bp == bp->next) return;
        bp = bp->next;
    }
}
/* v_ho_angle()
* this function sums up the potentials
* for the atoms defined in the angle data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int v_ho_angle( V, lambda )
float *V,lambda;
{
    ANGLE *bp;
    float r,x1,y1,z1,x2,y2,z2;
    float dp;
    ATOM *a1,*a2,*a3;
    float hol,get_f_variable(),target;

    hol = get_f_variable("lambda");
    if( hol < 0. ) hol = 0.;
    if( hol > 1. ) hol = 1.;

    bp = angle_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1; a2 = bp->atom2; a3 = bp->atom3;
        if( a1->active || a2->active || a3->active )
        {
            x1 = (a1->x -a2->x +lambda*(a1->dx-a2->dx));
            y1 = (a1->y -a2->y +lambda*(a1->dy-a2->dy));
            z1 = (a1->z -a2->z +lambda*(a1->dz-a2->dz));
            x2 = (a3->x -a2->x +lambda*(a3->dx-a2->dx));
            y2 = (a3->y -a2->y +lambda*(a3->dy-a2->dy));
            z2 = (a3->z -a2->z +lambda*(a3->dz-a2->dz));
            dp = x1*x2+y1*y2+z1*z2;
            r = (( x1*x1+y1*y1+z1*z1)*(x2*x2+y2*y2+z2*z2));
            if( r > 1.e-8){
                r = sqrt(r);
                dp = dp/r;  if( dp > 1.) dp = 1.; if( dp < -1.) dp = -1.;
                dp = acos(dp);
                target = hol*dp + (1.-hol)*bp->target;
                *V += bp->k * (target-dp)*(target-dp);
            }
        }
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* f_ho_angle()
*
* f_ho_angle increments the forces in the atom structures by the force
* due to the angle components.  NOTE THE WORD increment.
* the forces should first be zero'd.
* if not then this code will be invalid.  THIS IS DELIBERATE.
* on bigger (and better?) machines the different potential terms
* may be updated at random or in parrellel, if we assume that this routine
* will initialize the forces then we can't do this.
*/
int f_ho_angle(lambda)
float lambda;
/*  returns 0 if error, 1 if OK */
{
    ANGLE *bp;
    float r,k,ux1,uy1,uz1,ux2,uy2,uz2;
    ATOM *a1,*a2,*a3;
    float x1,y1,z1,x2,y2,z2;
    float r1,r2,dtheta,dp;
    float r11,r22,sdth;
    float hol,get_f_variable(),target;

    hol = get_f_variable("lambda");
    if( hol < 0. ) hol = 0.;
    if( hol > 1. ) hol = 1.;
    bp = angle_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        k = bp->k;
        a1 = bp->atom1; a2 = bp->atom2; a3 = bp->atom3;
        if( a1->active|| a2->active || a3->active )
        {
            x1 = (a1->x -a2->x +lambda*(a1->dx-a2->dx));
            y1 = (a1->y -a2->y +lambda*(a1->dy-a2->dy));
            z1 = (a1->z -a2->z +lambda*(a1->dz-a2->dz));
            x2 = (a3->x -a2->x +lambda*(a3->dx-a2->dx));
            y2 = (a3->y -a2->y +lambda*(a3->dy-a2->dy));
            z2 = (a3->z -a2->z +lambda*(a3->dz-a2->dz));
            dp = x1*x2+y1*y2+z1*z2;
            r1 = sqrt(x1*x1+y1*y1+z1*z1);
            r2 = sqrt(x2*x2+y2*y2+z2*z2);
            if( r1 < 1.e-5 || r2 < 1.e-5) goto SKIP;
            r = r1*r2;
            if( r > 1.e-8){

                dp = dp/r;  if( dp > 1.) dp = 1.; if( dp < -1.) dp = -1.;
                dtheta = acos(dp);
                target = hol*dtheta + (1.-hol)*bp->target;
                sdth = sin(dtheta); if( sdth < 1.e-3) sdth = 1.e-3;
                r11 = r1*sdth; r22 = r2*sdth;
                ux1 = x2/r2 - dp*x1/r1;
                uy1 = y2/r2 - dp*y1/r1;
                uz1 = z2/r2 - dp*z1/r1;
                ux2 = x1/r1 - dp*x2/r2;
                uy2 = y1/r1 - dp*y2/r2;
                uz2 = z1/r1 - dp*z2/r2;
                dtheta = -2.*k*(target - dtheta)*(1.-hol);
                ux1 = ux1*dtheta/r11;
                uy1 = uy1*dtheta/r11;
                uz1 = uz1*dtheta/r11;
                ux2 = ux2*dtheta/r22;
                uy2 = uy2*dtheta/r22;
                uz2 = uz2*dtheta/r22;
                if( a1->active)
                {
                    a1->fx += ux1;
                    a1->fy += uy1;
                    a1->fz += uz1;
                }

                if( a2->active)
                {
                    a2->fx += -ux1 - ux2;
                    a2->fy += -uy1 - uy2;
                    a2->fz += -uz1 - uz2;
                }

                if( a3->active)
                {
                    a3->fx += ux2;
                    a3->fy += uy2;
                    a3->fz += uz2;
                }

            }
        }
SKIP:
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
