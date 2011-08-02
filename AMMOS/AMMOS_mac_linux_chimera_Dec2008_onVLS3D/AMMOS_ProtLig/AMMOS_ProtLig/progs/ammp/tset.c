/* tset.c
*
* collection of routines to service  torsion searching
*
* POOP (Poor-mans Object Oriented Programming) using scope rules
*
*
*  tset forces a given torsion to be a given value
*  unlike tgroup tset is one deep and does not check energy.
*
*  tmin  searches one and only one torsion angle for an energy minimum
*
*
*  tmap i j k l ii jj kk ll nstep1 nstep2 ;
*   step through the two torsions and produce a map of the energies
*   the atoms are stored in the vx registers
*
*  tgroup/tsearch define a structure and use it
*  tset/tmin   do it on the fly.
*
*  tmin can be quite expensive
*
*/
/*
*  copyright 1993,1995 Robert W. Harrison
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
* but otherwise is self-contained. Note the hooks for Non-torsioned potentials
*/

#ifdef NEWCCALL
int tset( FILE *op, int echo ,int i1,int i2,int i3,int i4, float alpha )
#else
int tset(op,echo, i1,i2,i3,i4, alpha )
int i1,i2,i3,i4;
float alpha ;
FILE *op;
int echo;
#endif
{
    ATOM *ap,*ap1,*ap2,*ap3,*ap4, *a_m_serial(),*a_next();
    float get_torsion_value();
    float original,delta;
    int numatom,a_number();

    /* there must be atoms and they must be specified */
    numatom = a_number();
    if( numatom <= 0) return ;

    ap1 = a_m_serial(i1);
    if( ap1 == NULL){
        aaerror("An atom in tset is not defined");
        return ;}
    ap2 = a_m_serial(i2);
    if( ap2 == NULL){
        aaerror("An atom in tset is not defined");
        return ;}
    ap3 = a_m_serial(i3);
    if( ap3 == NULL){
        aaerror("An atom in tset is not defined");
        return ;}
    ap4 = a_m_serial(i4);
    if( ap4 == NULL){
        aaerror("An atom in tset is not defined");
        return ;}
    /* get the torsion value and check if we have to do the work */

    original = get_torsion_value( ap1,ap2,ap3,ap4);
    if( fabs(original  -alpha) < 1.e-3 ) return;

    /* set up the bond structure flags in ap->gx */
    tset_bond_build( ap1,ap2,ap3,ap4);

    delta = alpha -original;
    set_torsion( ap1,ap2,ap3,ap4,delta );
    if( echo )
        fprintf(op,"tset> original  %f delta %f final %f\n", original*180/3.14159,delta*180/3.14159,
                get_torsion_value( ap1,ap2,ap3,ap4)*180/3.14159);

}/* end of tset */
#ifdef NEWCCALL
int tmin( FILE *op, int echo ,int i1,int i2,int i3,int i4,int nstep, int (*vfs[])(),int nfs )
#else
int tmin(op,echo, i1,i2,i3,i4,nstep, vfs,nfs )
int i1,i2,i3,i4,nstep;
int (*vfs[])();
FILE *op;
int echo,nfs;
#endif
{
    ATOM *ap,*ap1,*ap2,*ap3,*ap4, *a_m_serial(),*a_next();
    float get_torsion_value();
    float original,delta;
    float vtemp,vmax;
    int numatom,a_number();
    int i,imax,ifs;

    /* there must be atoms and they must be specified */
    numatom = a_number();
    if( numatom <= 0) return ;

    ap1 = a_m_serial(i1);
    if( ap1 == NULL){
        aaerror("An atom in tset is not defined");
        return ;}
    ap2 = a_m_serial(i2);
    if( ap2 == NULL){
        aaerror("An atom in tset is not defined");
        return ;}
    ap3 = a_m_serial(i3);
    if( ap3 == NULL){
        aaerror("An atom in tset is not defined");
        return ;}
    ap4 = a_m_serial(i4);
    if( ap4 == NULL){
        aaerror("An atom in tset is not defined");
        return ;}

    if( nstep <= 0 ) nstep = 12;
    /* get the torsion value and check if we have to do the work */

    original = get_torsion_value( ap1,ap2,ap3,ap4);

    /* set up the bond structure flags in ap->gx */
    tset_bond_build( ap1,ap2,ap3,ap4);

    delta =  -original;
    set_torsion( ap1,ap2,ap3,ap4,delta );

    imax = -1;
    vmax = 10e20;
    delta = 2*3.141592653589793 /(float)nstep;

    for( i=0; i< nstep; i++)
    {
        vtemp = 0.;
        for( ifs = 0; ifs < nfs; ifs++)
            (*vfs[ifs])( &vtemp,0.);
        if( vtemp < vmax)
        {vmax = vtemp; imax = i;}
        set_torsion( ap1,ap2,ap3,ap4,delta );
    }
    set_torsion( ap1,ap2,ap3,ap4,imax*delta );


    if( echo )
        fprintf(op,"tset> original  %f  final %f\n", original*180/3.14159,
                get_torsion_value( ap1,ap2,ap3,ap4)*180/3.14159);

}/* end of tset */




int tset_bond_build( ap1,ap2,ap3,ap4)
ATOM *ap1,*ap2,*ap3,*ap4;
{
    ATOM *ap,*a_next(),*bonded[20];
    int i , numatom, a_number(),inbond;
    int j,tobonded;
    void get_bond();
    /* initialize the data */
    /* ap->gx < 0  not in bond structure
    *   == 0 in bond structure but check for neighbors 
    *  in bond structure, but neighbors are checked */
    numatom = a_number();
    for( i=0; i< numatom; i++)
    {  ap = a_next(i);
        ap->gx = -1.;
    }
    tobonded = 0;
    get_bond( ap3,bonded,20,&inbond);
    for( i=0; i< inbond; i++)
    {
        if( bonded[i] != ap1 && bonded[i]!= ap2)
        { bonded[i]->gx = 0.; tobonded++ ;}
    }
    get_bond( ap4,bonded,20,&inbond);
    for( i=0; i< inbond; i++)
    {
        if( bonded[i] != ap1 && bonded[i]!= ap2)
        { bonded[i]->gx = 0.; tobonded++ ;}
    }
    ap3->gx = 1.;
    ap4->gx = 1.;
    /* now for the meat */
    while( tobonded > 0 )
    {
        tobonded = 0;
        for( i=0; i< numatom;i++)
        {
            ap = a_next(i);
            if( ap != ap1 && ap != ap2) {
                if( ap->gx == 0.)
                {
                    ap->gx = 1.;
                    get_bond( ap,bonded,20,&inbond);
                    for( j=0; j< inbond; j++)
                    {
                        if( bonded[j]->gx < 1.)
                        { bonded[j]->gx = 0.; tobonded++ ;}
                    }
                } /*end of if ap->gx == 0 */
            }/* end of not ap1 not ap2 */
        } /* end for( i */
    }/*end while*/


}/* end of routine tset_bond_build */


float get_torsion_value( a1,a2,a3,a4)
ATOM *a1,*a2,*a3,*a4;
{
    float x1,x2,x3,y1,y2,y3,z1,z2,z3;
    float cx1,cy1,cz1,cx2,cy2,cz2;
    float dp,r;

    x1 = (a1->x -a2->x );
    y1 = (a1->y -a2->y );
    z1 = (a1->z -a2->z );
    x2 = (a3->x -a2->x );
    y2 = (a3->y -a2->y );
    z2 = (a3->z -a2->z );
    x3 = (a4->x -a3->x );
    y3 = (a4->y -a3->y );
    z3 = (a4->z -a3->z );
    /* 1 cross 2 */
    cx1 = y1*z2 - y2*z1;
    cy1 = -x1*z2 + x2*z1;
    cz1 = x1*y2 - x2*y1;
    r = cx1*cx1 + cy1*cy1 + cz1*cz1;
    if( r < 1.e-4) goto SKIP;
    r = 1./sqrt(r);
    cx1 = cx1*r;
    cy1 = cy1*r;
    cz1 = cz1*r;
    /* 3 cross 2 */
    cx2 = y3*z2 - y2*z3;
    cy2 = -x3*z2 + x2*z3;
    cz2 = x3*y2 - x2*y3;
    r = cx2*cx2 + cy2*cy2 + cz2*cz2;
    if( r < 1.e-4) goto SKIP;
    r = 1./sqrt(r);
    cx2 = cx2*r;
    cy2 = cy2*r;
    cz2 = cz2*r;
    /* if here everything is well determined */
    dp = cx1*cx2 + cy1*cy2 + cz1*cz2; /* cos( abs(theta)) */
    if( dp > 1.) dp = 1.; if( dp < -1.) dp = -1.;

    dp = acos(dp);
    /* determine the sign by triple product */
    r = cx1*x3 + cy1*y3 + cz1*z3;
    if( r > 0 ) dp =  -dp ;
    return dp;
SKIP:
    return 0;
} /*end of get_torsion_value() */


int set_torsion(ap1,ap2,ap3,ap4,howmuch)
ATOM *ap1,*ap2,*ap3,*ap4;
float howmuch;
{

    float nx,ny,nz;
    float phi,cphi,sphi,xphi;
    float rx,ry,rz, nnrx,nnry,nnrz, rnx,rny,rnz;
    ATOM *ap,*b1,*b2,*a_next();
    int numatom,a_number();
    int i;

    numatom = a_number();
    b1 = ap2; b2 = ap3;
    nx = b2->x - b1->x;
    ny = b2->y - b1->y;
    nz = b2->z - b1->z;
    rx = sqrt(nx*nx + ny*ny + nz*nz);
    if( rx < 1.e-6)
    {aaerror(" bad torsion radius in set_torsion \n"); return ;}
    nx = nx/rx;
    ny = ny/rx;
    nz = nz/rx;
    /*
            phi = 2.*3.141592653589793 /(float)tgp->ntry * (float)num;
    */
    cphi = cos(howmuch); sphi = sin(howmuch);
    phi = (1.-cphi);
    for( i=0; i< numatom; i++)
    {
        ap = a_next(i);
        if( ap->gx > 0. && ap != b2 )
        {
            rx = ap->x - b1->x;
            ry = ap->y - b1->y;
            rz = ap->z - b1->z;
            xphi = nx*rx + ny*ry + nz*rz;
            nnrx = xphi*nx;
            nnry = xphi*ny;
            nnrz = xphi*nz;
            rnx = ny*rz - nz*ry;
            rny = -nx*rz + nz*rx;
            rnz = nx*ry - ny*rx;
            rx = cphi*rx + phi*nnrx + sphi*rnx;
            ry = cphi*ry + phi*nnry + sphi*rny;
            rz = cphi*rz + phi*nnrz + sphi*rnz;
            ap->x = rx + b1->x;
            ap->y = ry + b1->y;
            ap->z = rz + b1->z;
        }
    }
}/* end of set_torsion */
#ifdef NEWCCALL
int tmap( FILE *op, int echo, int (*vfs[])(),int nfs,int i1,int i2,int i3,int i4,
          int j1 , int j2 ,int j3, int j4 , int nistep, int njstep)
#else
int tmap(op,echo,vfs,nfs,i1,i2,i3,i4,j1,j2,j3,j4,nistep,njstep )
int i1,i2,i3,i4;
int j1,j2,j3,j4;
int nistep,njstep;
int (*vfs[])(),nfs;
FILE *op;
int echo;
#endif
{
    ATOM *ap, *a_next();
    int numatm,a_number();
    int i,j,ifs;
    float V;
    float x,dx,y,dy;

    numatm = a_number();
    if( numatm <5 ) return; /* cant work if there aint enough atoms */

    /* save the atoms */
    for( i=0; i< numatm; i++)
    {
        ap = a_next(i);
        ap->vx = ap->x;
        ap->vy = ap->y;
        ap->vz = ap->z;
    }

    x = 0.;
    y = 0.;
    /* -180 to 180 convention is more popular */
    x = -3.141592653589793;
    dx = 0.;
    dy = 0.;
    if( nistep > 1)
        dx = 3.141592653589793*2./(nistep);
    if( njstep > 1)
        dy = 3.141592653589793*2./(njstep);
    fprintf(op,"%d %d Torsion map  %f %f steps\n",nistep,njstep,
            dx*180./3.141592653589793,dy*180./3.141592653589793);
    for( i=0; i< nistep; i++)
    {
        tset( op,0, i1,i2,i3,i4,x);
        y = 0.;
        y = -3.141592653589793;
        for( j=0; j< njstep; j++)
        {
            tset( op,0, j1,j2,j3,j4,y);
            V = 0.;
            for( ifs= 0; ifs< nfs; ifs++)
                (*vfs[ifs])(&V,0.);
            fprintf(op,"%f ",V);
            y = y + dy;
        }
        fprintf(op,"\n");
        x = x + dx;
    }
    /* restore the atoms */
    for( i=0; i< numatm; i++)
    {
        ap = a_next(i);
        ap->x = ap->vx;
        ap->y = ap->vy;
        ap->z = ap->vz;
    }

}/* end of tmap */
