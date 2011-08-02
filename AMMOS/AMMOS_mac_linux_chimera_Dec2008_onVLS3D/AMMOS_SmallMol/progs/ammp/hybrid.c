/* hybrid.c
*
* collection of routines to service bond hybrid  potentials
*
* these are the planar and pyrmid height potentials
*
*  
*
* POOP (Poor-mans Object Oriented Programming) using scope rules
*
* these routines hold a data base (in terms of array indeces)
* of hybrid, with the associated length and force constant
*
* (this could be table driven but what the hell memories cheap)
*
* the routines for potential value, force and (eventually) second
* derivatives are here also
*
* force and 2nd derivative routines assume zero'd arrays for output
* this allows for parralellization if needed (on a PC?)
*
* forces are bond wise symmetric - so we don't have to fuck around with
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
* but otherwise is self-contained. Note the hooks for Non-hybrided potentials
*/
/*
* order is origin, right left apex (center), force constant and target value 
*/
typedef struct{
    ATOM  *atom1, *atom2, *atom3, *atom4;
    float k,offset;
    void  *next;
}  HYBRID;
#define HLONG sizeof(HYBRID)

HYBRID  *hybrid_first = NULL;
HYBRID  *hybrid_last = NULL;
/* function hybrid adds a hybrid to the hybrid list
* returns 1 if ok
* returns 0 if not
*  is passed the array pointers, length and constant
* allocates the new memory, initializes it and
* returns
*/
int hybrid( p1,p2,p3,p4,fk,off)
int p1,p2,p3,p4;
float fk,off;
{
    HYBRID  *new;
    ATOM  *ap1, *ap2, *ap3, *ap4, *a_m_serial();
    char line[BUFSIZ];
    /* get the atom pointers for the two serial numbers */
    ap1 = a_m_serial( p1 );
    ap2 = a_m_serial( p2 );
    ap3 = a_m_serial( p3 );
    ap4 = a_m_serial( p4 );
    if( (ap1 == NULL) || (ap2 == NULL) || (ap3==NULL) || (ap4==NULL) )
    {
        sprintf( line,"undefined atom in hybrid %d %d %d %d \0",p1,p2,p3,p4);
        aaerror( line );
        return 0;
    }

    if( ( new = malloc((unsigned int) HLONG ) ) == NULL)
    {
        return 0;
    }
    /* initialize the pointers */
    if( hybrid_first == NULL) hybrid_first = new;
    if( hybrid_last == NULL) hybrid_last = new;
    new -> atom1 = ap1;
    new -> atom2 = ap2;
    new -> atom3 = ap3;
    new -> atom4 = ap4;
    new -> offset = off;
    new -> k = fk;
    new -> next = new;
    hybrid_last -> next = new;
    hybrid_last = new;
    return 1;
}


/* v_hybrid()
* this function sums up the potentials
* for the atoms defined in the hybrid data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int v_hybrid( V, lambda )
float *V,lambda;
{
    HYBRID  *bp;
    /* difference vectors */
    float x1,y1,z1,x2,y2,z2,x3,y3,z3;
    /* cross products and storage for normalizing */
    float r,cx1,cy1,cz1;
    float hite;
    ATOM  *a1, *a2, *a3, *a4;


    bp = hybrid_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1; a2 = bp->atom2; a3 = bp->atom3;
        a4 = bp->atom4;
        if( a1->active || a2->active || a3->active || a4->active){
            x1 = (a2->x -a1->x +lambda*(a2->dx-a1->dx));
            y1 = (a2->y -a1->y +lambda*(a2->dy-a1->dy));
            z1 = (a2->z -a1->z +lambda*(a2->dz-a1->dz));
            x2 = (a3->x -a1->x +lambda*(a3->dx-a1->dx));
            y2 = (a3->y -a1->y +lambda*(a3->dy-a1->dy));
            z2 = (a3->z -a1->z +lambda*(a3->dz-a1->dz));
            x3 = (a4->x -a1->x +lambda*(a4->dx-a1->dx));
            y3 = (a4->y -a1->y +lambda*(a4->dy-a1->dy));
            z3 = (a4->z -a1->z +lambda*(a4->dz-a1->dz));
            /* 1 cross 2 */
            cx1 = y1*z2 - y2*z1;
            cy1 = -x1*z2 + x2*z1;
            cz1 = x1*y2 - x2*y1;
            r = cx1*cx1 + cy1*cy1 + cz1*cz1;
            if( r < 1.e-16) goto SKIP;
            r = sqrt(r);
            /* height is x3 vectordot cx1( normalized) */
            hite = cx1*x3 + cy1*y3 + cz1*z3; hite = hite/r;
            *V += (bp->k)*(hite - bp->offset)*(hite - bp->offset) ;
        }
SKIP:
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}

/* f_hybrid()
*
* f_hybrid increments the forces in the atom structures by the force
* due to the hybrid components.  NOTE THE WORD increment.
* the forces should first be zero'd.
* if not then this code will be invalid.  THIS IS DELIBERATE.
* on bigger (and better?) machines the different potential terms
* may be updated at random or in parrellel, if we assume that this routine
* will initialize the forces then we can't do this.
*/
int f_hybrid(lambda)
float lambda;
/*  returns 0 if error, 1 if OK */
{
    HYBRID  *bp;
    /* difference vectors */
    float x1,y1,z1,x2,y2,z2,x3,y3,z3;
    /* cross products and storage for normalizing */
    float r,cx1,cy1,cz1;
    float dx,dy,dz;
    float hite;
    float df;
    float r3,c;
    ATOM  *a1, *a2, *a3,  *a4, *at;
    int i;


    bp = hybrid_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1; a2 = bp->atom2; a3 = bp->atom3;
        a4 = bp->atom4;
        if( a1->active || a2->active || a3->active || a4->active){
            for( i=0; i< 3; i++)
            {
                x1 = (a2->x -a1->x +lambda*(a2->dx-a1->dx));
                y1 = (a2->y -a1->y +lambda*(a2->dy-a1->dy));
                z1 = (a2->z -a1->z +lambda*(a2->dz-a1->dz));
                x2 = (a3->x -a1->x +lambda*(a3->dx-a1->dx));
                y2 = (a3->y -a1->y +lambda*(a3->dy-a1->dy));
                z2 = (a3->z -a1->z +lambda*(a3->dz-a1->dz));
                x3 = (a4->x -a1->x +lambda*(a4->dx-a1->dx));
                y3 = (a4->y -a1->y +lambda*(a4->dy-a1->dy));
                z3 = (a4->z -a1->z +lambda*(a4->dz-a1->dz));
                /* 1 cross 2 */
                cx1 = y1*z2 - y2*z1;
                cy1 = -x1*z2 + x2*z1;
                cz1 = x1*y2 - x2*y1;
                r = cx1*cx1 + cy1*cy1 + cz1*cz1;
                if( r < 1.e-16) goto SKIP;
                r = sqrt(r);r3 = r*r*r;
                /* height is x3 vectordot cx1( normalized) */
                hite = cx1*x3 + cy1*y3 + cz1*z3; hite = hite/r;
                df =  2*bp->k*(bp->offset - hite)/3;
                /* do the apex derivatives now (easy) */
                a4->fx += df/r*cx1;
                a4->fy += df/r*cy1;
                a4->fz += df/r*cz1;
                a1->fx -= df/r*cx1;
                a1->fy -= df/r*cy1;
                a1->fz -= df/r*cz1;
                /* now the side derivative (messy) */
                /* dot product of r3 and cx1 components */
                dx =  -cx1*x3/r3*df;
                dy =  -cy1*y3/r3*df;
                dz =  -cz1*z3/r3*df;
                c = df*( (-y3*z2+z3*y2)/r) ;
                c += dx*( y2*(x1*y2-x2*y1) - z2*(x2*z1-x1*z2));
                a2->fx += c; a1->fx -= c;
                c = df*( (-z3*x2+x3*z2)/r) ;
                c += dy*( z2*(y1*z2-y2*z1) - x2*(x1*y2-x2*y1));
                a2->fy += c; a1->fy -= c;
                c = df*( (-x3*y2+y3*x2)/r) ;
                c += dz*( x2*(x2*z1-x1*z2) - y2*(y1*z2-y2*z1));
                a2->fz += c; a1->fz -= c;
                c = df*( (-z3*y1+y3*z1)/r) ;
                c -= dx*( y1*(x1*y2-x2*y1) - z1*(x2*z1-x1*z2));
                a3->fx += c; a1->fx -= c;
                c = df*( (-x3*z1+z3*x1)/r) ;
                c -= dy*( z1*(y1*z2-y2*z1) - x1*(x1*y2-x2*y1));
                a3->fy += c; a1->fy -= c;
                c = df*( (-y3*x1+x3*y1)/r) ;
                c -= dz*( x1*(x2*z1-x1*z2) - y1*(y1*z2-y2*z1));
                a3->fz += c; a1->fz -= c;
                /* circularly shift the base atoms */
                at = a1; a1 = a2; a2 = a3; a3 = at;
            }
        if( a1->active == 0){ a1->fx = 0; a1->fy = 0.; a1->fz = 0;}
            if( a2->active == 0){ a2->fx = 0; a2->fy = 0.; a2->fz = 0;}
            if( a3->active == 0){ a3->fx = 0; a3->fy = 0.; a3->fz = 0;}
            if( a4->active == 0){ a4->fx = 0; a4->fy = 0.; a4->fz = 0;}
        }
SKIP:
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* function get_hybrid( a1,bonded,10,inbond);
* check the hybrid list for atoms 1-4 ed to a1
*/
void get_hybrid( a1,bonded,mbond,inbond)
ATOM  *a1,  *bonded[];
int mbond,*inbond ;
{
    HYBRID  *mine;
    mine = hybrid_first;
    *inbond = 0;
    while(1)
    {
        if( (mine == NULL) )
        {
            return;
        }
        if( mine->atom1 == a1)
        {
            bonded[(*inbond)++] = mine->atom4;
        }
        if( mine->atom4 == a1)
        {
            bonded[(*inbond)++] = mine->atom1;
        }
        if( mine == mine->next) return;
        mine = mine->next;
        if( *inbond == mbond ) return;
    }
}
/* routine dump_hybrids
* this function outputs the hybrid parameters
* and does it in a simple form
* hybrid ser1,ser2,ser3,k,theta (in degrees )
* the rest is just free format
*/
void dump_hybrids( where )
FILE *where;
{
    HYBRID  *b;
    ATOM  *a1, *a2, *a3, *a4;
    b = hybrid_first;
    if( b == NULL ) return;
    while( (b->next != b)  )
    {
        if( b->next == NULL) return;
        a1 = b->atom1; a2 = b->atom2;a3 = b->atom3; a4 = b->atom4;
        fprintf( where,"hybrid %d %d %d %d %f %f \;\n",
                 a1->serial,a2->serial,
                 a3-> serial,a4->serial,b->k,b->offset);
        b = b->next;
    }
    if( b->next == NULL) return;
    a1 = b->atom1; a2 = b->atom2;a3 = b->atom3; a4 = b->atom4;
    fprintf( where,"hybrid %d %d %d %d %f %f \;\n",
             a1->serial,a2->serial,
             a3-> serial,a4->serial,b->k,b->offset);
}

/* a_hybrid()
* this function sums up the potentials
* for the atoms defined in the hybrid data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int a_hybrid( V, lambda ,ilow,ihigh,op)
float *V,lambda;
int ilow,ihigh;
FILE *op;
{
    HYBRID  *bp;
    /* difference vectors */
    float x1,y1,z1,x2,y2,z2,x3,y3,z3;
    /* cross products and storage for normalizing */
    float r,cx1,cy1,cz1;
    float hite;
    ATOM  *a1, *a2, *a3, *a4;


    bp = hybrid_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1; a2 = bp->atom2; a3 = bp->atom3;
        a4 = bp->atom4;
        if( (a1->serial >= ilow && a1->serial <= ihigh)
                ||  (a2->serial >= ilow && a2->serial <= ihigh)
                ||  (a3->serial >= ilow && a3->serial <= ihigh)
                ||  (a4->serial >= ilow && a4->serial <= ihigh) )
        {

            x1 = (a2->x -a1->x +lambda*(a2->dx-a1->dx));
            y1 = (a2->y -a1->y +lambda*(a2->dy-a1->dy));
            z1 = (a2->z -a1->z +lambda*(a2->dz-a1->dz));
            x2 = (a3->x -a1->x +lambda*(a3->dx-a1->dx));
            y2 = (a3->y -a1->y +lambda*(a3->dy-a1->dy));
            z2 = (a3->z -a1->z +lambda*(a3->dz-a1->dz));
            x3 = (a4->x -a1->x +lambda*(a4->dx-a1->dx));
            y3 = (a4->y -a1->y +lambda*(a4->dy-a1->dy));
            z3 = (a4->z -a1->z +lambda*(a4->dz-a1->dz));
            /* 1 cross 2 */
            cx1 = y1*z2 - y2*z1;
            cy1 = -x1*z2 + x2*z1;
            cz1 = x1*y2 - x2*y1;
            r = cx1*cx1 + cy1*cy1 + cz1*cz1;
            if( r < 1.e-16) goto SKIP;
            r = sqrt(r);
            /* height is x3 vectordot cx1( normalized) */
            hite = cx1*x3 + cy1*y3 + cz1*z3; hite = hite/r;
            z2 = (bp->k)*(hite - bp->offset)*(hite - bp->offset) ;
            *V += z2;
            fprintf(op,"Hybrid %s %d %s %d %s %d %s %d E %f value %f error %f\n",
                    a1->name,a1->serial,a2->name,a2->serial,a3->name,a3->serial,a4->name,
                    a4->serial,z2,hite,hite- bp->offset);
        }
SKIP:
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* v_ho_hybrid()
* this function sums up the potentials
* for the atoms defined in the hybrid data structure.
*
* homotopy version
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int v_ho_hybrid( V, lambda )
float *V,lambda;
{
    HYBRID  *bp;
    /* difference vectors */
    float x1,y1,z1,x2,y2,z2,x3,y3,z3;
    /* cross products and storage for normalizing */
    float r,cx1,cy1,cz1;
    float hite;
    ATOM  *a1, *a2, *a3, *a4;
    float get_f_variable();
    float hol;


    hol = get_f_variable("lambda");
    if( hol >= 1.) return;
    if( hol <= 0.) hol = 0.;

    bp = hybrid_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1; a2 = bp->atom2; a3 = bp->atom3;
        a4 = bp->atom4;
        if( a1->active || a2->active || a3->active || a4->active){
            x1 = (a2->x -a1->x +lambda*(a2->dx-a1->dx));
            y1 = (a2->y -a1->y +lambda*(a2->dy-a1->dy));
            z1 = (a2->z -a1->z +lambda*(a2->dz-a1->dz));
            x2 = (a3->x -a1->x +lambda*(a3->dx-a1->dx));
            y2 = (a3->y -a1->y +lambda*(a3->dy-a1->dy));
            z2 = (a3->z -a1->z +lambda*(a3->dz-a1->dz));
            x3 = (a4->x -a1->x +lambda*(a4->dx-a1->dx));
            y3 = (a4->y -a1->y +lambda*(a4->dy-a1->dy));
            z3 = (a4->z -a1->z +lambda*(a4->dz-a1->dz));
            /* 1 cross 2 */
            cx1 = y1*z2 - y2*z1;
            cy1 = -x1*z2 + x2*z1;
            cz1 = x1*y2 - x2*y1;
            r = cx1*cx1 + cy1*cy1 + cz1*cz1;
            if( r < 1.e-16) goto SKIP;
            r = sqrt(r);
            /* height is x3 vectordot cx1( normalized) */
            hite = cx1*x3 + cy1*y3 + cz1*z3; hite = hite/r;
            r =  hite*(one+hol) - bp->offset*(one-hol)  ;
            *V += (bp->k)*r*r ;
        }
SKIP:
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}

/* f_ho_hybrid()
*
* f_ho_hybrid increments the forces in the atom structures by the force
* due to the hybrid components.  NOTE THE WORD increment.
* the forces should first be zero'd.
* if not then this code will be invalid.  THIS IS DELIBERATE.
* on bigger (and better?) machines the different potential terms
* may be updated at random or in parrellel, if we assume that this routine
* will initialize the forces then we can't do this.
*
* homotopy method version
*  
*/
int f_ho_hybrid(lambda)
float lambda;
/*  returns 0 if error, 1 if OK */
{
    HYBRID  *bp;
    /* difference vectors */
    float x1,y1,z1,x2,y2,z2,x3,y3,z3;
    /* cross products and storage for normalizing */
    float r,cx1,cy1,cz1;
    float dx,dy,dz;
    float hite;
    float df;
    float r3,c;
    ATOM  *a1, *a2, *a3,  *a4, *at;
    int i;
    float get_f_variable(),hol;

    hol = get_f_variable("lambda");
    if( hol >= 1.) return;
    if( hol <= 0.) hol = 0.;

    bp = hybrid_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1; a2 = bp->atom2; a3 = bp->atom3;
        a4 = bp->atom4;
        if( a1->active || a2->active || a3->active || a4->active){
            for( i=0; i< 3; i++)
            {
                x1 = (a2->x -a1->x +lambda*(a2->dx-a1->dx));
                y1 = (a2->y -a1->y +lambda*(a2->dy-a1->dy));
                z1 = (a2->z -a1->z +lambda*(a2->dz-a1->dz));
                x2 = (a3->x -a1->x +lambda*(a3->dx-a1->dx));
                y2 = (a3->y -a1->y +lambda*(a3->dy-a1->dy));
                z2 = (a3->z -a1->z +lambda*(a3->dz-a1->dz));
                x3 = (a4->x -a1->x +lambda*(a4->dx-a1->dx));
                y3 = (a4->y -a1->y +lambda*(a4->dy-a1->dy));
                z3 = (a4->z -a1->z +lambda*(a4->dz-a1->dz));
                /* 1 cross 2 */
                cx1 = y1*z2 - y2*z1;
                cy1 = -x1*z2 + x2*z1;
                cz1 = x1*y2 - x2*y1;
                r = cx1*cx1 + cy1*cy1 + cz1*cz1;
                if( r < 1.e-16) goto SKIP;
                r = sqrt(r);r3 = r*r*r;
                /* height is x3 vectordot cx1( normalized) */
                hite = cx1*x3 + cy1*y3 + cz1*z3; hite = hite/r;
                /*
                	df =  2*bp->k*(bp->offset - hite)/3;
                */
                df =  2*bp->k*(one-hol)*((one-hol)*bp->offset - (one+hol)*hite)/3;
                /* do the apex derivatives now (easy) */
                a4->fx += df/r*cx1;
                a4->fy += df/r*cy1;
                a4->fz += df/r*cz1;
                a1->fx -= df/r*cx1;
                a1->fy -= df/r*cy1;
                a1->fz -= df/r*cz1;
                /* now the side derivative (messy) */
                /* dot product of r3 and cx1 components */
                dx =  -cx1*x3/r3*df;
                dy =  -cy1*y3/r3*df;
                dz =  -cz1*z3/r3*df;
                c = df*( (-y3*z2+z3*y2)/r) ;
                c += dx*( y2*(x1*y2-x2*y1) - z2*(x2*z1-x1*z2));
                a2->fx += c; a1->fx -= c;
                c = df*( (-z3*x2+x3*z2)/r) ;
                c += dy*( z2*(y1*z2-y2*z1) - x2*(x1*y2-x2*y1));
                a2->fy += c; a1->fy -= c;
                c = df*( (-x3*y2+y3*x2)/r) ;
                c += dz*( x2*(x2*z1-x1*z2) - y2*(y1*z2-y2*z1));
                a2->fz += c; a1->fz -= c;
                c = df*( (-z3*y1+y3*z1)/r) ;
                c -= dx*( y1*(x1*y2-x2*y1) - z1*(x2*z1-x1*z2));
                a3->fx += c; a1->fx -= c;
                c = df*( (-x3*z1+z3*x1)/r) ;
                c -= dy*( z1*(y1*z2-y2*z1) - x1*(x1*y2-x2*y1));
                a3->fy += c; a1->fy -= c;
                c = df*( (-y3*x1+x3*y1)/r) ;
                c -= dz*( x1*(x2*z1-x1*z2) - y1*(y1*z2-y2*z1));
                a3->fz += c; a1->fz -= c;
                /* circularly shift the base atoms */
                at = a1; a1 = a2; a2 = a3; a3 = at;
            }
        if( a1->active == 0){ a1->fx = 0; a1->fy = 0.; a1->fz = 0;}
            if( a2->active == 0){ a2->fx = 0; a2->fy = 0.; a2->fz = 0;}
            if( a3->active == 0){ a3->fx = 0; a3->fy = 0.; a3->fz = 0;}
            if( a4->active == 0){ a4->fx = 0; a4->fy = 0.; a4->fz = 0;}
        }
SKIP:
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* gsdg_hybrid( atom *who)
*  if *who is active make it have the right hybrid value
*  to be used before other regularizations in both
* gsdg and bell
*
*/
int gsdg_hybrid( who)
ATOM *who;
{
    ATOM *a1,*a2,*a3,*a4;
    int i;
    float x3,y3,z3;
    float x1,y1,z1;
    float x2,y2,z2;
    float hite,r,cx1,cy1,cz1;
    HYBRID *hp;
    if( hybrid_first == NULL ) return 0;
    hp = hybrid_first;
    if( !who->active ) return 0;
    while( hp != NULL)
    {

        a1 = hp->atom1; a2 = hp->atom2; a3 = hp->atom3;
        a4 = hp->atom4;
        if( a4 == who){
            x1 = (a2->x -a1->x );
            y1 = (a2->y -a1->y );
            z1 = (a2->z -a1->z );
            x2 = (a3->x -a1->x );
            y2 = (a3->y -a1->y );
            z2 = (a3->z -a1->z );
            x3 = (a4->x -a1->x );
            y3 = (a4->y -a1->y );
            z3 = (a4->z -a1->z );
            /* 1 cross 2 */
            cx1 = y1*z2 - y2*z1;
            cy1 = -x1*z2 + x2*z1;
            cz1 = x1*y2 - x2*y1;
            r = cx1*cx1 + cy1*cy1 + cz1*cz1;
            if( r < 1.e-16) goto SKIP;
            r = sqrt(r);r = 1./r;
            cx1 *= r;
            cy1 *= r;
            cz1 *= r;
            /* height is x3 vectordot cx1( normalized) */
            hite = cx1*x3 + cy1*y3 + cz1*z3;
            hite = hp->offset - hite;
            cx1 *= hite;
            cy1 *= hite;
            cz1 *= hite;
            a4->x += cx1;
            a4->y += cy1;
            a4->z += cz1;
SKIP: ;
        }/* a4 is target atom */
        if(hp->next == hp) break;
        hp = hp->next;
    }
}/* end of gsdg_hybrid */
