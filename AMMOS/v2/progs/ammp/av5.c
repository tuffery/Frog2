/* av5.c
*
* collection of routines to service bond av5  potentials
*
*
*  av5 = b.b
*  where b = sum ( a2-a1,a3-a1, a4-a1,a5-a1) vectors
*
*  this term forces tetrahedral atoms to be tetrahedral
*  and introduces some mixing between angles and bonds
*  
*
* POOP (Poor-mans Object Oriented Programming) using scope rules
*
* these routines hold a data base (in terms of array indeces)
* of av5, with the associated length and force constant
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
* but otherwise is self-contained. Note the hooks for Non-av5ed potentials
*/
/*
* order is origin, right left apex (center), force constant and target value 
*/
typedef struct{
    ATOM  *atom1, *atom2, *atom3, *atom4 , *atom5;
    float k,offset;
    void  *next;
}  AV5;
#define AV5LONG sizeof(AV5)

AV5  *av5_first = NULL;
AV5  *av5_last = NULL;
/* function av5 adds a av5 to the av5 list
* returns 1 if ok
* returns 0 if not
*  is passed the array pointers, length and constant
* allocates the new memory, initializes it and
* returns
*/
int av5( p1,p2,p3,p4,p5,fk,off)
int p1,p2,p3,p4,p5;
float fk,off;
{
    AV5  *new;
    ATOM  *ap1, *ap2, *ap3, *ap4,*ap5, *a_m_serial();
    char line[80];
    /* get the atom pointers for the two serial numbers */
    ap1 = a_m_serial( p1 );
    ap2 = a_m_serial( p2 );
    ap3 = a_m_serial( p3 );
    ap4 = a_m_serial( p4 );
    ap5 = a_m_serial( p5 );
    if( (ap1 == NULL) || (ap2 == NULL) || (ap3==NULL) || (ap4==NULL) ||
            (ap5 == NULL) )
    {
        sprintf( line,"undefined atom in av5 %d %d %d %d %d\0",p1,p2,p3,p4,p5);
        aaerror( line );
        return 0;
    }

    if( ( new = malloc((unsigned int) AV5LONG ) ) == NULL)
    {
        return 0;
    }
    /* initialize the pointers */
    if( av5_first == NULL) av5_first = new;
    if( av5_last == NULL) av5_last = new;
    new -> atom1 = ap1;
    new -> atom2 = ap2;
    new -> atom3 = ap3;
    new -> atom4 = ap4;
    new -> atom5 = ap5;
    new -> offset = off;
    new -> k = fk;
    new -> next = new;
    av5_last -> next = new;
    av5_last = new;
    return 1;
}


/* v_av5()
* this function sums up the potentials
* for the atoms defined in the av5 data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int v_av5( V, lambda )
float *V,lambda;
{
    AV5  *bp;
    /* difference vectors */
    float x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;
    /* cross products and storage for normalizing */
    float r,cx,cy,cz;
    float hite;
    ATOM  *a1, *a2, *a3, *a4, *a5;


    bp = av5_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1; a2 = bp->atom2; a3 = bp->atom3;
        a4 = bp->atom4; a5 = bp->atom5;
        if( a1->active || a2->active || a3->active || a4->active
                || a5->active){
            x1 = (a2->x -a1->x +lambda*(a2->dx-a1->dx));
            y1 = (a2->y -a1->y +lambda*(a2->dy-a1->dy));
            z1 = (a2->z -a1->z +lambda*(a2->dz-a1->dz));
            x2 = (a3->x -a1->x +lambda*(a3->dx-a1->dx));
            y2 = (a3->y -a1->y +lambda*(a3->dy-a1->dy));
            z2 = (a3->z -a1->z +lambda*(a3->dz-a1->dz));
            x3 = (a4->x -a1->x +lambda*(a4->dx-a1->dx));
            y3 = (a4->y -a1->y +lambda*(a4->dy-a1->dy));
            z3 = (a4->z -a1->z +lambda*(a4->dz-a1->dz));
            x4 = (a5->x -a1->x +lambda*(a5->dx-a1->dx));
            y4 = (a5->y -a1->y +lambda*(a5->dy-a1->dy));
            z4 = (a5->z -a1->z +lambda*(a5->dz-a1->dz));

            cx = x1 + x2 + x3 + x4;
            cy = y1 + y2 + y3 + y4;
            cz = z1 + z2 + z3 + z4;

            hite = sqrt(cx*cx + cy*cy + cz*cz) - bp->offset;

            *V += bp->k*hite*hite;

        } /* if active */
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}

/* f_av5()
*
* f_av5 increments the forces in the atom structures by the force
* due to the av5 components.  NOTE THE WORD increment.
* the forces should first be zero'd.
* if not then this code will be invalid.  THIS IS DELIBERATE.
* on bigger (and better?) machines the different potential terms
* may be updated at random or in parrellel, if we assume that this routine
* will initialize the forces then we can't do this.
*/
int f_av5(lambda)
float lambda;
/*  returns 0 if error, 1 if OK */
{
    AV5  *bp;
    /* difference vectors */
    float x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;
    /* cross products and storage for normalizing */
    float r,cx1,cy1,cz1;
    float dx,dy,dz;
    float hite;
    float df;
    ATOM  *a1, *a2, *a3,  *a4,*a5 ;
    int i;


    bp = av5_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1; a2 = bp->atom2; a3 = bp->atom3;
        a4 = bp->atom4; a5 = bp->atom5;
        if( a1->active || a2->active || a3->active || a4->active
                || a5->active){
            x1 = (a2->x -a1->x +lambda*(a2->dx-a1->dx));
            y1 = (a2->y -a1->y +lambda*(a2->dy-a1->dy));
            z1 = (a2->z -a1->z +lambda*(a2->dz-a1->dz));
            x2 = (a3->x -a1->x +lambda*(a3->dx-a1->dx));
            y2 = (a3->y -a1->y +lambda*(a3->dy-a1->dy));
            z2 = (a3->z -a1->z +lambda*(a3->dz-a1->dz));
            x3 = (a4->x -a1->x +lambda*(a4->dx-a1->dx));
            y3 = (a4->y -a1->y +lambda*(a4->dy-a1->dy));
            z3 = (a4->z -a1->z +lambda*(a4->dz-a1->dz));
            x4 = (a5->x -a1->x +lambda*(a5->dx-a1->dx));
            y4 = (a5->y -a1->y +lambda*(a5->dy-a1->dy));
            z4 = (a5->z -a1->z +lambda*(a5->dz-a1->dz));

            cx1 = x1 + x2 + x3 + x4;
            cy1 = y1 + y2 + y3 + y4;
            cz1 = z1 + z2 + z3 + z4;

            hite = sqrt(cx1*cx1 + cy1*cy1 + cz1*cz1);

            df =  two*bp->k*(bp->offset - hite );
            /* do the  derivatives now  */
            r = two*sqrt(x1*x1 + y1*y1 + z1*z1);
            if( r > 1.e-6){
                a2->fx += df/r*(cx1+x1);
                a2->fy += df/r*(cy1+y1);
                a2->fz += df/r*(cz1+z1);
                a1->fx -= df/r*(cx1+x1);
                a1->fy -= df/r*(cy1+y1);
                a1->fz -= df/r*(cz1+z1);
            }
            r = two*sqrt(x2*x2 + y2*y2 + z2*z2);
            if( r > 1.e-6){
                a3->fx += df/r*(cx1+x2);
                a3->fy += df/r*(cy1+y2);
                a3->fz += df/r*(cz1+z2);
                a1->fx -= df/r*(cx1+x2);
                a1->fy -= df/r*(cy1+y2);
                a1->fz -= df/r*(cz1+z2);
            }
            r = two*sqrt(x3*x3 + y3*y3 + z3*z3);
            if( r > 1.e-6){
                a4->fx += df/r*(cx1+x3);
                a4->fy += df/r*(cy1+y3);
                a4->fz += df/r*(cz1+z3);
                a1->fx -= df/r*(cx1+x3);
                a1->fy -= df/r*(cy1+y3);
                a1->fz -= df/r*(cz1+z3);
            }
            r = two*sqrt(x4*x4 + y4*y4 + z4*z4);
            if( r > 1.e-6){
                a5->fx += df/r*(cx1+x4);
                a5->fy += df/r*(cy1+y4);
                a5->fz += df/r*(cz1+z4);
                a1->fx -= df/r*(cx1+x4);
                a1->fy -= df/r*(cy1+y4);
                a1->fz -= df/r*(cz1+z4);
            }
            if( a1->active == 0){ a1->fx = 0; a1->fy = 0.; a1->fz = 0;}
            if( a2->active == 0){ a2->fx = 0; a2->fy = 0.; a2->fz = 0;}
            if( a3->active == 0){ a3->fx = 0; a3->fy = 0.; a3->fz = 0;}
            if( a4->active == 0){ a4->fx = 0; a4->fy = 0.; a4->fz = 0;}
            if( a5->active == 0){ a5->fx = 0; a5->fy = 0.; a5->fz = 0;}

        } /* if active */

        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* function get_av5( a1,bonded,10,inbond);
* check the av5 list for atoms 1-4 ed to a1
*/
void get_av5( a1,bonded,mbond,inbond)
ATOM  *a1,  *bonded[];
int mbond,*inbond ;
{
    AV5  *mine;
    mine = av5_first;
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
/* routine dump_av5s
* this function outputs the av5 parameters
* and does it in a simple form
* av5 ser1,ser2,ser3,k,theta (in degrees )
* the rest is just free format
*/
void dump_av5s( where )
FILE *where;
{
    AV5  *b;
    ATOM  *a1, *a2, *a3, *a4, *a5;
    b = av5_first;
    if( b == NULL ) return;
    while( (b->next != b)  )
    {
        if( b->next == NULL) return;
        a1 = b->atom1; a2 = b->atom2;a3 = b->atom3; a4 = b->atom4;
        a5 = b->atom5;
        fprintf( where,"av5 %d %d %d %d %d %f %f \;\n",
                 a1->serial,a2->serial,
                 a3-> serial,a4->serial,a5->serial,b->k,b->offset);
        b = b->next;
    }
    if( b->next == NULL) return;
    a1 = b->atom1; a2 = b->atom2;a3 = b->atom3; a4 = b->atom4;
    a5 = b->atom5;
    fprintf( where,"av5 %d %d %d %d %d %f %f \;\n",
             a1->serial,a2->serial,
             a3-> serial,a4->serial,a5->serial,b->k,b->offset);
}

/* a_av5()
* this function sums up the potentials
* for the atoms defined in the av5 data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int a_av5( V, lambda ,ilow,ihigh,op)
float *V,lambda;
int ilow,ihigh;
FILE *op;
{
    AV5  *bp;
    /* difference vectors */
    float x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;
    /* cross products and storage for normalizing */
    float r,cx,cy,cz;
    float hite;
    ATOM  *a1, *a2, *a3, *a4, *a5;


    bp = av5_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1; a2 = bp->atom2; a3 = bp->atom3;
        a4 = bp->atom4; a5 = bp->atom5;
        if( (a1->serial >= ilow && a1->serial <= ihigh)
                ||  (a2->serial >= ilow && a2->serial <= ihigh)
                ||  (a3->serial >= ilow && a3->serial <= ihigh)
                ||  (a4->serial >= ilow && a4->serial <= ihigh)
                ||  (a5->serial >= ilow && a5->serial <= ihigh) )
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
            x4 = (a5->x -a1->x +lambda*(a5->dx-a1->dx));
            y4 = (a5->y -a1->y +lambda*(a5->dy-a1->dy));
            z4 = (a5->z -a1->z +lambda*(a5->dz-a1->dz));
            cx = x1 + x2 + x3 + x4;
            cy = y1 + y2 + y3 + y4;
            cz = z1 + z2 + z3 + z4;

            r = cx*cx + cy*cy + cz*cz;
            hite = sqrt(r);
            z2 = hite - bp->offset;


            z2 =  bp->k *z2*z2;
            *V += z2;
            fprintf(op,"AV5 %s %d %s %d %s %d %s %d %s %d E %f value %f error %f\n",
                    a1->name,a1->serial,a2->name,a2->serial,a3->name,a3->serial,a4->name,
                    a4->serial,a5->name,a5->serial,z2,hite,hite- bp->offset);
        }
SKIP:
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* ho_v_av5()
* this function sums up the potentials
* for the atoms defined in the av5 data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int v_ho_av5( V, lambda )
float *V,lambda;
{
    AV5  *bp;
    /* difference vectors */
    float x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;
    /* cross products and storage for normalizing */
    float r,cx,cy,cz;
    float hite;
    ATOM  *a1, *a2, *a3, *a4, *a5;
    float hol,get_f_variable();

    hol = get_f_variable("lambda");
    if( hol < 0.) hol = 0.;
    if( hol > 1.) hol = 1.;


    bp = av5_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1; a2 = bp->atom2; a3 = bp->atom3;
        a4 = bp->atom4; a5 = bp->atom5;
        if( a1->active || a2->active || a3->active || a4->active
                || a5->active){
            x1 = (a2->x -a1->x +lambda*(a2->dx-a1->dx));
            y1 = (a2->y -a1->y +lambda*(a2->dy-a1->dy));
            z1 = (a2->z -a1->z +lambda*(a2->dz-a1->dz));
            x2 = (a3->x -a1->x +lambda*(a3->dx-a1->dx));
            y2 = (a3->y -a1->y +lambda*(a3->dy-a1->dy));
            z2 = (a3->z -a1->z +lambda*(a3->dz-a1->dz));
            x3 = (a4->x -a1->x +lambda*(a4->dx-a1->dx));
            y3 = (a4->y -a1->y +lambda*(a4->dy-a1->dy));
            z3 = (a4->z -a1->z +lambda*(a4->dz-a1->dz));
            x4 = (a5->x -a1->x +lambda*(a5->dx-a1->dx));
            y4 = (a5->y -a1->y +lambda*(a5->dy-a1->dy));
            z4 = (a5->z -a1->z +lambda*(a5->dz-a1->dz));

            cx = x1 + x2 + x3 + x4;
            cy = y1 + y2 + y3 + y4;
            cz = z1 + z2 + z3 + z4;

            hite = (1.-hol)*(sqrt(cx*cx + cy*cy + cz*cz) - bp->offset);

            *V += bp->k*hite*hite;

        } /* if active */
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* f_av5()
*
* f_av5 increments the forces in the atom structures by the force
* due to the av5 components.  NOTE THE WORD increment.
* the forces should first be zero'd.
* if not then this code will be invalid.  THIS IS DELIBERATE.
* on bigger (and better?) machines the different potential terms
* may be updated at random or in parrellel, if we assume that this routine
* will initialize the forces then we can't do this.
*/
int f_ho_av5(lambda)
float lambda;
/*  returns 0 if error, 1 if OK */
{
    AV5  *bp;
    /* difference vectors */
    float x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;
    /* cross products and storage for normalizing */
    float r,cx1,cy1,cz1;
    float dx,dy,dz;
    float hite;
    float df;
    ATOM  *a1, *a2, *a3,  *a4,*a5 ;
    int i;
    float hol,get_f_variable();

    hol = get_f_variable("lambda");
    if( hol < 0.) hol = 0.;
    if( hol > 1.) hol = 1.;


    bp = av5_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1; a2 = bp->atom2; a3 = bp->atom3;
        a4 = bp->atom4; a5 = bp->atom5;
        if( a1->active || a2->active || a3->active || a4->active
                || a5->active){
            x1 = (a2->x -a1->x +lambda*(a2->dx-a1->dx));
            y1 = (a2->y -a1->y +lambda*(a2->dy-a1->dy));
            z1 = (a2->z -a1->z +lambda*(a2->dz-a1->dz));
            x2 = (a3->x -a1->x +lambda*(a3->dx-a1->dx));
            y2 = (a3->y -a1->y +lambda*(a3->dy-a1->dy));
            z2 = (a3->z -a1->z +lambda*(a3->dz-a1->dz));
            x3 = (a4->x -a1->x +lambda*(a4->dx-a1->dx));
            y3 = (a4->y -a1->y +lambda*(a4->dy-a1->dy));
            z3 = (a4->z -a1->z +lambda*(a4->dz-a1->dz));
            x4 = (a5->x -a1->x +lambda*(a5->dx-a1->dx));
            y4 = (a5->y -a1->y +lambda*(a5->dy-a1->dy));
            z4 = (a5->z -a1->z +lambda*(a5->dz-a1->dz));

            cx1 = x1 + x2 + x3 + x4;
            cy1 = y1 + y2 + y3 + y4;
            cz1 = z1 + z2 + z3 + z4;

            hite = sqrt(cx1*cx1 + cy1*cy1 + cz1*cz1);

            df =  two*bp->k*(1.-hol)*(1.-hol)*(bp->offset - hite );
            /* do the  derivatives now  */
            r = two*sqrt(x1*x1 + y1*y1 + z1*z1);
            if( r > 1.e-6){
                a2->fx += df/r*(cx1+x1);
                a2->fy += df/r*(cy1+y1);
                a2->fz += df/r*(cz1+z1);
                a1->fx -= df/r*(cx1+x1);
                a1->fy -= df/r*(cy1+y1);
                a1->fz -= df/r*(cz1+z1);
            }
            r = two*sqrt(x2*x2 + y2*y2 + z2*z2);
            if( r > 1.e-6){
                a3->fx += df/r*(cx1+x2);
                a3->fy += df/r*(cy1+y2);
                a3->fz += df/r*(cz1+z2);
                a1->fx -= df/r*(cx1+x2);
                a1->fy -= df/r*(cy1+y2);
                a1->fz -= df/r*(cz1+z2);
            }
            r = two*sqrt(x3*x3 + y3*y3 + z3*z3);
            if( r > 1.e-6){
                a4->fx += df/r*(cx1+x3);
                a4->fy += df/r*(cy1+y3);
                a4->fz += df/r*(cz1+z3);
                a1->fx -= df/r*(cx1+x3);
                a1->fy -= df/r*(cy1+y3);
                a1->fz -= df/r*(cz1+z3);
            }
            r = two*sqrt(x4*x4 + y4*y4 + z4*z4);
            if( r > 1.e-6){
                a5->fx += df/r*(cx1+x4);
                a5->fy += df/r*(cy1+y4);
                a5->fz += df/r*(cz1+z4);
                a1->fx -= df/r*(cx1+x4);
                a1->fy -= df/r*(cy1+y4);
                a1->fz -= df/r*(cz1+z4);
            }
            if( a1->active == 0){ a1->fx = 0; a1->fy = 0.; a1->fz = 0;}
            if( a2->active == 0){ a2->fx = 0; a2->fy = 0.; a2->fz = 0;}
            if( a3->active == 0){ a3->fx = 0; a3->fy = 0.; a3->fz = 0;}
            if( a4->active == 0){ a4->fx = 0; a4->fy = 0.; a4->fz = 0;}
            if( a5->active == 0){ a5->fx = 0; a5->fy = 0.; a5->fz = 0;}

        } /* if active */

        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
