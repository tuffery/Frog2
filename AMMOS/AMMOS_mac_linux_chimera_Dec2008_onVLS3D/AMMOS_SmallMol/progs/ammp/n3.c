/* n3.c
*
* collection of routines to service stepped
*
* step
*
*  stepped bonds
*  step i j d1 d2 d3 k1 k2 k3;
*  dv/dx = k1(x-d1)  x<d1
* dv/dx =  0 d1<x<d2
* dv/dx = k2(x-d2) x> d2
* dv/dx = k3(x-d3) x> d3
*
* POOP (Poor-mans Object Oriented Programming) using scope rules
*
* these routines hold a data base (in terms of array indeces)
* of STEPts, with the associated length and force constant
* These are updateable - unlike bonds which are "permanent"
*
* (this could be table driven but what the hell memories cheap)
*
* the routines for potential value, force and (eventually) second
* derivatives are here also
*
* force and 2nd derivative routines assume zero'd arrays for output
* this allows for parralellization if needed (on a PC?)
*
* forces are symmetric 
*/
/*
*  copyright 1992,1994 Robert W. Harrison
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
typedef struct{
    ATOM *atom1,*atom2;
    float d1,d2,d3,k1,k2,k3;
    void *next;
}  STEP;
#define STEPLONG sizeof(STEP)

STEP *STEP_first = NULL;
STEP *STEP_last = NULL;
/* function STEP adds a STEP to the STEP list
* returns 1 if ok
* returns 0 if not
*  is passed the atom serial numbers, length and constant
* allocates the new memory, initializes it and
* returns
*/
int step( p1,p2,d1,d2,d3,k1,k2,k3)
int p1,p2;
float d1,d2,d3,k1,k2,k3 ;
{
    ATOM *ap1,*ap2,*a_m_serial();
    STEP *new;
    float x;
    /*	char line[BUFSIZ]; */
    /* get the atom pointers for the two serial numbers */
    ap1 = a_m_serial( p1 );
    ap2 = a_m_serial( p2 );
    /* check the distance order so that d1<d2<d3 */
    if( !(d1 <=d2 && d2<=d3))
    {
        /* hard coded shell sort (uugh) */
        if( d1 >d3){ x = d1; d1 = d3; d3 = x; }
        if( d2 >d3){ x = d2; d2 = d3; d3 = x; }
        if( d1 >d2){ x = d1; d1 = d2; d2 = x; }
        if( d1 >d3){ x = d1; d1 = d3; d3 = x; }
        if( d2 >d3){ x = d2; d2 = d3; d3 = x; }
        if( d1 >d2){ x = d1; d1 = d2; d2 = x; }
        if( d1 >d3){ x = d1; d1 = d3; d3 = x; }
        if( d2 >d3){ x = d2; d2 = d3; d3 = x; }
        if( d1 >d2){ x = d1; d1 = d2; d2 = x; }
    }
    if( (ap1 == NULL) || (ap2 == NULL) )
    {
        /*      sprintf( line,"undefined atom in STEP %d %d \0",p1,p2);
        	aaerror( line );
        */
        return 1;
    }
    /* check to see if a STEPt is already defined */
    new = STEP_first;
    if( new != NULL)
    {
        while(1)
        {
            if( new == NULL) break;
            if( (new->atom1 == ap1 && new->atom2 == ap2) ||
                    (new->atom1 == ap2 && new->atom2 == ap1) )
            {
                new->d1=d1;
                new->d2=d2;
                new->d3=d3;
                new->k1=k1;
                new->k2=k2;
                new->k3=k3;
                return 1;
            }
            if( new == new->next) break;
            new = new->next;
        }
    }
    if( ( new = malloc( STEPLONG ) ) == NULL)
    {
        return 0;
    }
    /* initialize the pointers */
    if( STEP_first == NULL) STEP_first = new;
    if( STEP_last == NULL) STEP_last = new;
    new -> atom1 = ap1;
    new -> atom2 = ap2;
    new->d1=d1;
    new->d2=d2;
    new->d3=d3;
    new->k1=k1;
    new->k2=k2;
    new->k3=k3;
    new -> next = new;
    STEP_last -> next = new;
    STEP_last = new;
    return 1;
}




/* v_step()
* this function sums up the potentials
* for the atoms defined in the STEP data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int v_step( V, lambda )
float *V,lambda;
{
    STEP *bp;
    float r,xt,yt,zt;
    ATOM *a1,*a2;


    bp = STEP_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1; a2 = bp->atom2;
        if( a1->active || a2->active ){
            if( lambda == 0.)
            {
                r = (a1->x - a2->x)*(a1->x - a2->x);
                r = r + (a1->y - a2->y)*(a1->y - a2->y);
                r = r + (a1->z - a2->z)*(a1->z - a2->z);
            } else
            {
                xt = (a1->x -a2->x +lambda*(a1->dx-a2->dx));
                yt = (a1->y -a2->y +lambda*(a1->dy-a2->dy));
                zt = (a1->z -a2->z +lambda*(a1->dz-a2->dz));
                r = xt*xt+yt*yt+zt*zt;
            }
            r = sqrt(r);
            if( r < bp->d1)
            {
                r = r - bp->d1;
                *V += bp->k1 * r*r;
            } else if( r > bp->d2 && r < bp->d3) {
                r = r - bp->d2;
                *V += bp->k2 * r*r;
            } else if( r >= bp->d3 ){
                r = r - bp->d3;
                *V += bp->k3 * r*r + bp->k2*bp->d2*bp->d2;
            }
        }
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* f_step()
*
* f_step increments the forces in the atom structures by the force
* due to the STEP components.  NOTE THE WORD increment.
* the forces should first be zero'd.
* if not then this code will be invalid.  THIS IS DELIBERATE.
* on bigger (and better?) machines the different potential terms
* may be updated at random or in parrellel, if we assume that this routine
* will initialize the forces then we can't do this.
*/
int f_step(lambda)
float lambda;
/*  returns 0 if error, 1 if OK */
{
    STEP *bp;
    float r,t,ux,uy,uz;
    ATOM *a1,*a2;


    bp = STEP_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1; a2 = bp->atom2;
        if( a1->active || a2->active){
            if( lambda == 0.)
            {
                ux = (a2->x - a1->x);
                uy = (a2->y - a1->y);
                uz = (a2->z - a1->z);
            }else{
                ux = (a2->x -a1->x +lambda*(a2->dx-a1->dx));
                uy = (a2->y -a1->y +lambda*(a2->dy-a1->dy));
                uz = (a2->z -a1->z +lambda*(a2->dz-a1->dz));
            }
            r = ux*ux + uy*uy + uz*uz ;
            /* watch for FP errors*/
            if( r <= 1.e-5)
            { r = 0; ux = 1.; uy = 0.; uz = 0.;}else{
                r = sqrt(r); t = 1/r; ux = ux*t; uy = uy*t; uz = uz*t;
            }
            if( r < bp->d1)
            {
                r = r - bp->d1;
                ux = 2*bp->k1 * r *ux;
                uy = 2*bp->k1 * r *uy;
                uz = 2*bp->k1 * r *uz;
            } else if( r > bp->d2 && r < bp->d3) {
                r = r - bp->d2;
                ux = 2*bp->k2 * r *ux;
                uy = 2*bp->k2 * r *uy;
                uz = 2*bp->k2 * r *uz;
            }else if(r >= bp->d3){
                r = r - bp->d3;
                ux = 2.*bp->k3 *r *ux;
                uy = 2.*bp->k3 *r *uy;
                uz = 2.*bp->k3 *r *uz;
            }else{
                ux = 0.; uy = 0.; uz = 0.;
            }
            if( a1->active){
                a1->fx += ux;
                a1->fy += uy;
                a1->fz += uz;
            }

            if( a2->active) {
                a2->fx -= ux;
                a2->fy -= uy;
                a2->fz -= uz;
            }
        }
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* function get_step( a1,STEPed,10,inSTEP);
* check the STEPS list for atoms STEPed to a1
*/
void get_step( a1,STEPed,mSTEP,inSTEP)
ATOM *a1, *STEPed[];
int mSTEP,*inSTEP ;
{
    STEP *mine;
    mine = STEP_first;
    *inSTEP = 0;
    while(1)
    {
        if( (mine == NULL) )
        {
            return;
        }
        if( mine->atom1 == a1)
        {
            STEPed[(*inSTEP)++] = mine->atom2;
        }
        if( mine->atom2 == a1)
        {
            STEPed[(*inSTEP)++] = mine->atom1;
        }
        if( mine == mine->next) return;
        mine = mine->next;
        if( *inSTEP == mSTEP ) return;
    }
}
/* function get_step( a1,STEPed,10,inSTEP);
* check the STEPS list for atoms STEPed to a1
*/
void get_step_and_length( a1,STEPed,r,mSTEP,inSTEP)
ATOM *a1, *STEPed[];
int mSTEP,*inSTEP ;
float r[];
{
    STEP *mine;
    mine = STEP_first;
    *inSTEP = 0;
    while(1)
    {
        if( (mine == NULL) )
        {
            return;
        }
        if( mine->atom1 == a1)
        {
            r[*inSTEP] = mine->d1;
            STEPed[(*inSTEP)++] = mine->atom2;
        }
        if( mine->atom2 == a1)
        {
            r[*inSTEP] = mine->d1;
            STEPed[(*inSTEP)++] = mine->atom1;
        }
        if( mine == mine->next) return;
        mine = mine->next;
        if( *inSTEP == mSTEP ) return;
    }
}
/* STEP_next()
*  like bond_next() but for STEP structures
*/
int STEP_next( int i, ATOM **n1, ATOM **n2  )
{
    static STEP *np = NULL ;
    *n1 = NULL ; *n2 = NULL;
    if( STEP_first == NULL ) return (1==0);
if( np == NULL || i <= 0 ) { np = STEP_first;}
    else{ np = np->next; }
    *n1 = np->atom1; *n2 = np->atom2;
    if( np->next != np)return (1==1);
    return ( 1== 0);
}
/* routine dump_steps
* this function outputs the STEP parameters
* and does it in a simple form
* STEP ser1,ser2,k,req
* the rest is just free format
*/
void dump_steps( where )
FILE *where;
{
    STEP *b;
    ATOM *a1,*a2;
    b = STEP_first;
    if( b == NULL ) return;
    while( (b->next != b) )
    {
        if( b->next == NULL) return;
        a1 = b->atom1; a2 = b->atom2;
        fprintf( where,"step %d %d %f %f %f %f %f %f\;\n",a1->serial,a2->serial,
                 b->d1,b->d2,b->d3,b->k1,b->k2,b->k3);
        b = b->next;
    }
    if( b->next == NULL) return;
    a1 = b->atom1; a2 = b->atom2;
    fprintf( where,"step %d %d %f %f %f %f %f %f\;\n",a1->serial,a2->serial,
             b->d1,b->d2,b->d3,b->k1,b->k2,b->k3);
}

/* a_step()
* this function sums up the potentials
* for the atoms defined in the STEP data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int a_step( V, lambda,ilow,ihigh,op )
float *V,lambda;
int ilow,ihigh;
FILE *op;
{
    STEP *bp;
    float r,xt,yt,zt;
    ATOM *a1,*a2;


    bp = STEP_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1; a2 = bp->atom2;
        if(( a1->serial >= ilow && a1->serial <=ihigh)
                ||( a2->serial >= ilow && a2->serial <=ihigh))
        {
            if( lambda == 0.)
            {
                r = (a1->x - a2->x)*(a1->x - a2->x);
                r = r + (a1->y - a2->y)*(a1->y - a2->y);
                r = r + (a1->z - a2->z)*(a1->z - a2->z);
            } else
            {
                xt = (a1->x -a2->x +lambda*(a1->dx-a2->dx));
                yt = (a1->y -a2->y +lambda*(a1->dy-a2->dy));
                zt = (a1->z -a2->z +lambda*(a1->dz-a2->dz));
                r = xt*xt+yt*yt+zt*zt;
            }
            r = sqrt(r);
            zt = 0;
            if( r < bp->d1)
                zt= bp->k1*( r - bp->d1)*(r - bp->d1);
            if( r > bp->d2 && r < bp->d3)
                zt= bp->k2*( r - bp->d2)*(r - bp->d2);
            if( r >= bp->d3)
                zt= bp->k3*( r - bp->d3)*(r - bp->d3)+bp->k2*bp->d2*bp->d2;
            *V += zt;
            fprintf(op,"step %s %d %s %d E %f value %f\n"
                    ,a1->name,a1->serial,a2->name,a2->serial,zt,r);
        }
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* gsdg_step( ATOM *ap )
*  
* setup the distances for STEP terms
*/
int gsdg_step( ap)
ATOM *ap;
{
    ATOM *bp;
    STEP *np;

    np = STEP_first;
    while(1)
    { if( np == NULL ) return 0;
        if( np->atom1 == ap )
        {  bp = np->atom2; bp->vx = (np->d1*np->d2 );
            bp->vy = np->k1; }
        if( np->atom2 == ap )
        {  bp = np->atom1; bp->vx = (np->d1*np->d2 );
            bp->vy = np->k1; }

        if( np == np->next ) return 0;
        np = np->next;
    }
}
/* get step and bounds
*  used mostly by kohonen to get information 
* about all the atoms related to me by a step potential
*/
void get_step_and_bounds( a1,steped,r,rmid,rmax,mstep,instep)
ATOM *a1, *steped[];
int mstep,*instep ;
float r[],rmid[],rmax[];
{
    STEP *mine;
    mine = STEP_first;
    *instep = 0;
    while(1)
    {
        if( (mine == NULL) )
        {
            return;
        }
        if( mine->atom1 == a1)
        {
            r[*instep] = mine->d1;
            rmid[*instep] = mine->d2;
            rmax[*instep] = mine->d3;
            steped[(*instep)++] = mine->atom2;
        }
        if( mine->atom2 == a1)
        {
            r[*instep] = mine->d1;
            rmid[*instep] = mine->d2;
            rmax[*instep] = mine->d3;
            steped[(*instep)++] = mine->atom1;
        }
        if( mine == mine->next) return;
        mine = mine->next;
        if( *instep == mstep ) return;
    }
}