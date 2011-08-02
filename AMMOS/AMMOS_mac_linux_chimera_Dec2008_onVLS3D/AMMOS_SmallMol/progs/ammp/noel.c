/* noel.c
*
* collection of routines to service noel potentials
*
* noel  == NOE Length
*
*  stepped bonds
*  noel i j d d- d+ kl kh
*  v =  (r - ( d - d-) )**2 *kl r < d-
*  v =  (r - ( d + d+) )**2 *kh r > d+
*
*  kl,kh can be of any value 
*
* POOP (Poor-mans Object Oriented Programming) using scope rules
*
* these routines hold a data base (in terms of array indeces)
* of noelts, with the associated length and force constant
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
    float d,dm,dh,km,kh;
    int ismet;
    void *next;
}  NOEL;
#define NOELLONG sizeof(NOEL)

NOEL *noel_first = NULL;
NOEL *noel_last = NULL;
/* function noel adds a noel to the noel list
* returns 1 if ok
* returns 0 if not
*  is passed the atom serial numbers, length and constant
* allocates the new memory, initializes it and
* returns
*/
int noel( p1,p2,d,dm,dh,km,kh)
int p1,p2;
float d,dm,dh,km,kh ;
{
    ATOM *ap1,*ap2,*a_m_serial();
    NOEL *new;
    /*	char line[BUFSIZ]; */
    /* get the atom pointers for the two serial numbers */
    ap1 = a_m_serial( p1 );
    ap2 = a_m_serial( p2 );
    if( (ap1 == NULL) || (ap2 == NULL) )
    {
        /*      sprintf( line,"undefined atom in noel %d %d \0",p1,p2);
        	aaerror( line );
        */
        return 1;
    }
    /* check to see if a noelt is already defined */
    new = noel_first;
    if( new != NULL)
    {
        while(1)
        {
            if( new == NULL) break;
            if( (new->atom1 == ap1 && new->atom2 == ap2) ||
                    (new->atom1 == ap2 && new->atom2 == ap1) )
            {
                new->ismet = (1==1);
                new->d = d; new->dm = dm; new->dh = dh;
                new->km = km; new->kh = kh; return 1;
            }
            if( new == new->next) break;
            new = new->next;
        }
    }
    if( ( new = malloc( NOELLONG ) ) == NULL)
    {
        return 0;
    }
    /* initialize the pointers */
    if( noel_first == NULL) noel_first = new;
    if( noel_last == NULL) noel_last = new;
    new -> atom1 = ap1;
    new -> atom2 = ap2;
    new->ismet = (1==1);
    new->d = d; new->dm = dm; new->dh = dh;
    new->km = km; new->kh = kh;
    new -> next = new;
    noel_last -> next = new;
    noel_last = new;
    return 1;
}

int noel_generate( dm,dh)
float dm,dh;
{
    NOEL* np;
    ATOM *ap1,*ap2;
    float r,dx,dy,dz;

    if( noel_first == NULL) return 1;
    np = noel_first;
    while( np->next != np)
    {
        ap1 = np->atom1;  ap2 = np->atom2;
        if( ap1->active && ap2->active ){
            if( ap1->x != 0. && ap1->y != 0. && ap1->z != 0.){
                if( ap2->x != 0. && ap2->y != 0. && ap2->z != 0.){
                    dx = ap1->x - ap2->x;
                    dy = ap1->y - ap2->y;
                    dz = ap1->z - ap2->z;
                    r = dx*dx + dy*dy + dz*dz;
                    r = sqrt(r);
                    np->d = r;
                    if( dm > 0. ) np->dm = dm;
                    if( dh > 0. ) np->dh = dh;
                }}}
        np = np->next;
        if( np == NULL )break;
    }
    if( np == NULL ) return 0;

    ap1 = np->atom1;  ap2 = np->atom2;
    if( ap1->active && ap2->active ){
        if( ap1->x != 0. && ap1->y != 0. && ap1->z != 0.){
            if( ap2->x != 0. && ap2->y != 0. && ap2->z != 0.){
                dx = ap1->x - ap2->x;
                dy = ap1->y - ap2->y;
                dz = ap1->z - ap2->z;
                r = dx*dx + dy*dy + dz*dz;
                r = sqrt(r);
                np->d = r;
                if( dm > 0. ) np->dm = dm;
                if( dh > 0. ) np->dh = dh;
            }}}

    return 1;}



/* v_noel()
* this function sums up the potentials
* for the atoms defined in the NOEL data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int v_noel( V, lambda )
float *V,lambda;
{
    NOEL *bp;
    float r,xt,yt,zt;
    ATOM *a1,*a2;


    bp = noel_first;
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
            bp->ismet = (1==1);
            if( r < bp->d -bp->dm)
            {
                bp->ismet = (1==0);
                r = r - bp->d + bp->dm;
                *V += bp->km * r*r;
            } else if( r > bp->d+ bp->dh) {
                bp->ismet = (1==0);
                r = r - bp->d - bp->dh;
                *V += bp->kh * r*r;
            }
        }
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
int v_noe_lin( V, lambda )
float *V,lambda;
{
    NOEL *bp;
    float r,xt,yt,zt;
    ATOM *a1,*a2;


    bp = noel_first;
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
            bp->ismet = (1==1);
            if( r < bp->d -bp->dm)
            {
                bp->ismet = (1==0);
                r = r - bp->d + bp->dm;
                if( r < -bp->dm) r = -bp->dm;
                *V += bp->km * r*r;
            } else if( r > bp->d+ bp->dh) {
                bp->ismet = (1==0);
                r = r - bp->d - bp->dh;
                if( r > bp->dh) r =  bp->dh;
                *V += bp->kh * r*r;
            }
        }
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* f_noel()
*
* f_noel increments the forces in the atom structures by the force
* due to the noel components.  NOTE THE WORD increment.
* the forces should first be zero'd.
* if not then this code will be invalid.  THIS IS DELIBERATE.
* on bigger (and better?) machines the different potential terms
* may be updated at random or in parrellel, if we assume that this routine
* will initialize the forces then we can't do this.
*/
int f_noel(lambda)
float lambda;
/*  returns 0 if error, 1 if OK */
{
    NOEL *bp;
    float r,t,ux,uy,uz;
    ATOM *a1,*a2;


    bp = noel_first;
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
            bp->ismet = (1==1);
            if( r < bp->d -bp->dm)
            {
                bp->ismet = (1==0);
                r = r - bp->d + bp->dm;
                ux = 2*bp->km * r *ux;
                uy = 2*bp->km * r *uy;
                uz = 2*bp->km * r *uz;
            } else if( r > bp->d+ bp->dh) {
                bp->ismet = (1==0);
                r = r - bp->d - bp->dh;
                ux = 2*bp->kh * r *ux;
                uy = 2*bp->kh * r *uy;
                uz = 2*bp->kh * r *uz;
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
int f_noe_lin(lambda)
float lambda;
/*  returns 0 if error, 1 if OK */
{
    NOEL *bp;
    float r,t,ux,uy,uz;
    ATOM *a1,*a2;


    bp = noel_first;
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
            bp->ismet = (1==1);
            if( r < bp->d -bp->dm)
            {
                bp->ismet = (1==0);
                r = r - bp->d + bp->dm;
                if( r < -bp->dm) r = -bp->dm;
                ux = 2.*bp->km*r  *ux;
                uy = 2.*bp->km*r  *uy;
                uz = 2.*bp->km*r  *uz;
            } else if( r > bp->d+ bp->dh) {
                bp->ismet = (1==0);
                r = r - bp->d - bp->dh;
                if( r > bp->dh) r =  bp->dh;
                ux = 2.*bp->kh*r  *ux;
                uy = 2.*bp->kh*r  *uy;
                uz = 2.*bp->kh*r  *uz;
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
/* function get_noel( a1,noeled,10,innoel);
* check the NOELS list for atoms noeled to a1
*/
void get_noel( a1,noeled,mnoel,innoel)
ATOM *a1, *noeled[];
int mnoel,*innoel ;
{
    NOEL *mine;
    mine = noel_first;
    *innoel = 0;
    while(1)
    {
        if( (mine == NULL) )
        {
            return;
        }
        if( mine->atom1 == a1)
        {
            noeled[(*innoel)++] = mine->atom2;
        }
        if( mine->atom2 == a1)
        {
            noeled[(*innoel)++] = mine->atom1;
        }
        if( mine == mine->next) return;
        mine = mine->next;
        if( *innoel == mnoel ) return;
    }
}
/* function get_noel( a1,noeled,10,innoel);
* check the NOELS list for atoms noeled to a1
*/
void get_noel_and_length( a1,noeled,r,mnoel,innoel)
ATOM *a1, *noeled[];
int mnoel,*innoel ;
float r[];
{
    NOEL *mine;
    mine = noel_first;
    *innoel = 0;
    while(1)
    {
        if( (mine == NULL) )
        {
            return;
        }
        if( mine->atom1 == a1)
        {
            r[*innoel] = mine->d;
            noeled[(*innoel)++] = mine->atom2;
        }
        if( mine->atom2 == a1)
        {
            r[*innoel] = mine->d;
            noeled[(*innoel)++] = mine->atom1;
        }
        if( mine == mine->next) return;
        mine = mine->next;
        if( *innoel == mnoel ) return;
    }
}
void get_noel_and_bounds( a1,noeled,r,rmax,mnoel,innoel)
ATOM *a1, *noeled[];
int mnoel,*innoel ;
float r[],rmax[];
{
    NOEL *mine;
    mine = noel_first;
    *innoel = 0;
    while(1)
    {
        if( (mine == NULL) )
        {
            return;
        }
        if( mine->atom1 == a1)
        {
            r[*innoel] = mine->d-mine->dm;
            rmax[*innoel] = mine->d+mine->dh;
            noeled[(*innoel)++] = mine->atom2;
        }
        if( mine->atom2 == a1)
        {
            r[*innoel] = mine->d-mine->dm;
            rmax[*innoel] = mine->d+mine->dh;
            noeled[(*innoel)++] = mine->atom1;
        }
        if( mine == mine->next) return;
        mine = mine->next;
        if( *innoel == mnoel ) return;
    }
}
/* noel_next()
*  like bond_next() but for NOEl structures
*/
int noel_next( int i, ATOM **n1, ATOM **n2, int *ismet )
{
    static NOEL *np = NULL ;
    *n1 = NULL ; *n2 = NULL;
    if( noel_first == NULL ) return (1==0);
if( np == NULL || i <= 0 ) { np = noel_first;}
    else{ np = np->next; }
    *n1 = np->atom1; *n2 = np->atom2;
    *ismet = np->ismet;
    if( np->next != np)return (1==1);
    return ( 1== 0);
}
/* routine dump_noels
* this function outputs the noel parameters
* and does it in a simple form
* noel ser1,ser2,k,req
* the rest is just free format
*/
void dump_noels( where )
FILE *where;
{
    NOEL *b;
    ATOM *a1,*a2;
    b = noel_first;
    if( b == NULL ) return;
    while( (b->next != b) )
    {
        if( b->next == NULL) return;
        a1 = b->atom1; a2 = b->atom2;
        fprintf( where,"noel %d %d %f %f %f %f %f;\n",a1->serial,a2->serial,
                 b->d, b->dm, b->dh, b->km, b->kh);
        b = b->next;
    }
    if( b->next == NULL) return;
    a1 = b->atom1; a2 = b->atom2;
    fprintf( where,"noel %d %d %f %f %f %f %f ;\n",a1->serial,a2->serial,
             b->d, b->dm, b->dh, b->km, b->kh);
}

/* a_noel()
* this function sums up the potentials
* for the atoms defined in the NOEL data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int a_noel( V, lambda,ilow,ihigh,op )
float *V,lambda;
int ilow,ihigh;
FILE *op;
{
    NOEL *bp;
    float r,xt,yt,zt;
    ATOM *a1,*a2;


    bp = noel_first;
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
            if( r < bp->d -bp->dm)
                zt= bp->km*( r - bp->d+ bp->dm)*(r - bp->d+ bp->dm);
            if( r > bp->d +bp->dh)
                zt= bp->kh*( r - bp->d- bp->dh)*(r - bp->d- bp->dh);
            *V += zt;
            fprintf(op,"NOEl %s %d %s %d E %f value %f error %f\n"
                    ,a1->name,a1->serial,a2->name,a2->serial,zt,r,r-bp->d);
        }
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* gsdg_noel( ATOM *ap )
*  
* setup the distances for NOEL terms
*/
int gsdg_noel( ap)
ATOM *ap;
{
    ATOM *bp;
    NOEL *np;

    np = noel_first;
    while(1)
    { if( np == NULL ) return;
        if( np->atom1 == ap )
        {  bp = np->atom2; bp->vx = (np->d*np->d );
            bp->vy = np->km; }
        if( np->atom2 == ap )
        {  bp = np->atom1; bp->vx = (np->d*np->d );
            bp->vy = np->km; }

        if( np == np->next ) return;
        np = np->next;
    }
}
/* v_ho_noel()
* this function sums up the potentials
* for the atoms defined in the NOEL data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int v_ho_noel( V, lambda )
float *V,lambda;
{
    NOEL *bp;
    float r,xt,yt,zt;
    ATOM *a1,*a2;
    float hol,get_f_variable(),target;

    hol = get_f_variable("lambda");

    bp = noel_first;
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
            bp->ismet = (1==1);
            if( r < bp->d -bp->dm)
            {
                bp->ismet = (1==0);
                target = hol*r + (1.-hol)*(bp->d-bp->dm);
                r = r - target;
                *V += bp->km * r*r;
            } else if( r > bp->d+ bp->dh) {
                bp->ismet = (1==0);
                target = hol*r + (1.-hol)*(bp->d+bp->dh);
                r = r - target;
                *V += bp->kh * r*r;
            }
        }
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* f_ho_noel()
*
* f_ho_noel increments the forces in the atom structures by the force
* due to the noel components.  NOTE THE WORD increment.
* the forces should first be zero'd.
* if not then this code will be invalid.  THIS IS DELIBERATE.
* on bigger (and better?) machines the different potential terms
* may be updated at random or in parrellel, if we assume that this routine
* will initialize the forces then we can't do this.
*/
int f_ho_noel(lambda)
float lambda;
/*  returns 0 if error, 1 if OK */
{
    NOEL *bp;
    float r,k,ux,uy,uz;
    ATOM *a1,*a2;
    float hol,get_f_variable(),target;

    hol = get_f_variable("lambda");


    bp = noel_first;
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
            r = ux*ux + uy*uy + uz*uz;
            /* watch for FP errors*/
            bp->ismet = (1==1);
            if( r <= 1.e-5)
            { r = 0; ux = 1.; uy = 0.; uz = 0.; }else{
                r = sqrt(r); ux = ux/r; uy = uy/r; uz = uz/r;
            }
            if( r < bp->d -bp->dm)
            {
                bp->ismet = (1==0);
                target = hol*r + (1.-hol)*(bp->d-bp->dm);
                r = r - target;
                ux = 2*bp->km * r*(1.-hol) *ux;
                uy = 2*bp->km * r*(1.-hol) *uy;
                uz = 2*bp->km * r*(1.-hol) *uz;
            } else if( r > bp->d+ bp->dh) {
                bp->ismet = (1==0);
                target = hol*r + (1.-hol)*(bp->d+bp->dh);
                r = r - target;
                ux = 2*bp->kh * r*(1.-hol) *ux;
                uy = 2*bp->kh * r*(1.-hol) *uy;
                uz = 2*bp->kh * r*(1.-hol) *uz;
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
