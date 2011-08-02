/* bonds.c
*
* collection of routines to service bond length potentials
*
* POOP (Poor-mans Object Oriented Programming) using scope rules
*
* these routines hold a data base (in terms of array indeces)
* of bonds, with the associated length and force constant
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
*  copyright 1992,1993,1994,1995 Robert W. Harrison
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
    float length,k;
    float dlength,dk; /* for abc terms */
    int ndiff;
    void *next;
}  BOND;
#define BLONG sizeof(BOND)

BOND *bond_first = NULL;
BOND *bond_last = NULL;
/* function bond adds a bond to the bond list
* returns 1 if ok
* returns 0 if not
*  is passed the atom serial numbers, length and constant
* allocates the new memory, initializes it and
* returns
*/
int bond( p1,p2,bl,fk)
int p1,p2;
float bl,fk ;
{
    ATOM *ap1,*ap2,*a_m_serial();
    BOND *new;
    int i;
    char line[BUFSIZ];
    /* get the atom pointers for the two serial numbers */
    ap1 = a_m_serial( p1 );
    ap2 = a_m_serial( p2 );
    if( (ap1 == NULL) || (ap2 == NULL) )
    {
        sprintf( line,"undefined atom in bond %d %d \0",p1,p2);
        aaerror( line );
        return 0;
    }

    if( ( new = malloc( BLONG ) ) == NULL)
    {
        return 0;
    }
    /* initialize the pointers */
    if( bond_first == NULL) bond_first = new;
    if( bond_last == NULL) bond_last = new;
    new -> atom1 = ap1;
    new -> atom2 = ap2;
    new -> length = bl;
    new -> k = fk;
    new -> next = new;
    /* update the exclude list in the atoms structure */
    if( ap1->dontuse < NEXCLUDE)
    {
        for( i=0; i< ap1->dontuse; i++)
            if( ap1->excluded[i] == ap2) goto excluded1;
        ap1->excluded[ap1->dontuse] = ap2; (ap1->dontuse)++;
    }else{
        aaerror(" too many bonds to an atom increase NEXCLUDE in ammp.h");
        exit(0);
    }
excluded1:
    if( ap2->dontuse < NEXCLUDE)
    {
        for( i=0; i< ap2->dontuse; i++)
            if( ap2->excluded[i] == ap1) goto excluded2;
        ap2->excluded[ap2->dontuse] = ap1; (ap2->dontuse)++;
    }else{
        aaerror(" too many bonds to an atom increase NEXCLUDE in ammp.h");
        exit(0);
    }
excluded2:
    bond_last -> next = new;
    bond_last = new;
    return 1;
}


/* v_bond()
* this function sums up the potentials
* for the atoms defined in the BOND data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int v_bond( V, lambda )
float *V,lambda;
{
    BOND *bp;
    float r,xt,yt,zt;
    ATOM *a1,*a2;


    bp = bond_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1; a2 = bp->atom2;
        if( a1->active || a2->active )
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
            r = sqrt(r); *V += bp->k*( r - bp->length)*(r - bp->length);
        }
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* f_bond()
*
* f_bond increments the forces in the atom structures by the force
* due to the bond components.  NOTE THE WORD increment.
* the forces should first be zero'd.
* if not then this code will be invalid.  THIS IS DELIBERATE.
* on bigger (and better?) machines the different potential terms
* may be updated at random or in parrellel, if we assume that this routine
* will initialize the forces then we can't do this.
*/
int f_bond(lambda)
float lambda;
/*  returns 0 if error, 1 if OK */
{
    BOND *bp;
    float r,k,ux,uy,uz;
    ATOM *a1,*a2;


    bp = bond_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        k = bp->k;
        a1 = bp->atom1; a2 = bp->atom2;
        if( a1->active || a2->active )
        {
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
            if( r <= 1.e-15)
            { r = 0; ux = 1.; uy = 0.; uz = 0.; }else{
                r = sqrt(r); ux = ux/r; uy = uy/r; uz = uz/r;
            }
            ux = 2*k*(r-bp->length)*ux;
            uy = 2*k*(r-bp->length)*uy;
            uz = 2*k*(r-bp->length)*uz;
            if( a1->active ){
                a1->fx += ux;
                a1->fy += uy;
                a1->fz += uz;
            }
            if( a2->active ){
                a2->fx -= ux;
                a2->fy -= uy;
                a2->fz -= uz;
            }
        }
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* function get_bond( a1,bonded,10,inbond);
* check the BONDS list for atoms bonded to a1
*/
void get_bond( a1,bonded,mbond,inbond)
ATOM *a1, *bonded[];
int mbond,*inbond ;
{
    BOND *mine;
    mine = bond_first;
    *inbond = 0;
    while(1)
    {
        if( (mine == NULL) )
        {
            return;
        }
        if( mine->atom1 == a1)
        {
            bonded[(*inbond)++] = mine->atom2;
        }
        if( mine->atom2 == a1)
        {
            bonded[(*inbond)++] = mine->atom1;
        }
        if( mine == mine->next) return;
        mine = mine->next;
        if( *inbond == mbond ) return;
    }
}
/* function get_bond_data( a1,bonded,10,inbond);
* check the BONDS list for atoms bonded to a1
*/
void get_bond_and_length( a1,bonded,r,mbond,inbond)
ATOM *a1, *bonded[];
float r[];
int mbond,*inbond ;
{
    BOND *mine;
    mine = bond_first;
    *inbond = 0;
    while(1)
    {
        if( (mine == NULL) )
        {
            return;
        }
        if( mine->atom1 == a1)
        {
            r[*inbond] = mine->length;
            bonded[(*inbond)++] = mine->atom2;
        }
        if( mine->atom2 == a1)
        {
            r[*inbond] = mine->length;
            bonded[(*inbond)++] = mine->atom1;
        }
        if( mine == mine->next) return;
        mine = mine->next;
        if( *inbond == mbond ) return;
    }
}
/* routine dump_bonds
* this function outputs the bond parameters
* and does it in a simple form
* bond ser1,ser2,k,req
* the rest is just free format
*/
void dump_bonds( where )
FILE *where;
{
    BOND *b;
    ATOM *a1,*a2;
    b = bond_first;
    if( b == NULL ) return;
    while( (b->next != b) )
    {
        if( b->next == NULL) return;
        a1 = b->atom1; a2 = b->atom2;
        fprintf( where,"bond %d %d %f %f \;\n",a1->serial,a2->serial,
                 b->length,b->k);
        b = b->next;
    }
    if( b->next == NULL) return;
    a1 = b->atom1; a2 = b->atom2;
    fprintf( where,"bond %d %d %f %f \;\n",a1->serial,a2->serial,
             b->length,b->k);
}

/* v_mmbond()
* this function sums up the potentials
* for the atoms defined in the BOND data structure.
*/
/* mm3 bond formula */
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int v_mmbond( V, lambda )
float *V,lambda;
{
    BOND *bp;
    float r,xt,yt,zt;
    ATOM *a1,*a2;


    bp = bond_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1; a2 = bp->atom2;
        if( a1->active || a2->active )
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
            /*	*V += bp->k*( r - bp->length)*(r - bp->length); */
            r = r - bp->length;
            *V += bp->k*r*r*(1.-2.55*r+7.*2.55/12*r*r);
        }
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* f_mmbond()
*
* mm3 bond formula
* f_mmbond increments the forces in the atom structures by the force
* due to the bond components.  NOTE THE WORD increment.
* the forces should first be zero'd.
* if not then this code will be invalid.  THIS IS DELIBERATE.
* on bigger (and better?) machines the different potential terms
* may be updated at random or in parrellel, if we assume that this routine
* will initialize the forces then we can't do this.
*/
int f_mmbond(lambda)
float lambda;
/*  returns 0 if error, 1 if OK */
{
    BOND *bp;
    float r,k,ux,uy,uz;
    ATOM *a1,*a2;


    bp = bond_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        k = bp->k;
        a1 = bp->atom1; a2 = bp->atom2;
        if( a1->active || a2->active )
        {
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
            if( r <= 1.e-5)
            { r = 0; ux = 1.; uy = 0.; uz = 0.; }else{
                r = sqrt(r); ux = ux/r; uy = uy/r; uz = uz/r;
            }
            /*	ux = 2*k*(r-bp->length)*ux;
            	uy = 2*k*(r-bp->length)*uy; 
            	uz = 2*k*(r-bp->length)*uz;
            	*V += bp->k*r*r*(1.-2.55*r+7.*2.55/12*r*r);
            */
            r = r - bp->length;
            ux = k*r*(2.-3*2.55*r+4*7*2.55*r*r/12)*ux;
            uy = k*r*(2.-3*2.55*r+4*7*2.55*r*r/12)*uy;
            uz = k*r*(2.-3*2.55*r+4*7*2.55*r*r/12)*uz;
            if( a1->active ){
                a1->fx += ux;
                a1->fy += uy;
                a1->fz += uz;
            }
            if( a2->active ){
                a2->fx -= ux;
                a2->fy -= uy;
                a2->fz -= uz;
            }
        }
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* a_bond()
* this function sums up the potentials
* for the atoms defined in the BOND data structure.
* only does bonds in the given range
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int a_bond( V, lambda,ilow,ihigh,op )
float *V,lambda;
int ilow,ihigh;
FILE *op;
{
    BOND *bp;
    float r,xt,yt,zt;
    ATOM *a1,*a2;


    bp = bond_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1; a2 = bp->atom2;
        if( a1->active || a2->active )
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
                r = sqrt(r);  zt = bp->k*( r - bp->length)*(r - bp->length);
                *V += zt;
                fprintf(op,"Bond %s %d %s %d E %f value %f error %f\n"
                        ,a1->name,a1->serial,a2->name,a2->serial,zt,r,r-bp->length);
            }
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* a_mmbond()
* this function sums up the potentials
* for the atoms defined in the BOND data structure.
*/
/* mm3 bond formula */
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int a_mmbond( V, lambda,ilow,ihigh,op )
float *V,lambda;
int ilow,ihigh;
FILE *op;
{
    BOND *bp;
    float r,xt,yt,zt;
    ATOM *a1,*a2;


    bp = bond_first;
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
            /*	*V += bp->k*( r - bp->length)*(r - bp->length); */
            r = r - bp->length;
            zt = bp->k*r*r*(1.-2.55*r+7.*2.55/12*r*r);
            *V += zt;
            fprintf(op,"mmBond %s %d %s %d E %f value %f error %f\n"
                    ,a1->name,a1->serial,a2->name,a2->serial,zt,r+bp->length,r);
            /*
            	fprintf(op,"mmBond %d %d  E %f error %f\n",a1->serial,a2->serial,zt,r);
            */
        }
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* bond_length()
* given two ATOMS find the bond length between them
* return -1. if not bonded
*/
float bond_length( a1,a2)
ATOM *a1, *a2;
{
    BOND *bp;
    bp = bond_first;
    while(1)
    {
        if( bp == NULL ) return -1.;
        if( a1 == bp->atom1 && a2 == bp->atom2) return bp->length;
        if( a2 == bp->atom1 && a1 == bp->atom2) return bp->length;
        if( bp == bp->next) return -1.;
        bp = bp->next;
    }
}

gsdg_bond( who)
ATOM *who;
{
    BOND *bp;
    ATOM *ap;

    bp = bond_first;
    while(1)
    { if( bp == NULL ) return;
        if( bp->atom1 == who )
        { ap = bp->atom2; ap->vx = bp->length*bp->length; ap->vy = bp->k;}
        if( bp->atom2 == who )
        { ap = bp->atom1; ap->vx = bp->length*bp->length; ap->vy = bp->k;}
        if( bp == bp->next) return;
        bp = bp->next;

    }
}


BOND *get_bond_pointer( a1,a2)
ATOM *a1, *a2;
{
    BOND *bp;

    bp = bond_first;
    while(1)
    {
        if( bp == NULL ) return;
        if( a1 == bp->atom1 && a2 == bp->atom2) return bp;
        if( a1 == bp->atom2 && a2 == bp->atom1) return bp;
        if( bp == bp->next) return NULL;
        bp = bp->next;
    }
    return NULL;
}

int bond_next( i,ap1,ap2)
int i;
ATOM **ap1,**ap2;
{
    static BOND *bp;
    if( i <= 0){ bp = bond_first;
        if( bp == NULL ){*ap1 = NULL; *ap2 == NULL; return 0; }}
    *ap1 = bp->atom1;
    *ap2 = bp->atom2;
    if( bp->next == bp || bp == NULL) return 0;
    if( bp->next != NULL)bp = bp->next;

    return 1;
}
/* v_ho_bond()
* this function sums up the potentials
* for the atoms defined in the BOND data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int v_ho_bond( V, lambda )
float *V,lambda;
{
    BOND *bp;
    float r,xt,yt,zt;
    ATOM *a1,*a2;
    float hol, get_f_variable();
    float target;

    hol = get_f_variable( "lambda");
    if( hol < 0. ) hol = 0.;
    if( hol > 1. ) hol = 1.;
    bp = bond_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1; a2 = bp->atom2;
        if( a1->active || a2->active )
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
            target = hol*r + (1.-hol)*bp->length;
            *V += bp->k*( r - target)*(r - target);
        }
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* f_ho_bond()
*
* f_bond increments the forces in the atom structures by the force
* due to the bond components.  NOTE THE WORD increment.
* the forces should first be zero'd.
* if not then this code will be invalid.  THIS IS DELIBERATE.
* on bigger (and better?) machines the different potential terms
* may be updated at random or in parrellel, if we assume that this routine
* will initialize the forces then we can't do this.
*/
int f_ho_bond(lambda)
float lambda;
/*  returns 0 if error, 1 if OK */
{
    BOND *bp;
    float r,k,ux,uy,uz;
    ATOM *a1,*a2;
    float hol, get_f_variable();
    float target;

    hol = get_f_variable( "lambda");
    if( hol < 0. ) hol = 0.;
    if( hol > 1. ) hol = 1.;
    bp = bond_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        k = bp->k;
        a1 = bp->atom1; a2 = bp->atom2;
        if( a1->active || a2->active )
        {
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
            if( r <= 1.e-15)
            { r = 0; ux = 1.; uy = 0.; uz = 0.; }else{
                r = sqrt(r); ux = ux/r; uy = uy/r; uz = uz/r;
            }
            target = hol*r + (1.-hol)*bp->length;
            ux = 2*k*(r-target)*(1.- hol)*ux;
            uy = 2*k*(r-target)*(1.- hol)*uy;
            uz = 2*k*(r-target)*(1.- hol)*uz;
            if( a1->active ){
                a1->fx += ux;
                a1->fy += uy;
                a1->fz += uz;
            }
            if( a2->active ){
                a2->fx -= ux;
                a2->fy -= uy;
                a2->fz -= uz;
            }
        }
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
