/* morse.c
*
* collection of routines to service morse length potentials
*
* POOP (Poor-mans Object Oriented Programming) using scope rules
*
* these routines hold a data base (in terms of array indeces)
* of morse bonds, with the associated length and force constants
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
* but otherwise is self-contained. Note the hooks for Non-morseed potentials
*/
typedef struct{
    ATOM *atom1,*atom2;
    float length,k,order;
    void *next;
}  MORSE;
#define MLONG sizeof(MORSE)

MORSE *morse_first = NULL;
MORSE *morse_last = NULL;
/* function morse adds a morse to the morse list
* returns 1 if ok
* returns 0 if not
*  is passed the atom serial numbers, length and constant
* allocates the new memory, initializes it and
* returns
*/
int morse( p1,p2,bl,fk,order)
int p1,p2;
float bl,fk ,order;
{
    ATOM *ap1,*ap2,*a_m_serial();
    MORSE *new;
    char line[BUFSIZ];
    /* get the atom pointers for the two serial numbers */
    ap1 = a_m_serial( p1 );
    ap2 = a_m_serial( p2 );
    if( (ap1 == NULL) || (ap2 == NULL) )
    {
        sprintf( line,"undefined atom in morse %d %d \0",p1,p2);
        aaerror( line );
        return 0;
    }

    if( ( new = malloc( MLONG ) ) == NULL)
    {
        return 0;
    }
    /* initialize the pointers */
    if( morse_first == NULL) morse_first = new;
    if( morse_last == NULL) morse_last = new;
    new -> atom1 = ap1;
    new -> atom2 = ap2;
    new -> length = bl;
    new -> k = fk;
    new -> order = order;
    new -> next = new;
    morse_last -> next = new;
    morse_last = new;
    return 1;
}


/* v_morse()
* this function sums up the potentials
* for the atoms defined in the MORSE data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int v_morse( V, lambda )
float *V,lambda;
{
    MORSE *bp;
    float r,xt,yt,zt;
    ATOM *a1,*a2;


    bp = morse_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1; a2 = bp->atom2;
        if( a1->active || a2->active){
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
            /* *V += bp->k*( r - bp->length)*(r - bp->length);
            */
            xt = 1.- exp( -(bp->order)*(r - bp->length)) ;
            *V += bp->k*xt*xt;
        }
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* f_morse()
*
* f_morse increments the forces in the atom structures by the force
* due to the morse components.  NOTE THE WORD increment.
* the forces should first be zero'd.
* if not then this code will be invalid.  THIS IS DELIBERATE.
* on bigger (and better?) machines the different potential terms
* may be updated at random or in parrellel, if we assume that this routine
* will initialize the forces then we can't do this.
*/
int f_morse(lambda)
float lambda;
/*  returns 0 if error, 1 if OK */
{
    MORSE *bp;
    float r,k,ux,uy,uz;
    ATOM *a1,*a2;


    bp = morse_first;
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
                ux = (a1->x -a2->x +lambda*(a1->dx-a2->dx));
                uy = (a1->y -a2->y +lambda*(a1->dy-a2->dy));
                uz = (a1->z -a2->z +lambda*(a1->dz-a2->dz));
            }
            r = ux*ux + uy*uy + uz*uz;
            /* watch for FP errors*/
            if( r <= 1.e-5)
            { r = 0; ux = 1.; uy = 0.; uz = 0.; }else{
                r = sqrt(r); ux = ux/r; uy = uy/r; uz = uz/r;
            }
            k =  exp( -(bp->order)*(r - bp->length)) ;
            k = -2*bp->order *(1.-k)*k* bp->k;
            ux = 2*k*(r-bp->length)*ux;
            uy = 2*k*(r-bp->length)*uy;
            uz = 2*k*(r-bp->length)*uz;
            if( a1->active){
                a1->fx += ux;
                a1->fy += uy;
                a1->fz += uz;
            }
            if( a2->active){
                a2->fx -= ux;
                a2->fy -= uy;
                a2->fz -= uz;
            }
        }
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* function get_morse( a1,morseed,10,inmorse);
* check the MORSES list for atoms morseed to a1
*/
void get_morse( a1,morseed,mmorse,inmorse)
ATOM *a1, *morseed[];
int mmorse,*inmorse ;
{
    MORSE *mine;
    mine = morse_first;
    *inmorse = 0;
    while(1)
    {
        if( (mine == NULL) )
        {
            return;
        }
        if( mine->atom1 == a1)
        {
            morseed[(*inmorse)++] = mine->atom2;
        }
        if( mine->atom2 == a1)
        {
            morseed[(*inmorse)++] = mine->atom1;
        }
        if( mine == mine->next) return;
        mine = mine->next;
        if( *inmorse == mmorse ) return;
    }
}
/* routine dump_morses
* this function outputs the morse parameters
* and does it in a simple form
* morse ser1,ser2,k,req
* the rest is just free format
*/
void dump_morse( where )
FILE *where;
{
    MORSE *b;
    ATOM *a1,*a2;
    b = morse_first;
    if( b == NULL ) return;
    while( (b->next != b) )
    {
        if( b->next == NULL) return;
        a1 = b->atom1; a2 = b->atom2;
        fprintf( where,"morse %d %d %f %f %f \;\n",a1->serial,a2->serial,
                 b->length,b->k,b->order);
        b = b->next;
    }
    if( b->next == NULL) return;
    a1 = b->atom1; a2 = b->atom2;
    fprintf( where,"morse %d %d %f %f %f\;\n",a1->serial,a2->serial,
             b->length,b->k,b->order);
}
