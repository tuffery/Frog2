/* swarm.c
*
*
*  A set of distances is restrained to the mean value
*  not quite sure yet wether to use sigma weights
*  
*  this is another way of enforcing regularity in polymers
*  and may well be easier and better than the SPIN ideas.
*  It is mathematically a similar concept, but a rather
*  different implementation.
*
*
* collection of routines to service SWARM length potentials
*
* POOP (Poor-mans Object Oriented Programming) using scope rules
*
* these routines hold a data base (in terms of array indeces)
* of SWARM bonds, with the associated length and force constants
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
*  copyright 1999 Robert W. Harrison
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
* but otherwise is self-contained. Note the hooks for Non-SWARMed potentials
*/
#define NSWARM 8
typedef struct{
    ATOM *atoms[NSWARM];
    float k;
    int n;
    void *next;
}  SWARM;
#define SLONG sizeof(SWARM)

SWARM *SWARM_first = NULL;
SWARM *SWARM_last = NULL;
/* function SWARM adds a SWARM to the SWARM list
* returns 1 if ok
* returns 0 if not
*  is passed the atom serial numbers, length and constant
* allocates the new_p memory, initializes it and
* returns
*/
int swarm( k,n,p1,p2,p3,p4,p5,p6,p7,p8)
int n,p1,p2,p3,p4,p5,p6,p7,p8;
float k ;
{
    ATOM *a_m_serial();
    SWARM *new_p;
    int i,ip[8];

    /* get the atom pointers for the two serial numbers */
    if( ( new_p = malloc( SLONG ) ) == NULL)
    {
        return 0;
    }

    ip[0] = p1;
    ip[1] = p2;
    ip[2] = p3;
    ip[3] = p4;
    ip[4] = p5;
    ip[5] = p6;
    ip[6] = p7;
    ip[7] = p8;
    for( i=0; i< n; i++)
    {
        new_p->atoms[i] = a_m_serial(ip[i]);
        if( new_p->atoms[i] == NULL)
        {
            aaerror("undefined atom in SWARM term;\n");
            free(new_p);
            return 0;
        }
    }

    /* initialize the pointers */
    if( SWARM_first == NULL) SWARM_first = new_p;
    if( SWARM_last == NULL) SWARM_last = new_p;

    new_p -> k = k;
    new_p -> n = n;
    new_p -> next = new_p;
    SWARM_last -> next = new_p;
    SWARM_last = new_p;
    return 1;
}


/* v_swarm()
* this function sums up the potentials
* for the atoms defined in the SWARM data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int v_swarm( V, lambda )
float *V,lambda;
{
    SWARM *bp;
    float r,xt,yt,zt;
    float rmean;
    float dx,dy,dz;
    ATOM *a1,*a2;
    int i;


    bp = SWARM_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        xt = 0.; yt = 0.; zt = 0.;
        rmean = 0.;
        for( i=0; i< bp->n; i+=2)
        {
            a1 = bp->atoms[i];
            a2 = bp->atoms[i+1];
            if( lambda == 0.)
            {
                dx = a2->x - a1->x;
                dy = a2->y - a1->y;
                dz = a2->z - a1->z;
            } else	{
                dx = (a2->x -a1->x +lambda*(a2->dx-a1->dx));
                dy = (a2->y -a1->y +lambda*(a2->dy-a1->dy));
                dz = (a2->z -a1->z +lambda*(a2->dz-a1->dz));
            }
            rmean = rmean+ sqrt(dx*dx + dy*dy + dz*dz);
        }/* for i */
        rmean = (rmean+rmean)/bp->n;
        for( i=0; i< bp->n ; i+= 2)
        {
            a1 = bp->atoms[i];
            a2 = bp->atoms[i+1];
            if( lambda == 0.)
            {
                dx = a2->x - a1->x;
                dy = a2->y - a1->y;
                dz = a2->z - a1->z;
            } else	{
                dx = (a2->x -a1->x +lambda*(a2->dx-a1->dx));
                dy = (a2->y -a1->y +lambda*(a2->dy-a1->dy));
                dz = (a2->z -a1->z +lambda*(a2->dz-a1->dz));
            }

            r = sqrt(dx*dx +dy*dy + dz*dz);
            *V += bp->k*(rmean-r)*(rmean -r);
        }/* i */
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}

/* f_swarm()
*
* f_swarm increments the forces in the atom structures by the force
* due to the SWARM components.  NOTE THE WORD increment.
* the forces should first be zero'd.
* if not then this code will be invalid.  THIS IS DELIBERATE.
* on bigger (and better?) machines the different potential terms
* may be updated at random or in parrellel, if we assume that this routine
* will initialize the forces then we can't do this.
*/
int f_swarm(lambda)
float lambda;
/*  returns 0 if error, 1 if OK */
{
    SWARM *bp;
    float xt,yt,zt,dx,dy,dz;
    float r,k,ux,uy,uz,dV,rmean;
    ATOM *a1,*a2;
    int i;


    bp = SWARM_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;




        xt = 0.; yt = 0.; zt = 0.;
        rmean = 0.;
        for( i=0; i< bp->n; i+=2)
        {
            a1 = bp->atoms[i];
            a2 = bp->atoms[i+1];
            if( lambda == 0.)
            {
                dx = a2->x - a1->x;
                dy = a2->y - a1->y;
                dz = a2->z - a1->z;
            } else	{
                dx = (a2->x -a1->x +lambda*(a2->dx-a1->dx));
                dy = (a2->y -a1->y +lambda*(a2->dy-a1->dy));
                dz = (a2->z -a1->z +lambda*(a2->dz-a1->dz));
            }
            rmean += sqrt(dx*dx + dy*dy + dz*dz);
        }/* i */
        rmean = (rmean + rmean )/bp->n;
        k = 2.*bp->k*(1.-1./(float)bp->n);
        for( i=0; i < bp->n; i+= 2)
        {
            a1 = bp->atoms[i];
            a2 = bp->atoms[i+1];
            if( lambda == 0.)
            {
                dx = a2->x - a1->x;
                dy = a2->y - a1->y;
                dz = a2->z - a1->z;
            } else	{
                dx = (a2->x -a1->x +lambda*(a2->dx-a1->dx));
                dy = (a2->y -a1->y +lambda*(a2->dy-a1->dy));
                dz = (a2->z -a1->z +lambda*(a2->dz-a1->dz));
            }
            xt = dx; yt = dy; zt = dz;
            /*	*V += bp->k*( target -r)*(target -r); */
            r = sqrt(xt*xt + yt*yt + zt*zt);
            dV = k *(r-rmean );
            if( r > 1.e-7) r = 1./r;
            ux = dV*xt*r;
            uy = dV*yt*r;
            uz = dV*zt*r;

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

        }/* i */

        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}

/* routine dump_swarms
* this function outputs the SWARM parameters
* and does it in a simple form
* SWARM ser1,ser2,k,req
* the rest is just free format
*/
void dump_swarms( where )
FILE *where;
{
    SWARM *b;
    ATOM *a1;
    int i;
    b = SWARM_first;
    if( b == NULL ) return;
    while( (1==1) )
    {
        if( b->next == NULL) return;

        fprintf( where,"SWARM %f %d ",b->k,b->n);
        for( i=0; i< b->n; i++)
        {
            a1 = b->atoms[i];
            fprintf(where,"%d ",a1->serial);
        }
        fprintf(where,";\n");
        if( b == b->next) return;
        b = b->next;
    }
}

int a_swarm(V,lambda,ilow,ihigh,op)
float *V,lambda;
int ilow,ihigh;
FILE *op;
{
    SWARM *bp;
    ATOM *a1,*a2;
    float xt,yt,zt,dx,dy,dz,rmean,r;
    int i,j;

    bp = SWARM_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        j = (1==0);
        for( i=0; i< bp->n; i++)
        {
            a1 = bp->atoms[i];
            if( a1->serial >= ilow && a1->serial <= ihigh)
            { j = (1==1); break;}
        }
        if( j) {
            xt = 0.; yt = 0.; zt = 0.;
            rmean = 0.;
            for( i=0; i< bp->n; i+=2)
            {
                a1 = bp->atoms[i];
                a2 = bp->atoms[i+1];
                if( lambda == 0.)
                {
                    dx = a2->x - a1->x;
                    dy = a2->y - a1->y;
                    dz = a2->z - a1->z;
                } else	{
                    dx = (a2->x -a1->x +lambda*(a2->dx-a1->dx));
                    dy = (a2->y -a1->y +lambda*(a2->dy-a1->dy));
                    dz = (a2->z -a1->z +lambda*(a2->dz-a1->dz));
                }
                rmean += sqrt(dx*dx + dy*dy + dz*dz);
            }/* for i */

            rmean = (rmean + rmean)/bp->n;
            for( i=0; i< bp->n; i+=2)
            {
                if( lambda == 0.)
                {
                    dx = a2->x - a1->x;
                    dy = a2->y - a1->y;
                    dz = a2->z - a1->z;
                } else	{
                    dx = (a2->x -a1->x +lambda*(a2->dx-a1->dx));
                    dy = (a2->y -a1->y +lambda*(a2->dy-a1->dy));
                    dz = (a2->z -a1->z +lambda*(a2->dz-a1->dz));
                }
                r = sqrt(dx*dx + dy*dy + dz*dz);
                yt = bp->k *( r- rmean)*(r - rmean);
                *V += yt;
                xt += yt;
            }/* i */
            fprintf(op,"SWARM ");
            for( i=0; i< bp->n; i++)
            { a1 = bp->atoms[i]; fprintf(op,"%d ",a1->serial); }
            fprintf(op,"Mean distance %f energy %f\n",rmean,xt);

        }/* if (j) */
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }

}
