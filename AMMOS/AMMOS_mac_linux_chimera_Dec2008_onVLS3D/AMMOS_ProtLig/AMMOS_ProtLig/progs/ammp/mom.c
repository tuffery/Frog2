/* mom.c
*
* evaluate charges by the method of moments
*
* largely cadged from Rappe and Goddard JPC 95 3358-3363
* cleaned up the solver 
*
* modifies the atom structure to have a jaa and chi field
*  (self colomb energy and electronegativity )
*/
/*
*  copyright 1993 Robert W. Harrison
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
#define MAXMOM 100
ATOM *mom_list[MAXMOM];
int in_mom_list = 0;

void mom( op,tq,niter)
FILE *op;
float tq;
int niter;
{
    float (*v)[], (*Z)[],(*Zk)[];
    float r,dx,dy,dz;
    float dv,tv;
    int i,j,row;
    float mom_jab();

    if( in_mom_list == 0 ) return;
    if( niter <= 0 ) niter = 20; /* default number of trials */
    v = malloc( in_mom_list * sizeof(float)) ;
    Z = malloc( in_mom_list*in_mom_list * sizeof(float)) ;
    Zk = malloc( in_mom_list*in_mom_list * sizeof(float)) ;
    if( v == NULL || Z == NULL|| Zk == NULL)
    { aaerror(" cannot allocate memory in mom\n"); return;}
    row = in_mom_list; /* this is just to save typing */
    /* build the impedance matrix in Zk */
    for( i=0; i< row; i++)
        for( j=i+1; j< row; j++)
        {
            dx = mom_list[j]->x - mom_list[i]->x;
            dy = mom_list[j]->y - mom_list[i]->y;
            dz = mom_list[j]->z - mom_list[i]->z;
            r = sqrt(dx*dx + dy*dy + dz*dz);
            (*Zk)[i + j*row] = mom_jab(r,mom_list[i]->jaa,mom_list[j]->jaa);
            (*Zk)[j + i*row] = (*Zk)[i + j*row];
        }
    for( i=0; i< row; i++)
    { (*Zk)[i+i*row] = mom_list[i]->jaa/2 ; }
    /* now we're ready to do it */

    dv = 0; /* the offset charge */

    for( j=0; j< niter; j++)
    {
        for( i=0;i< row*row ; i++)
            (*Z)[i] = (*Zk)[i] ;


        for( i=0; i< row; i++)
        { (*v)[i]  =  -mom_list[i]->chi; }

        tv = 0.;
        for( i=0; i< row; i++)
        { tv += (*v)[i]; }
        tv = tv/row +dv;
        for( i=0; i< row; i++)
            (*v)[i] -= tv;

        mom_solve( Z,v,row,row);

        tv = 0.;
        for( i=0; i< row; i++)
        {
            tv += (*v)[i];
        }



        /*	dv += 2*(tv - tq); /*
        /*	dv += 14.4*(tv - tq)/row*1.5;
        */
        dv += 14.4*(tv - tq)/row*.25;
        fprintf( op,"MoM iter %d error %e\n",j,tv-tq);

    }/* end of j for */

    for( i=0; i< row; i++)
    {
        mom_list[i]->q = (*v)[i];
    }
    /* cleanup */
    free(Zk) ;free( Z); free(v); in_mom_list = 0;
}

void mom_add( s1,s2)
int s1,s2;
{
    int i,j;
    ATOM *ap,*a_m_serial(),*a_next();
    int a_number(),numatm;
    numatm = a_number();
    if( numatm == 0 ) return ;
    if( s2 > 0 && s1 > s2)
    {  i = s1; s1 = s2; s2 = i; }
    if( s2 > s1)
    {
        for( i=0; i< numatm; i++)
        {
            ap = a_next(i);
            if(  ap->serial >=s1 && ap->serial <= s2 )
            {
                if( ap->chi > 0. && ap->jaa > 0.)
                {
                    for( j=0; j< in_mom_list; j++)
                    { if( mom_list[j] == ap) goto THERE_NOW; }
                    mom_list[in_mom_list++ ] = ap;
THERE_NOW:
                    if( in_mom_list == MAXMOM) in_mom_list --;
                }
            }
        }
        return;
    }
    if( (ap = a_m_serial(s1)) != NULL)
    {
        if( ap->chi > 0. && ap->jaa > 0.)
        {
            for( j=0; j< in_mom_list; j++)
            {if( ap == mom_list[j]) goto THERE_AGAIN;}
            mom_list[in_mom_list++ ] = ap;
THERE_AGAIN:
            if( in_mom_list == MAXMOM) in_mom_list --;
        }
    }
}

void mom_param( serial,chi,jaa )
int serial;
float chi,jaa;
{
    ATOM *ap,*a_m_serial();

    if( (ap = a_m_serial(serial)) == NULL)
    { aaerror(" MOM> cannot modify non-extant atom "); return;}
    ap->chi = chi;
    ap->jaa = jaa;
}

float mom_jab( r,j1,j2 )
float r;
float j1,j2;
{
    float a,b,b2,b3;
    /* fit by guess to the repulsion curve
    *  tested on methane 
    */
    /*	if( r < 4.) return 9/(1+ .09375*r*r);
    */
    /*
    	if( r < 4.) return 8.5/(1+ .085069444*r*r)/2;
    */
    if( r < 30.) {
        /* taken from wallace h,h interaction with a small fudge */
        /* the energy of h(1s) h(1s) is given, and we scale it
           by a term depending on atom type */
        /*		a = (j1+j2)/4.;
        */
        a = sqrt(j1*j2)/2.;
        /*		r = 1.08*r;
        */
        /* effective radius is adjusted here  (1. is ok this is a little better) */
        r = 1.1*r;
        b = exp( -r );
        /*		b2 = 1. +  (33*r + 9*r*r + r*r*r)/48.;
        */
        /*		b2 = 1. +  (20*r - 9*r*r - r*r*r)/48.;
        */
        b2 = 1. +  (15*r - 9*r*r - r*r*r)/48.;
        return a/r*(b*b2);

    }

    return 14.4/r/2;
}


/*  this is a routine to solve a linear equation by
    guassian elimination.  (basically solve.for translated) */
/* in order to have the  array matrix be of any length it must be passed as
   a linear array.  Since C has the opposite convention for array packing from 
   FORTRAN ( row fastest rather than column fastest) the leading dimension
   ilead is the row size of the array to which matrix points */
mom_solve( matrix,vector,irow,ilead )
int irow,ilead;
float (*matrix)[];
float (*vector)[];

{
    float quotient;
    int i,j,k;
    int  mpi,mpj,mpk;
    mpi = 0;
    for ( i=0 ; i < irow - 1 ; i++ )
    {  /* for each row */
        j = i ;
        mpj = mpi;
        while ( (*matrix) [mpi + i] == 0)
        {
            if( j == irow )
            { return (-1); }
            j ++;
            mpj += ilead;
            (*vector)[i] +=  (*vector)[j];
            for (k = i; k <irow  ; k++)
            {(*matrix)[mpi + k] += (*matrix)[mpj +k]; }
        }
        /* if here then the diagonal element is not zero so we can do the division*/
        mpj = mpi +ilead ;
        for ( j= i+1; j < irow ; j++ )
        {
            if( (*matrix)[mpj + i] != 0 )
            {
                quotient = (*matrix)[mpj + i]/(*matrix)[mpi + i];
                (*vector)[j] -= (*vector)[i]*quotient;
                for ( k=i ; k< irow ; k++ )
                { (*matrix)[mpj + k] -= (*matrix)[mpi + k]*quotient; }
            }  /* if */
            mpj += ilead;
        } /* for j */
        mpi += ilead;
    } /* for i */

    /* now start the back substitution loop */
    mpi = 0;
    for ( i = 0; i < irow - 1 ; i++ )
    {
        k = irow - i - 1;
        mpj= 0;
        mpk =  k*ilead;
        for ( j = 0; j < k ; j++)
        { (*vector)[j] -=(*matrix)[mpj+k]/(*matrix)[mpk+k]*(*vector)[k];
            mpj +=ilead; }
    } /* i */

    /* and finally divide by the diagonal elements */
    mpi = 0;
    for ( i=0; i <irow ; i++ )
    { (*vector)[i] /= (*matrix)[mpi + i];
        mpi += ilead;    }
    return (0);
}
