/*  gdock.c
*
*  genetic algorithm docker
*
* AMMP version 
*  given a range of atom ID's and uses them to define hull
*  also uses the total potential, 
*
*
*/  
/*
*  copyright 1998 Robert W. Harrison
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

/* specific definitions to pre-tune algorithm
* change these and re-compile
*/
#define MAX_PROGENY  4
#define MIN_GENE 8

/* genetic algorithm/quarternion rigid body solver */
void gdock( toler,ngene,niter,vara,varx, potent, inpotent, imin, imax)
float vara,varx,toler; /* variation in angle and xyz */
int ngene, niter, imin,imax, inpotent;
int (*potent[])();
{
    ATOM *(*ap)[],  *a_m_serial();
    int i,j,k,l,natom,nargs,mycase;
    /* quarternion_rot_tran() normalizes the quarternion and applies the
    * translation */
    void quarternion_rot_tran();

    /* store the genome here */
    float (*gene)[];
    float (*genet)[];
    /* store the energy and "fitness" */
    float (*fvals)[],(*fitness)[];
    float fbest,best[7];
    /* stud register*/
    int (*studbook)[];
    /* use these to analyze the values */
    int ibest, simplex_get_best();
    int worst, simplex_get_worst();
    float sigma, simplex_get_var();
    float randf(); /* rng assumed to be seeded */
    float x,y,z;
    int v_nonbon(),u_v_nonbon();

    if( ngene < MIN_GENE) ngene = MIN_GENE;
    /* allocate the genome  q1..q4,x,y,z X ngene storage */
    studbook = malloc(ngene*sizeof(int));
    fvals= malloc( ngene*sizeof(float));
    fitness = malloc( ngene*sizeof(float));
    gene = malloc( ngene*7*sizeof(float));
    genet = malloc(ngene*7*sizeof(float));

    if( genet == NULL || gene==NULL || fitness == NULL || fvals == NULL)
    {
        aaerror("gene cannot allocate memory in gdock\n");
        return ;
    }
    /* 	initialize the geneome
    *  always include the current position */
    (*gene)[0] = 0.;
    (*gene)[1] = 0.;
    (*gene)[2] = 0.;
    (*gene)[3] = 1.;
    (*gene)[4] = 0.;
    (*gene)[5] = 0.;
    (*gene)[6] = 0.;
    for (i=1; i< ngene; i++)
    {
        for( j=0; j< 3; j++)
        {
            k = i*7 + j;
            (*gene)[k] = (*gene)[j] + (2.*randf()-1.)*varx;
        }
        for( j=3; j< 7; j++)
        {
            k = i*7 + j;
            (*gene)[k] = (*gene)[j] + (2.*randf()-1.)*vara;
        }
    }
    /* figure out how many atoms there are in the range given */
    if( imin > imax ){ i= imin; imin = imax ; imax = i;}
    natom = 0;
    for( i=imin; i <= imax; i ++)
    {
        if( a_m_serial(i) != NULL) natom ++;
    }
    if( natom == 0) return ;

    nargs = 7;

    ap = malloc( (natom )*sizeof( ATOM *));
    if( ap == NULL)
    { fprintf(stderr," ap cannot allocate memory in gdock \n");
        return ;}
    /* now gather up the atoms */
    for( i=0; i< natom; i++) (*ap)[i] = NULL;
    j = 0;
    for( i=imin; i <= imax; i ++)
    {
        if( ((*ap)[j] =a_m_serial(i)) != NULL)
        {
            (*ap)[j]->gx = (*ap)[j]->x;
            (*ap)[j]->gy = (*ap)[j]->y;
            (*ap)[j]->gz = (*ap)[j]->z;
            j++;
        }
        if( j == natom) break;
    }

    /* initialize the fvals data structure */
    fbest = 10.e10;
    best[0] = 0.;
    best[1] = 0.;
    best[2] = 0.;
    best[3] = 1.;
    best[4] = 0.;
    best[5] = 0.;
    best[6] = 0.;
    for( i=0; i< ngene; i++)
    {
        /* how to do the potential */
        quarternion_rot_tran( &(*gene)[i*7],ap,natom);
        (*fvals)[i] = 0.;
        for( j=0; j < inpotent; j++)
        {
            if( (*potent[j]) != v_nonbon && (*potent[j]) != u_v_nonbon)
                (*potent[j])(&(*fvals)[i],0.);
            else
                zone_nonbon(&(*fvals)[i],0., ap, natom);
        }
        if( (*fvals)[i] < fbest){
            fbest = (*fvals)[i];
            for( j=0; j< 7; j++)
                best[j] = (*gene)[7*i+j];
        }
    }

    for( k=0; k< niter; k++ )
    {

        /* scane the range of values */
        /* the ngene-1 hack is because the simplex routines add 1 */
        ibest = simplex_get_best( fvals,ngene-1 );
        if( (sigma =simplex_get_var(fvals,ngene-1 )) < toler ) goto DONE;
        worst = simplex_get_worst(fvals,ngene-1 );

        printf(" %d best %d energy %f\n", k,ibest,(*fvals)[ibest]);
        printf(" %d worst %d energy %f\n", k,worst,(*fvals)[worst]);
        /* make up the new genome */
        /* evaluate "fitness" */
        x = (*fvals)[ibest];
        y = (*fvals)[worst];
        z = fabs(x-y);
        if( z > 0.) z = 1./z;
        /* linearly scale */
        x = 0;
        for( i=0; i< ngene; i++)
        {
            (*fitness)[i] = fabs((*fvals)[i] - y)*z;
            x += (*fitness)[i];
            (*studbook)[i] = 0; /* might as well here */
        }
        /* normalize */
        if( x > 0.) x = 1./x;
        for( i=0; i< ngene; i++)
            (*fitness)[i] *= x;
        /* convert to cumulative */
        /* dont , we just need deltas
        	x = 0.;
        	for( i=0; i< ngene; i++)
        	{
        		x += (*fitness)[i];
        		(*fitness)[i] = x;
        	}
        */
        /* now select new trials */
        for( i=0; i< ngene; i++)
        {
            /* setup with replacement */
REDO_THE_GENE: ;
            worst = 0;
            for( j=0; j< ngene; j++)
            {
                if( (*studbook)[j] < MAX_PROGENY)
                {
                    if( randf() <= (*fitness)[j])
                    {
                        worst ++;
                        (*genet)[i*7] = (*gene)[j*7] + (2.*randf()-1.)*varx*.1;
                        (*genet)[i*7+1] = (*gene)[j*7+1] + (2.*randf()-1.)*varx*.1;
                        (*genet)[i*7+2] = (*gene)[j*7+2] + (2.*randf()-1.)*varx*.1;
                        (*genet)[i*7+3] = (*gene)[j*7+3] + (2.*randf()-1.)*vara*.1;
                        (*genet)[i*7+4] = (*gene)[j*7+4] + (2.*randf()-1.)*vara*.1;
                        (*genet)[i*7+5] = (*gene)[j*7+5] + (2.*randf()-1.)*vara*.1;
                        (*genet)[i*7+6] = (*gene)[j*7+6] + (2.*randf()-1.)*vara*.1;
                        (*studbook)[j] += 1;
                        goto REPLACED;
                    }}
            } /* j */
            if( worst == 0) goto REDO_THE_GENE;
REPLACED: ;
        }/* i */
        /* now recombine a few */
        for( worst=0; worst < 2; worst ++)
        {
            i = randf()*ngene;
            j = randf()*ngene;
            if( i == j) continue; /* don't bother if i'm me */
            /* only recombine orientations vs translations */
            for( ibest =3; ibest< 7; ibest++)
            {
                x = (*genet)[i*7+ibest];
                (*genet)[i*7+ibest] = (*genet)[j*7+ibest];
                (*genet)[j*7+ibest] = x;
            }
        }
        /* copy the best over - this aids the convergence */
        ibest = simplex_get_best(fvals,ngene-1 );
        for( i=0; i< 7; i++)
            (*genet)[ibest*7+i] = (*gene)[ibest*7+i];
        /* lazy,lazy,lazy programming - copy the new genome over */
        for( i=0; i< 7*ngene; i++)
            (*gene)[i] = (*genet)[i];

        for( i=0; i< ngene; i++)
        {
            /* how to do the potential */
            quarternion_rot_tran( &(*gene)[i*7],ap,natom);
            (*fvals)[i] = 0.;
            for( j=0; j < inpotent; j++)
            {
                if( (*potent[j]) != v_nonbon && (*potent[j]) != u_v_nonbon)
                    (*potent[j])(&(*fvals)[i],0.);
                else
                    zone_nonbon(&(*fvals)[i],0., ap, natom);
            }
            if( (*fvals)[i] < fbest){
                fbest = (*fvals)[i];
                for( j=0; j< 7; j++)
                    best[j] = (*gene)[7*i+j];
            }
        }

    }/*k */

DONE:;
    printf(" putting best %f into coordinates \n",fbest);
    fflush(stdout);
    quarternion_rot_tran( best,ap,natom);
    free( ap );
    free( gene );
    free( fitness );
    free( fvals );
    free( studbook);
}/* end of the routine */

