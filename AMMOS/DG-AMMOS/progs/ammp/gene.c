/* gene.c
*
* genetic optimizer for AMMP
*
* given potentials use the genetic algorithm to find an optimum
*
*/
/*
*  copyright 1993,1994 Robert W. Harrison
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

void gene(op, vfs,ffs,nfs, nstep, ndeep, sigma, target,nos)
int nfs,(*vfs[])(),(*ffs[])();
int nstep,ndeep,nos;
float sigma,target;
FILE *op;
{
    float (*xvs)[];
    float (*yvs)[];
    float (*zvs)[];
    float (*value)[];
    float (*cumvalue)[];
    float vmax,vmin;
    float randg(),randf();
    float sigo,mean,x,y;
    int thebest,theworst;
    int numatm,a_number();
    ATOM *ap,*a_next();
    int i,j,k,l;
    void gene_valid_coordinate();

    /* initialize */
    numatm = a_number();
    if( numatm < 1) {
        aaerror(" no atoms defined ?\n"); return ;
    }
    xvs = malloc( ndeep*numatm*sizeof(float));
    yvs = malloc( ndeep*numatm*sizeof(float));
    zvs = malloc( ndeep*numatm*sizeof(float));
    value = malloc( ndeep*sizeof(float));
    cumvalue = malloc( ndeep*sizeof(float));
    if( xvs == NULL )
    { aaerror(" cannot allocate memory in gene\n"); exit(0);}
    if( yvs == NULL )
    { aaerror(" cannot allocate memory in gene\n"); exit(0);}
    if( zvs == NULL )
    { aaerror(" cannot allocate memory in gene\n"); exit(0);}
    if( value == NULL )
    { aaerror(" cannot allocate memory in gene\n"); exit(0);}
    if( cumvalue == NULL )
    { aaerror(" cannot allocate memory in gene\n"); exit(0);}

    for( i=0;i< numatm; i++)
    { ap = a_next(i);
        k = 0;
        for( j=0; j< ndeep; j++)
        {
            (*xvs)[k+i] = ap->x ;
            if( ap->active) (*xvs)[k+i] += randg()*sigma;
            (*yvs)[k+i] = ap->y ;
            if( ap->active) (*yvs)[k+i] += randg()*sigma;
            (*zvs)[k+i] = ap->z ;
            if( ap->active) (*zvs)[k+i] += randg()*sigma;
            k += numatm;
        }/* end loop j */
    }/* end of loop i */

    k = 0;
    for( j=0; j< ndeep; j++)
    {
        (*value)[j] = 0;
        for( i=0;i< numatm; i++)
        { ap = a_next(i);
            ap->x = (*xvs)[k+i];
            ap->y = (*yvs)[k+i];
            ap->z = (*zvs)[k+i];
        }
        k+= numatm;
        gene_valid_coordinate();
        for( i=0; i< nfs; i++)
        { (*vfs[i])( &(*value)[j],0.);}
        fprintf(op," Genetic initialization %d to %f\n", j, (*value)[j]);
    }
    /* do the work */

    for( l=0; l< nstep; l++)
    {
        sigo= 0.; mean = 0.;
        vmin = 10e10;
        thebest = -1;
        vmax = -10e10;
        theworst = -1;
        for( j=0; j< ndeep; j++)
        {
            if( (*value)[j] < vmin ) {
                vmin = (*value)[j];
                thebest = j;
            }
            if( (*value)[j] > vmax ) {
                vmax = (*value)[j];
                theworst = j;
            }
            sigo += (*value)[j]*(*value)[j];
            mean += (*value)[j];
        }
        mean = mean/ndeep;
        sigo = sqrt( sigo/ndeep - mean*mean);
        if( sigo < target) goto DONE ;
        x = 0;
        for( j=0; j< ndeep; j++)
        {
            y = vmax -(*value)[j] ;
            /*
            		if( y > 4*sigo ) y = 4*sigo;
            		if( y < -4*sigo ) y = -4*sigo;
            */
            if( y < vmax - mean - 2*sigo) y = vmax - mean - 2*sigo;
            x += y ;
            (*cumvalue)[j] = x;
        }
        if( x <= 1.e-4) goto DONE;
        for( j=0; j< ndeep; j++)
        {
            (*cumvalue)[j] /= x;
        }
        mean = randf();
        for( j=1; j< ndeep ; j++)
        {
            if( (*cumvalue)[j-1] < mean && (*cumvalue)[j] > mean)
            { thebest = j; break; }
        }
        /* make up the coords */
        k = thebest*numatm;
        if( nos > 0 && l < ndeep) k = l*numatm;
        for( i=0; i< numatm; i++)
        { ap = a_next(i);
            mean = randf();
            if( mean < 2./(float)numatm){
                mean = randf();
                if( (*cumvalue)[0] > mean)
                {
                    k = 0;
                }else{
                    for( j=1; j< ndeep ; j++)
                    {
                        if( (*cumvalue)[j-1] < mean && (*cumvalue)[j] > mean)
                        { k = j*numatm; break; }
                    }
                }
                if( nos > 0 && l < ndeep) k = l*numatm;
            }
            /* make up the coords */
            if( ap->active){
                ap->x = (*xvs)[ k + i]/* + .4*randf()-.2 */;
                ap->y = (*yvs)[ k + i]/* + .4*randf()-.2*/;
                ap->z = (*zvs)[ k + i]/* + .4*randf()-.2*/;
            }else{
                ap->x = (*xvs)[ k + i];
                ap->y = (*yvs)[ k + i];
                ap->z = (*zvs)[ k + i];
            }
        }
        /* minimize a little */
        /*
        steep( vfs,ffs, nfs, nos, 0.);
        */	
        gene_valid_coordinate();
        cngdel( vfs,ffs, nfs, nos,nos, 0.,0);
        /* insert into the queue */
        x = (*value)[0];
        thebest = 0;
        for( i=1; i< ndeep; i++)
        { if ( (*value)[i] < x) {
                x = (*value)[i]; thebest = i;
            }
        }
        x = 0.;
        for( i=0; i< nfs; i++)
        { (*vfs[i])( &x,0.);}
        fprintf( op," genetic> %f %f %f\n",
                 (*value)[theworst],x,(*value)[thebest]);
        if( x < vmax ){
            k = theworst*numatm;
            (*value)[theworst] = x;
            for( i=0; i< numatm; i++)
            {ap = a_next(i);
                (*xvs)[k+i] = ap->x;
                (*yvs)[k+i] = ap->y;
                (*zvs)[k+i] = ap->z; }
        }
    }

    /* clean up */
DONE:
    vmin = 10e10;
    thebest = -1;
    for( j=0; j< ndeep; j++)
    {  if( (*value)[j] < vmin ) {
            vmin = (*value)[j];
            thebest = j;
        }
    }
    thebest *= numatm;
    for( i=0; i< numatm; i++)
    { ap = a_next(i);
        ap->x = (*xvs)[thebest + i];
        ap->y = (*yvs)[thebest + i];
        ap->z = (*zvs)[thebest + i];
    }
    free( cumvalue);
    free( value);
    free( zvs );
    free( yvs);
    free( xvs);


}/* end of routine */


void gene_valid_coordinate()
{

    int na,a_number();
    ATOM *ap1,*ap2,*a_next();
    int i,j;
    float x,y,z;

    na = a_number();
    if( na < 1) return ;

    ap1 = a_next(-1);
    ap1 = ap1->next;
    for( i=1; i< na; i++)
    {
        for( j=0; j< i; j++)
        {
            ap2 = a_next(j);

            x = fabs(ap1->x -ap2->x);
            if( x < 1.e-5)
            {
                y = fabs(ap1->y -ap2->y);
                if( y < 1.e-5)
                {
                    z = fabs(ap1->z -ap2->z);
                    if( z< 1.e-5)
                    {
                        ap2->x += 1.e-4;
                        ap2->y += 1.e-4;
                        ap2->z += 1.e-4;
                    }}}
        }
        ap1 = ap1->next;
    }

}
