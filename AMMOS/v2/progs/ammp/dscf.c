/* direct scf calculation routines
*
*/
/*
*  copyright 1993,1994,1995 Robert W. Harrison
*  
*  This notice may not be removed
*  This program may be copied for scientific use
*  It may not be sold for profit without explicit
*  permission of the author(s) who retain any
*  commercial rights including the right to modify 
*  this notice
*/

#include <assert.h>

#define ANSI 1
/* misc includes - ANSI and some are just to be safe */
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#ifdef ANSI
#include <stdlib.h>
#endif
#include "ammp.h"
#include "orbit.h"

#define MAXPARAM 2000 

static float (*dscf_overlap_inverse)[] = NULL;

int direct_scf( output, nstep ,toler, what)
FILE *output;
int nstep;
float toler;
char *what;
{
    int no,o_number();
    ORBITAL *op, *op1 ,*o_next(), *o_m_serial();
    ATOM *ap;
    int iter,i,j;
    char whattotry;
    float ecurrent,enew;
    float initstep;
    /* space to krylov up to MAXPARAM params at once
    *  this is slightly too many 
    */
    int inparam;
    static int use_brute = (1==1);
    float *addresses[MAXPARAM];
    ORBITAL *opps[MAXPARAM];

    int dscf_krymin(); /* gradient search routine */
    int dscf_brute(); /*  brute force routine */
    float Ritz_product(),Loose_Ritz_product();
    float INDO_Ritz_product();
    float tight_INDO_Ritz_product();
    float Ritz_orth();
    float (*touse)();
    void (*forgrad)();
    void Ritz_orth_grad();
    void Ritz_orth_grad_virial();
    void dscf_make_map();

    no = o_number();
    if( no <= 0 ) return ;

    whattotry = 0; /* just do coefficients */
    if( what != NULL && *what != '\0' )
    {
        if ( strncmp( what,"coef",4)  == 0 ) whattotry = 0;  /* coefficients */
        if ( strncmp( what,"expo",4)  == 0 ) whattotry = 1;  /* exponents */
        if ( strncmp( what,"xyz",3)   == 0 ) whattotry = 2; /* atomic centers */
        if ( strncmp( what,"geom",4)  == 0 ) whattotry = 3; /*orbital geometry*/
        if ( strncmp( what,"ana",3)  == 0 ) whattotry = 4; /* just report */
        if ( strncmp( what,"fre",3)  == 0 ) whattotry = 5; /* freeze */
        if ( strncmp( what,"tha",3)  == 0 ) whattotry = 6; /* thaw */
        if ( strncmp( what,"gan",3)  == 0 ) whattotry = 7; /* gang together */
        if ( strncmp( what,"poli",4) == 0 ) whattotry = 8;
        if ( strncmp( what,"edmap",3) == 0 ) whattotry = 9;
        if ( strncmp( what,"phimap",3) == 0 ) whattotry = 10;
        if ( strncmp( what,"empir",3) == 0 ) whattotry = 11;
        if ( strncmp( what,"charge",3) == 0 ) whattotry = 12;
        if ( strncmp( what,"indo",3) == 0 ) whattotry = 13;
        if ( strncmp( what,"igeom",3) == 0 ) whattotry = 14;
        if ( strncmp( what,"ipoli",3) == 0 ) whattotry = 15;
        if ( strncmp( what,"couple",3) == 0 ) whattotry = 16;
        if ( strncmp( what,"overlap",3) == 0 ) whattotry = 17;
        if ( strncmp( what,"brute",3) == 0 ) whattotry = 18;
        if ( strncmp( what,"gradient",3) == 0 ) whattotry = 19;
    }

    /* now inside of krymin
    	for( iter=0; iter< nstep; iter++)
    */
    for( iter=0; iter< 1;iter++)
    {
        op = o_next(-1);
        inparam = 0;
        for( i=0; i< no; i++)
        {

            switch( whattotry )
            {
            case 0:
                if( !op->active) goto SKIP;
                if( op->gang != NULL) goto SKIP;
                /*	if( op->ncouple < 1) goto SKIP;
                */
                if( op->ncouple >= 1){
                    for(j=1; j< op->n; j++)
                    {
                        addresses[inparam+j-1] = &(op->a[j]);
                        opps[inparam+j-1] = op;
                    }
                    inparam += op->n-1;
                } else {
                    for(j=0; j< op->n; j++)
                    {
                        addresses[inparam+j] = &(op->a[j]);
                        opps[inparam+j] = op;
                    }
                    inparam += op->n;
                }
                touse = Loose_Ritz_product;
                forgrad = Ritz_orth_grad;
                initstep = 1.;
                break;
            case 1:
                if( !op->active) goto SKIP;
                if( op->gang != NULL) goto SKIP;
                /*	if( op->ncouple < 1) goto SKIP;
                */
                for(j=0; j< op->n; j++)
                {
                    addresses[inparam+j] = &(op->rl[j]);
                    opps[inparam+j] = op;
                }
                inparam += op->n;
                for(j=1; j< op->n; j++)
                {
                    addresses[inparam+j-1] = &(op->a[j]);
                    opps[inparam+j-1] = op;
                }
                inparam += op->n-1;
                /*
                	forgrad = Ritz_orth_grad;
                	touse = Loose_Ritz_product;
                */
                forgrad = Ritz_orth_grad;
                touse = Ritz_orth;
                initstep = .1;
                break;
            case 2:
                /*	if( !op->active) goto SKIP;
                */
                ap = op->myatom;
                if( ap->active ){
                    addresses[inparam] = &(ap->x);
                    opps[inparam] = op;
                    addresses[inparam+1] = &(ap->y);
                    opps[inparam+1] = op;
                    addresses[inparam+2] = &(ap->z);
                    opps[inparam+2] = op;
                    inparam += 3;
                    /*
                    	if( op->gang == NULL){
                    	for(j=0; j< op->n; j++)
                    	{
                    		addresses[inparam+j] = &(op->rl[j]);
                    		opps[inparam+j] = op;
                    	}
                    	inparam += op->n;
                    	for(j=1; j< op->n; j++)
                    	{
                    		addresses[inparam+j-1] = &(op->a[j]);
                    		opps[inparam+j-1] = op;
                    	}
                    	inparam += op->n-1;
                    			}
                    */
                    initstep = .01;
                }
                touse = Ritz_product;
                touse = Ritz_orth;
                forgrad = Ritz_orth_grad;
                break;
            case 3:
                if( !op->active) goto SKIP;
                if( op->gang != NULL) goto SKIP;
                touse = Ritz_product;
                touse = Loose_Ritz_product;
                forgrad = Ritz_orth_grad;
                if( op->gang == NULL){
                    for(j=0; j< op->n; j++)
                    {
                        addresses[inparam+j] = &(op->rl[j]);
                        opps[inparam+j] = op;
                    }
                    inparam += op->n;
                    for(j=1; j< op->n; j++)
                    {
                        addresses[inparam+j-1] = &(op->a[j]);
                        opps[inparam+j-1] = op;
                    }
                    inparam += op->n-1;
                }
                switch( op->type )
                {
                case Or1:
                    goto SKIP;
                    break;
                case Or1o:
                    addresses[inparam] = &op->x;
                    opps[inparam] = op;
                    addresses[inparam+1] = &op->y;
                    opps[inparam+1] = op;
                    addresses[inparam+2] = &op->z;
                    opps[inparam+2] = op;
                    inparam += 3;
                    initstep = .01;
                    break;
                case Or2:
                    addresses[inparam] = &op->along;
                    opps[inparam] = op;
                    inparam += 1;
                    initstep = .01;
                    break;
                case Or2p:
                    addresses[inparam] = &op->along;
                    opps[inparam] = op;
                    inparam += 1;
                    initstep = .01;
                    break;
                case Or3:
                    addresses[inparam] = &op->along;
                    opps[inparam] = op;
                    addresses[inparam+1] = &op->along2;
                    opps[inparam+1] = op;
                    inparam += 2;
                    initstep = .01;
                    break;
                case Or4p:
                case Or4s:
                    addresses[inparam] = &op->along;
                    opps[inparam] = op;
                    inparam += 1;
                    initstep = .01;
                    break;
                case Orm:
                    goto SKIP;
                    break;

                }
                break;
            case 4:
                o_update_normals();
                renormalize();
                Ritz_product(op);
                report(output);
                return;
                break;
            case 5:
                op = o_m_serial(nstep);
                if( op == NULL ) return;
                op->active = (1== 0);
                return ;
                break;
            case 6:
                op = o_m_serial(nstep);
                if( op == NULL ) return;
                op->active = (1== 1);
                return ;
                break;
            case 7:
                op = o_m_serial(nstep);
                op1 = o_m_serial( (int) toler);
            if(op == NULL ) { return ;}
                if(op1 == NULL ) {op->gang = NULL; return ;}
                /* check for indirect gangs ang gang loops !!! */
                while(op1->gang != NULL && op1->gang != op ) op1 = op1->gang;
                if( op->type != op1->type ) return;
                op->gang = op1;
                return;
                break;
            case 8:
                if( op->gang != NULL) goto SKIP;
                touse = Ritz_product;
                forgrad = Ritz_orth_grad;
                forgrad = Ritz_orth_grad_virial;

                if( op->gang == NULL){
                    for(j=0; j< op->n; j++)
                    {
                        addresses[inparam+j] = &(op->rl[j]);
                        opps[inparam+j] = op;
                    }
                    inparam += op->n;
                    if( op->ncouple >= 1){
                        for(j=1; j< op->n; j++)
                        {
                            addresses[inparam+j-1] = &(op->a[j]);
                            opps[inparam+j-1] = op;
                        }
                        inparam += op->n-1;
                    } else {
                        for(j=0; j< op->n; j++)
                        {
                            addresses[inparam+j] = &(op->a[j]);
                            opps[inparam+j] = op;
                        }
                        inparam += op->n;

                    }
                }
                break;
            case 9:
                dscf_make_map( 0,output,toler,(float)nstep);
                return ;
                break;
            case 10:
                dscf_make_map( 1,output,toler,(float)nstep);
                return ;
                break;
            case 11:
                dscf_make_map( 2,output,toler,(float)nstep);
                return ;
                break;
            case 12:
                dscf_fit_q( (float)nstep);
                return ;
                break;
            case 13:
                if( !op->active) goto SKIP;
                if( op->gang != NULL) goto SKIP;
                for(j=0; j< op->n; j++)
                {
                    addresses[inparam+j] = &(op->rl[j]);
                    opps[inparam+j] = op;
                }
                inparam += op->n;
                for(j=1; j< op->n; j++)
                {
                    addresses[inparam+j-1] = &(op->a[j]);
                    opps[inparam+j-1] = op;
                }
                inparam += op->n-1;
                forgrad = Ritz_orth_grad;
                touse = INDO_Ritz_product;
                initstep = .1;
                break;
            case 14:
                if( !op->active) goto SKIP;
                if( op->gang != NULL) goto SKIP;
                touse = Ritz_product;
                touse = tight_INDO_Ritz_product;
                forgrad = Ritz_orth_grad;
                /* touse = Loose_Ritz_product; */

                if( op->gang == NULL){
                    for(j=0; j< op->n; j++)
                    {
                        addresses[inparam+j] = &(op->rl[j]);
                        opps[inparam+j] = op;
                    }
                    inparam += op->n;
                    for(j=1; j< op->n; j++)
                    {
                        addresses[inparam+j-1] = &(op->a[j]);
                        opps[inparam+j-1] = op;
                    }
                    inparam += op->n-1;
                }
                switch( op->type )
                {
                case Or1:
                    goto SKIP;
                    break;
                case Or1o:
                    addresses[inparam] = &op->x;
                    opps[inparam] = op;
                    addresses[inparam+1] = &op->y;
                    opps[inparam+1] = op;
                    addresses[inparam+2] = &op->z;
                    opps[inparam+2] = op;
                    inparam += 3;
                    initstep = .01;
                    break;
                case Or2:
                    addresses[inparam] = &op->along;
                    opps[inparam] = op;
                    inparam += 1;
                    initstep = .01;
                    break;
                case Or3:
                    addresses[inparam] = &op->along;
                    opps[inparam] = op;
                    addresses[inparam+1] = &op->along2;
                    opps[inparam+1] = op;
                    inparam += 2;
                    initstep = .01;
                    break;
                case Or4p:
                case Or4s:
                    addresses[inparam] = &op->along;
                    opps[inparam] = op;
                    inparam += 1;
                    initstep = .01;
                    break;
                case Orm:
                    goto SKIP;
                    break;

                }
                break;

            case 15:
                if( op->gang != NULL) goto SKIP;
                touse = Ritz_product;
                forgrad = Ritz_orth_grad;
                touse = tight_INDO_Ritz_product;

                if( op->gang == NULL){
                    for(j=0; j< op->n; j++)
                    {
                        addresses[inparam+j] = &(op->rl[j]);
                        opps[inparam+j] = op;
                    }
                    inparam += op->n;
                    for(j=1; j< op->n; j++)
                    {
                        addresses[inparam+j-1] = &(op->a[j]);
                        opps[inparam+j-1] = op;
                    }
                    inparam += op->n-1;
                }
                break;

            case 16: /* couple two orbitals */
                op = o_m_serial(nstep);
                op1 = o_m_serial( (int) toler);
            if(op == NULL || op->ncouple <1) { return ;}
                if(op1 == NULL ) { return ;}
                j = op->ncouple;
                op->couple[j] = op1;
                op->ncouple += 1;
                op1->ncouple = -1;
                op1->couple[1] = op;
                op1->couple[0] = op1;
                return;
                break;
            case 17: /* overlap analysis */
                o_update_normals();
                renormalize();
                overlap( output);
                return;
                break;
            case 18:
                use_brute = (1==1);
                return;
                break;
            case 19:
                use_brute = (0==1);
                return;
                break;

            }/* end switch( whattotry) */

SKIP:  ;

            op = op->next;

        }/* i */
        o_update_normals();
        renormalize();
        ecurrent = (*touse)();
        initstep = .01;
        if( inparam > 0 ){
            if( use_brute)
                dscf_brute( addresses,opps,inparam,initstep, touse,forgrad,output,nstep);
            else
                dscf_krymin( addresses,opps,inparam,initstep, touse,forgrad,output,nstep);
            renormalize();
        }

        enew = (*touse)();
        report( output );
        fprintf(output,"iteration %d start %f finish %f delta %e\n",
                iter, ecurrent,enew,enew-ecurrent);
        /*
        	if( enew > ecurrent) break;
        */
        if( fabs( enew-ecurrent) < toler ) break;
        ecurrent = enew;

    }/* iter */

}/* end of directscf */
/* krymin is now a finite difference engine !!! */
int dscf_brute( where,whatO,n,is,what,forgrad,output,nstep )
FILE *output;
int n,nstep;
float is;
float *where[];
ORBITAL *whatO[];
float (*what)();
void (*forgrad)();
{
    void *parameter;
    static	float search[MAXPARAM],ograd[MAXPARAM];
    static	float params[MAXPARAM],grad[MAXPARAM];
    float dscf_line_brute();
    float Just_Ritz_product();
    float standard;
    int i,j;
    int iter;
    int renormalize();
    float beta,betad;
    void Ritz_orth_grad();
    assert( n < MAXPARAM);

    if( n== 0 ) return;

    for( i=0; i<n; i++)
    {  grad[i] = 0.; search[i] = 0.; }


    for( iter=0;iter<nstep; iter++)
    {

        renormalize();
        for( i=0; i< n; i++)
        { params[i] = *where[i]; ograd[i] =grad[i];}

        for(j=0; j< n; j++)
        {
            for( i=0; i< n; i++)
            { params[i] = *where[i]; ograd[i] =grad[i];}
            for(i=0; i< n; i++)
                search[i] = 0.;
            search[j] = fabs(0.1*(*where[j]));
            if( search[j] < .01) search[j] = 0.01;
            parameter = whatO[j];
            beta = dscf_line_brute( search,where,params,n,what,parameter);
        }
        for(j=n; j> 0; j--)
        {
            for( i=0; i< n; i++)
            { params[i] = *where[i]; ograd[i] =grad[i];}
            for(i=0; i< n; i++)
                search[i] = 0.;
            search[j-1] = -fabs(0.1*(*where[j-1])) ;
            if( search[j-1] > -.01) search[j-1] = -0.01;
            parameter = whatO[j-1];
            beta = dscf_line_brute( search,where,params,n,what,parameter);
        }
        report( output );

    } /* iter */
}

/* krymin is now a finite difference engine !!! */
int dscf_krymin( where,whatO,n,is,what,forgrad,output,nstep )
FILE *output;
int n,nstep;
float is;
float *where[];
ORBITAL *whatO[];
float (*what)();
void (*forgrad)();
{
    void *parameter;
    static	float search[MAXPARAM],ograd[MAXPARAM];
    static	float params[MAXPARAM],grad[MAXPARAM];
    float dscf_line();
    float Just_Ritz_product();
    float standard;
    int i,j;
    int iter;
    int renormalize();
    float beta,betad;
    void Ritz_orth_grad();
    assert( n < MAXPARAM);

    if( n== 0 ) return;

    for( i=0; i<n; i++)
    {  grad[i] = 0.; search[i] = 0.; }


    for( iter=0;iter<nstep; iter++)
    {

        renormalize();
        for( i=0; i< n; i++)
        { params[i] = *where[i]; ograd[i] =grad[i];}
        printf("gradient started");
        fflush(stdout);
        (*forgrad)(grad,where,params,whatO,n);
        beta = 0.; betad = 0.;
        for(i=0; i<n; i++)
        {
            beta += (grad[i] - ograd[i])*grad[i];
            betad += ograd[i]*ograd[i];
        }
        if( fabs(betad) > 0.1* beta )
        { beta = beta/betad; }else { beta = 0. ;}
        for( i=0;i< n;i++)
        {
            search[i] = (grad[i] + beta*search[i]);
        }

        for( i=0; i< n; i++)
        {  *where[i] = params[i];}

        printf(" and done\n");
        /* don't care about this return value */
        for( i=0; i< n; i++)
        { params[i] = *where[i]; ograd[i] =grad[i];}

        beta = dscf_line( search,where,params,n,what,parameter);
        report( output );

    } /* iter */
}

/* dscf_line is called by dscf_krymin()
*
* it searches for the closest local minimum
*
*/
float dscf_line_brute( search,where,params,n, what ,whom)
int n;
float search[],params[],*where[];
float (*what)();
void *whom;
{
    float beststep,step,dstep,myval,bestval;
#define MAXBUF 100
    float steps[MAXBUF],vals[MAXBUF];
    /* standard non-cached calls */
    float Loose_Ritz_product(),Ritz_orth(),Ritz_product(),
    INDO_Ritz_product(),tight_INDO_Ritz_product();
    /* Cached version(s) */
    float Cached_Ritz_orth();
    int inbuf;
    int nfunc;
    int i,irun,itry;

    if( what == Ritz_orth || what == Loose_Ritz_product )
        what = Cached_Ritz_orth;

    inbuf = 1;
    steps[0] = 0.;
    if( what == Cached_Ritz_orth){
        vals[0] = (*what)(whom,-1); /* initialize the cache */
    } else {
        vals[0] = (*what)(whom);
    }
    bestval = vals[0];
    nfunc = 1;

    step = 0.;
    beststep = step;
    dstep = 1.;
    dstep = .1;
    for( itry = 0; itry < 7; itry ++)
    {
        /* use irun to search  in terms of dstep */
        for( irun = 0; irun < 100; irun ++)
        {
            step += dstep;
            /* first look in the table */

            for( i=0; i< inbuf; i++)
            {
                if( step == steps[i] ) {
                    myval = vals[i];
                    goto STEP_FOUND;
                }
            }
            /* if here then we have to calculate the function */
            for( i=0; i< n; i++)
            {
                *where[i] = params[i] + search[i]*step;
            }
            if( what == Cached_Ritz_orth){
                myval = (*what)(whom,0);
            }else{
                myval = (*what)(whom);
            }
            nfunc ++;
            if (inbuf == MAXBUF)
            {
                if( what == Cached_Ritz_orth){ (*what)(whom,1); }
                return beststep;
            }
            assert( inbuf < MAXBUF );
            if( inbuf < MAXBUF)
            {vals[inbuf] = myval;steps[inbuf++] = step;}
            else { inbuf -= 10;   }
STEP_FOUND:
            if( myval >= bestval)
            {
                step -= dstep;
                dstep = -dstep*.1;
                break;
            } else {
                bestval = myval;
                beststep = step;
                /*
                step += dstep;
                */
            }

        }/* irun */
    }/* itry */

    /* update to the best found */
    for( i=0; i< n; i++)
    {
        *where[i] = params[i] + search[i]*beststep;
    }
    printf("%d function calls \n",nfunc);
    if( what == Cached_Ritz_orth){ (*what)(whom,1); }
    return beststep;
}
/* dscf_line is called by dscf_krymin()
*
* it searches for the closest local minimum
*
*/
float dscf_line( search,where,params,n, what ,whom)
int n;
float search[],params[],*where[];
float (*what)();
void *whom;
{
    float beststep,step,dstep,myval,bestval;
#define MAXBUF 100
    float steps[MAXBUF],vals[MAXBUF];
    int inbuf;
    int nfunc;
    int i,irun,itry;

    inbuf = 1;
    steps[0] = 0.;
    vals[0] = (*what)(whom);
    bestval = vals[0];
    nfunc = 1;

    step = 0.;
    beststep = step;
    dstep = 1.;
    dstep = .1;
    for( itry = 0; itry < 10; itry ++)
    {
        /* use irun to search  in terms of dstep */
        for( irun = 0; irun < 10; irun ++)
        {
            step += dstep;
            /* first look in the table */

            for( i=0; i< inbuf; i++)
            {
                if( step == steps[i] ) {
                    myval = vals[i];
                    goto STEP_FOUND;
                }
            }
            /* if here then we have to calculate the function */
            for( i=0; i< n; i++)
            {
                *where[i] = params[i] + search[i]*step;
            }
            myval = (*what)(whom);
            nfunc ++;
            if (inbuf == MAXBUF) return beststep;
            assert( inbuf < MAXBUF );
            if( inbuf < MAXBUF)
            {vals[inbuf] = myval;steps[inbuf++] = step;}
            else { inbuf -= 10;   }
STEP_FOUND:
            if( myval >= bestval)
            {
                step -= dstep;
                dstep = -dstep*.5;
                break;
            } else {
                bestval = myval;
                beststep = step;
                /*
                step += dstep;
                */
            }

        }/* irun */
    }/* itry */

    /* update to the best found */
    for( i=0; i< n; i++)
    {
        *where[i] = params[i] + search[i]*beststep;
    }
    printf("%d function calls \n",nfunc);
    return beststep;
}

/* map.c
*
* routines to use floating gaussian orbitals for QM calculations
*
*/
/*
*  copyright 1993,1994,1995 Robert W. Harrison
*  
*  This notice may not be removed
*  This program may be copied for scientific use
*  It may not be sold for profit without explicit
*  permission of the author(s) who retain any
*  commercial rights including the right to modify 
*  this notice
*/


/* do the map */

void dscf_make_map( type,where,grid,guard)
int type;
FILE *where;
float grid,guard ; /* the grid (i.e., .5A) and guard radii */
{
    ATOM *ap,*a_next();
    int i,j,k,na,a_number();
    int ix,imin,imax;
    int jx,jmin,jmax;
    int kx,kmin,kmax;
    int irow;
    float xa,ya,za,xmax,xmin,ymax,ymin,zmax,zmin;
    void dscf_map_setup();
    float dscf_map_ED(),dscf_map_field(),dscf_map_empirical();

    na = a_number();
    if( na <= 0){ fprintf(where,"No atoms - No map\n"); return ;}

    xmax = -10e10;
    ymax = -10e10;
    zmax = -10e10;
    xmin =  10e10;
    ymin =  10e10;
    zmin =  10e10;
    for( i=0; i< na; i++)
    {
        ap = a_next(i);
        if( ap->na > 0.) {
            if( xmax < ap->x) xmax = ap->x;
            if( ymax < ap->y) ymax = ap->y;
            if( zmax < ap->z) zmax = ap->z;
            if( xmin > ap->x) xmin = ap->x;
            if( ymin > ap->y) ymin = ap->y;
            if( zmin > ap->z) zmin = ap->z;
        }
    }
    /* add guard box */
    xmax += guard;
    ymax += guard;
    zmax += guard;
    xmin -= guard;
    ymin -= guard;
    zmin -= guard;
    /* force to integer grids */
    i = xmax/grid + 1; imax = i; xmax = i*grid;
    i = ymax/grid + 1; jmax = i; ymax = i*grid;
    i = zmax/grid + 1; kmax = i; zmax = i*grid;
    i = xmin/grid - 1; imin = i; xmin = i*grid;
    i = ymin/grid - 1; jmin = i; ymin = i*grid;
    i = zmin/grid - 1; kmin = i; zmin = i*grid;
    /* figure out the fake cell (we'll keep it 90 90 90 ) */
    ix = imax -imin; jx = jmax -jmin; kx = kmax -kmin;
    /* write the header */
    fprintf(where,"\n       2 !NTITLE\nREMARKS AMMP dscf map                                                                                  \nREMARKS                                                                                                 \n");
    fprintf(where,"%8i%8i%8i%8i%8i%8i%8i%8i%8i\n",
            ix,imin,imax,jx,jmin,jmax,kx,kmin,kmax);
    fprintf(where,"%12.5e%12.5e%12.5e%12.5e%12.5e%12.5e\nZXY",
            xmax-xmin,ymax-ymin,zmax-zmin,90.,90.,90.);
    /* now do the work we loop with z outer then x then y
    * and fake a FORTRAN implied do */
    dscf_map_setup(1==1);
    for( k=0; k<= kx; k++)
    {
        fprintf(where,"\n%8d\n",k);
        irow = 0;
        za = k*grid + zmin;
        for( j=0; j<= jx; j++)
            for( i=0; i<= ix; i++)
            {
                xa = i*grid + xmin;
                ya = j*grid + ymin;
                if(irow ==6) {fprintf(where,"\n"); irow = 0;}
                if( type == 0 ) fprintf(where,"%12.5e",dscf_map_ED(xa,ya,za));
                if( type == 1 ) fprintf(where,"%12.5e",dscf_map_field(xa,ya,za));
                if( type == 2 ) fprintf(where,"%12.5e",dscf_map_empirical(xa,ya,za));
                irow ++;
            }}
    fprintf(where,"\n%d\n",-9999);
    dscf_map_setup(0==1);
}/* end of routine */


/*
* prepare normalizers so that the electrostatic potential can be 
* calculated
*/
void dscf_map_setup(flag)
int flag; /* flag true sets up, false cleans up */
{
    int o_update_normals(),no,o_number();
    int renormalize();
    void orth_invert();
    float (*ol)[],(*ok)[];
    /* simply update the orbital coefficients and normalize
    * so that phiphi(op,op) = 1. 
    */
    no = o_number();
    if( no < 1 ) return;
    if( flag){
        if( dscf_overlap_inverse != NULL){ free(dscf_overlap_inverse);
            dscf_overlap_inverse = NULL;}
        o_update_normals();
        renormalize();
        dscf_overlap_inverse = malloc(no*no*sizeof(float));
        ol = malloc(no*no*sizeof(float));
        ok = malloc(no*no*sizeof(float));
        orth_invert( dscf_overlap_inverse,ol,ok,no);
        free( ok);
        free( ol);
    }else{
        free(dscf_overlap_inverse);
        dscf_overlap_inverse = NULL;
    }

}
/* generate the electron density as a function of x,y,z */
float dscf_map_ED( x,y,z)
float x,y,z; /* where to map */
{
    float ED;
    ORBITAL *op,*oop,*op1,*oop1,*o_next();
    int i,ii,j,k,l,ll,m,n,no,o_number();
    float dx,dy,dz,r,rs,rab,re,rx,accum;
    float temp;
    float s1,s2,sign;

    ED = 0.;
    no = o_number();
    if( no <= 0 ) return ED;
    oop = o_next(0);
    for( i=0; i< no; i++)
    {
        for( ii=0; ii< oop->ncouple; ii++)
        {
            op = oop->couple[ii];
            for( l = 0; l < no; l++)
            {
                oop1 = o_next(l);
                for( ll=0; ll< oop1->ncouple ; ll++)
                {
                    op1 = oop1->couple[ll];

                    accum = 0.;
                    for( j=0; j< op->n; j++)
                    {
                        for( n=0; n< op1->n; n++)
                        {
                            rs = (op->r[j]+ op1->r[n]);
                            r = op->r[j]*op1->r[n]/rs;
                            /*
                            	temp = op->a[j]*op1->a[n] *pow( PI/rs,1.5);
                            */
                            temp = op->a[j]*op1->a[n];

                            s1 = 1;
                            for( k=0; k< op->ncenter; k++)
                            {
                                if( op->type == Or4p) s1 = -s1;
                                s2 = 1;
                                for( m=0; m< op1->ncenter; m++)
                                {
                                    if( op1->type == Or4p) s2 = -s2;
                                    sign = s1*s2;
                                    dx = (op->rx[k] - op1->rx[m]);
                                    dy = (op->ry[k] - op1->ry[m]);
                                    dz = (op->rz[k] - op1->rz[m]);
                                    rab = dx*dx + dy*dy + dz*dz;
                                    rab = rab*INVBOHR*INVBOHR;
                                    dx = (op->r[j]*op->rx[k] + op1->r[n]*op1->rx[m])/rs;
                                    dy = (op->r[j]*op->ry[k] + op1->r[n]*op1->ry[m])/rs;
                                    dz = (op->r[j]*op->rz[k] + op1->r[n]*op1->rz[m])/rs;
                                    dx -= x;
                                    dy -= y;
                                    dz -= z;
                                    rx = (dx*dx + dy*dy + dz*dz)*INVBOHR*INVBOHR;

                                    dx = rs*rx +r*rab;
                                    if( dx > 50.) dx = 50.;
                                    if( dx < -50.) dx = -50.;
                                    accum += sign*temp*exp(-dx)*(*dscf_overlap_inverse)[i*no+l];
                                } /*m */
                            } /* k */
                        } /* n */
                    } /*j */
                    if( op->ipair == 2) ED += accum;
                    ED += accum;
                } /* ll */
            } /* l */
        } /*ii */
        oop = oop->next;
    } /*i */

    return ED;
}

float dscf_map_field( x,y,z)
float x,y,z;
{
    float phi;
    ATOM *ap,*a_next();
    ORBITAL *op,*oop,*op1,*oop1,*o_next();
    int i,ii,j,k,l,ll,m,n,no,o_number(),na,a_number();
    float dx,dy,dz,rab,r,ntemp,rs,re;
    float Fzero(),phiphi();
    float s1,s2,sign;

    phi = 0.;
    /* do the nuclear terms */
    na = a_number();
    if( na <= 0 ) return phi;
    for( i=0; i< na; i++)
    {
        ap = a_next(i);
        if( ap->na != 0)
        {
            dx = ap->x - x;
            dy = ap->y - y;
            dz = ap->z - z;
            r = sqrt(dx*dx + dy*dy + dz*dz)*INVBOHR;
            if( r < 1.e-1) r = 1.e-1; /* keep the map in bounds */
            phi += ap->na/r;
        }}
    /* now the electrons */
    no = o_number();
    /* return the potential == (by def) -(integral)E */
    if( no <=0 ) return -phi;
    oop = o_next(0);
    for( i=0; i< no; i++)
    {
        for( ii=0; ii< oop->ncouple; ii++)
        { op = oop->couple[ii];

            oop1 = o_next(-1);
            for( l=0; l < no; l++)
            {
                for( ll=0; ll< oop1->ncouple; ll++)
                {
                    op1 = oop1->couple[ll];
                    for( j=0; j< op->n; j++)
                    {
                        for( m=0; m< op1->n; m++)
                        {
                            rs = (op->r[j] + op1->r[m]);
                            r = op->r[j]*op1->r[m]/rs;
                            ntemp = op->a[j]*op1->a[m]*TWOPI/rs
                                    *(*dscf_overlap_inverse)[l*no+i];

                            s1 = 1.;
                            for( k=0; k< op->ncenter; k++)
                            {
                                if( op->type == Or4p) s1 = -s1;
                                s2 = 1.;
                                for( n=0; n< op1->ncenter; n++)
                                {
                                    if( op1->type == Or4p) s2 = -s2;
                                    sign = s1*s2;
                                    dx = (op->rx[k] - op1->rx[n]);
                                    dy = (op->ry[k] - op1->ry[n]);
                                    dz = (op->rz[k] - op1->rz[n]);
                                    rab = (dx*dx + dy*dy + dz*dz)*INVBOHR*INVBOHR;
                                    re = exp(-r*rab);
                                    dx = (op->r[j]*op->rx[k] + op1->r[m]*op1->rx[n])/rs;
                                    dy = (op->r[j]*op->ry[k] + op1->r[m]*op1->ry[n])/rs;
                                    dz = (op->r[j]*op->rz[k] + op1->r[m]*op1->rz[n])/rs;
                                    dx -= x;
                                    dy -= y;
                                    dz -= z;
                                    rab = (dx*dx + dy*dy + dz*dz)*INVBOHR*INVBOHR;
                                    if( rab < 1.e-4) rab = 1.e-4;


                                    /*
                                    	only if rab is not the sqrt
                                    	in phiH4phi ... rs*rab is the square of here   
                                    		phi -= ntemp*Fzero( rs*rab)*re;
                                    		phi -= sign*2*ntemp * 0.5 *ROOTPI/(rab)*erf(rab)*re;
                                    */
                                    phi -= sign*2*ntemp*Fzero( rs*rab)*re;
                                    /* the above sign is + this gives the same values as others
                                    * but is not neccesarily correct */
                                }
                            }
                        }
                    }} oop1 = oop1->next; } /*l */
        }
        oop = oop->next;
    }
    /* return the potential == (by def) -(integral)E */
    return -phi;
}
float dscf_map_empirical( x,y,z)
float x,y,z;
{
    float phi;
    ATOM *ap,*a_next();
    ORBITAL *op,*op1,*o_next();
    int i,j,k,l,m,n,no,o_number(),na,a_number();
    float dx,dy,dz,rab,r,ntemp,rs,re;
    float Fzero(),phiphi();
    float s1,s2,sign;

    phi = 0.;
    /* do the nuclear terms only use the fit charges */
    na = a_number();
    if( na <= 0 ) return phi;
    for( i=0; i< na; i++)
    {
        ap = a_next(i);
        if( ap->na != 0)
        {
            dx = ap->x - x;
            dy = ap->y - y;
            dz = ap->z - z;
            r = sqrt(dx*dx + dy*dy + dz*dz)*INVBOHR;
            if( r < 1.e-2) r = 1.e-2; /* keep the map in bounds */
            phi += ap->q/r;

        }}

    return phi;
}


#define Kollman
/* fit the charges total charge = total */
int dscf_fit_q(total)
float total;
{
    int na,a_number(),no,o_number();
    int i,j,k;
    int ncharge;
    float (*matrix)[],(*vector)[],(*radii)[];
    float x,y,z,xc,yc,zc,r,t,randf();
    float xmax,ymax,zmax,xmin,ymin,zmin;
    int ix,iy,iz;
    void rand3();
    int mom_solve();
    float dscf_map_field();
    float dscf_map_empirical();
    void dscf_map_setup();
    ATOM *ap,*a_next(),*(*todo)[];
    ORBITAL *op,*o_next();
    float weight;

    na = a_number();
    if( na == 0 ) return 0 ;
    /* first figure out how many atoms and get their addresses */
    ncharge = 0;
    xc = 0.;
    yc = 0.;
    zc = 0.;
    xmax = -10e10;
    ymax = -10e10;
    zmax = -10e10;
    xmin = 10e10;
    ymin = 10e10;
    zmin = 10e10;
    for( i=0; i< na; i++)
    { ap = a_next(i);
        if( ap->na >= 0.1 ){ ncharge += 1;
            xc += ap->x;
            yc += ap->y;
            zc += ap->z;
            if( ap->x > xmax ) xmax = ap->x;
            if( ap->y > ymax ) ymax = ap->y;
            if( ap->z > zmax ) zmax = ap->z;
            if( ap->x < xmin ) xmin = ap->x;
            if( ap->y < ymin ) ymin = ap->y;
            if( ap->z < zmin ) zmin = ap->z;
        }
    }/* i */
    xc /= ncharge;
    yc /= ncharge;
    zc /= ncharge;
    todo = malloc( ncharge*sizeof(ATOM *) );
    if( todo == NULL ) { aaerror("cannot allocate memory in dscf_fit_q");
        return 1; }
    matrix = malloc( (ncharge+1)*(ncharge+1)*sizeof(float) );
    if( matrix == NULL ) { aaerror("cannot allocate memory in dscf_fit_q");
        return 1; }
    vector = malloc( (ncharge+1)*sizeof(float) );
    if( vector == NULL ) { aaerror("cannot allocate memory in dscf_fit_q");
        return 1; }
    radii = malloc( (ncharge)*sizeof(float) );
    if( radii == NULL ) { aaerror("cannot allocate memory in dscf_fit_q");
        return 1; }
    j = 0;
    for( i=0; i< na; i++)
    { ap = a_next(i);
        if( ap->na >= 0.1) (* todo)[j++] = ap ;
    }/* i */
    /* setup storage with constraint equations */
    for( i=0; i< ncharge; i++)
    {

        for( j=0; j< ncharge; j++)
            (*matrix)[i*ncharge + i  + j ] = 0.;

#ifdef Kollman
        (*matrix)[i*(ncharge+1 ) + ncharge] = 1.;
#else
(*matrix)[i*(ncharge+1 ) + ncharge] = 1.;
#endif
        (*matrix)[ncharge*(ncharge+1)  + i] = 1.;
        (*vector)[i] = 0.;
    }
    (*vector)[ncharge] = total;
#ifdef Kollman
    (*matrix)[(ncharge+1)*(ncharge+1) -1 ] = 0.;
#else
    (*matrix)[(ncharge+1)*(ncharge+1) -1 ] = 1.;
#endif


    dscf_map_setup(1==1);
    /*
    	for( k=0; k< ncharge*20; k++)
    	{
    */
    /* pick a point */
    /*
    REPICK:	rand3(&x,&y,&z);
    	r = randf()*20.; 	
    	x = x*r + xc;
    	y = y*r + yc;
    	z = z*r + zc;
    */
    for( iz = (int)zmin-8; iz< 8+(int)zmax;iz++)
        for( iy = (int)ymin-8; iy< 8+(int)ymax;iy++)
            for( ix = (int)xmin-8; ix< 8+(int)xmax;ix++)
            {
                x = (float)ix + xc;
                y = (float)iy + yc;
                z = (float)iz + zc;
                for( i=0; i< ncharge; i++)
                { ap = (*todo)[i];
                    if( ap == NULL ) aaerror("ap is NULL");;
                    t = x-ap->x;
                    r = t*t;
                    t = y-ap->y;
                    r += t*t;
                    t = z-ap->z;
                    r += t*t;
                    /*
                    		if( r < 4  ) goto SKIP;
                    		if( r < 3.0 ) goto REPICK;
                    		if( r < 2.5 ) goto SKIP;
                    		if( r < 1  ) goto SKIP;
                    */
                    if( r < 9.  ) goto SKIP;
                    /*
                    		(*radii)[i] = one/sqrt(r)*INVBOHR;
                    		(*radii)[i] = one/sqrt(r);
                    */
                    (*radii)[i] = one/sqrt(r);
                } /* i */
                /* if here we have a good point */
                t = dscf_map_field(x,y,z)*BOHR;

                for( i=0;i< ncharge; i++)
                {
                    for( j=0; j< ncharge; j++)
                        (*matrix)[i*ncharge + i  + j] += (*radii)[i]*(*radii)[j];
                    (*vector)[i] += (*radii)[i]*t;
                }
SKIP: ;
            }/*k */
    /*
    	for( i= 0; i< ncharge+1; i++){
    	for( j=0; j< ncharge+1; j++)
    		printf("%f ",(*matrix)[i*ncharge + i + j]);
    	printf("\n");}
    */

    mom_solve( matrix,vector,ncharge+1,ncharge+1);

    for( i=0; i< ncharge; i++)
    { ap = (*todo)[i];
#ifdef Kollman
        ap->q = (*vector)[i] ;
#else
        ap->q = (*vector)[i] + (*vector)[ncharge]/ncharge ;
#endif
    }
    /* testing */
    r = 0.;
    t = 0.;
    for( i= 0; i< 100; i++)
    {
        rand3(&x,&y,&z);
        x = 2*(x-.5) + xc;
        y = 2*(y-.5) + yc;
        z = 2*(z-.5) + zc;
        (*vector)[0] = dscf_map_field(x,y,z)*BOHR;
        (*vector)[1] = dscf_map_empirical(x,y,z);
        r = r + fabs((*vector)[0]-(*vector)[1]);
        t = t + fabs((*vector)[0])+fabs((*vector)[1]);

        printf("%f %f %f %f %f\n",(*vector)[0],(*vector)[1],
               (*vector)[0]-(*vector)[1],r,t);
    }

    printf("%f\n",r/t);
    dscf_map_setup(1==0);
    free(radii);
    free(vector);
    free(matrix);
    free(todo);
    return 0;
} /* end of dscf_fit_q() */
