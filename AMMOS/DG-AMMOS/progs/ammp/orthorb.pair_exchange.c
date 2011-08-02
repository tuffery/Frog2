/*
* orthorb.c
*
* perform the Ritz sum on orthonormalized combinations
* of the basis set.
*/

/*
*  copyright 1997 Robert W. Harrison
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
#include "orbit.h"

#define orth_orb_dstep 1.e-4

extern ORBITAL *firstORBITAL;
extern float total_ex,total_col,total_nuc,total_enuc,total_kinet,total_norm;

float Ritz_orth()
{
    float (*gscoef)[];
    float (*o_cache)[];
    float (*overlaps)[];
    float H,ritzd,ritzn;
    int no,o_number();
    int i,j,k,l,m,n;
    float cv,ev,x,y;
    float orth_kinet();
    float orth_nuclear();
    float orth_exchange();
    float orth_coloumb();
    float nuclear(),phiphi();
    void o_update_normals(),orth_invert();
    ORBITAL *o1,*o2,*o3,*o4, *o_next();

    no = o_number();
    if(no < 1) return 0.;
    i = no*no;
    gscoef = malloc( i*sizeof(float));
    if( gscoef == NULL ) {
        aaerror("cannot allocate memory in Ritz_orth");
        exit(0); }
    overlaps = malloc( i*sizeof(float));
    if( overlaps == NULL ) {
        aaerror("cannot allocate memory in Ritz_orth");
        exit(0); }
    o_cache = malloc( i*sizeof(float));
    if( o_cache == NULL ) {
        aaerror("cannot allocate memory in Ritz_orth");
        exit(0); }
    H = 0.;
    total_ex = 0.;
    total_col = 0.;
    total_nuc = 0.;
    total_enuc = 0.;
    total_kinet = 0.;
    total_norm = 0.;
    /* force the expansion to be current and normal */
    o_update_normals();
    renormalize();
    ritzn = 0.;
    ritzd = 1.;
    total_norm = ritzd;
    total_nuc = nuclear();
    H = total_nuc;

    /* orthogonalize the bases */
    orth_invert( gscoef,overlaps,o_cache,no);
    /* calculate the trace */
    o4 = o_next(-1);
    for( l=0; l<no; l++)
    {
        o1 = o_next(-1);
        for( i=0; i<= l; i++)
        {
            x = 2.;
            if( i==l) x = 1.;

            cv = x*orth_nuclear(o1,o4);
            ev = x*orth_kinet(o1,o4);
            total_enuc += cv*(*gscoef)[l*no+i];
            total_kinet += ev*(*gscoef)[l*no+i];

            o2 = o_next(-1);
            for( j=0; j< no; j++)
            {
                o3 = o_next(-1);
                for( k=0; k <no  ; k++)
                {
                    cv = orth_coloumb(o1,o2,o3,o4);
                    total_col += x*cv*2.*(*gscoef)[l*no+i]*(*gscoef)[k*no+j];
                    y = (*gscoef)[l*no+k]*(*gscoef)[i*no+j];
                    if( i != l)
                        y += (*gscoef)[i*no+k]*(*gscoef)[l*no+j];
                    /*
                    	y = (*gscoef)[l*no+j]*(*gscoef)[k*no+i];
                    */
                    total_ex  -= cv*y;
                    o3 = o3->next;
                }
                o2 = o2->next;
            }/* j */
            o1 = o1->next;
        }/* i */
        o4 = o4->next;
    }/* l */
    H = total_nuc+total_kinet+total_enuc+total_col+total_ex;

    free( o_cache);
    free( overlaps);
    free( gscoef);
    return H;
}/* end of Ritz_orth() */

float Ritz_orth_INDO(cutoff)
float cutoff;
{
    float (*gscoef)[];
    float (*o_cache)[];
    float (*overlaps)[];
    float H,ritzd,ritzn;
    int no,o_number();
    int i,j,k,l,m,n;
    float cv,ev,x;
    float orth_kinet();
    float orth_nuclear();
    float orth_exchange();
    float orth_coloumb();
    float nuclear(),phiphi();
    void o_update_normals(),orth_invert();
    ORBITAL *o1,*o2,*o3,*o4, *o_next();

    no = o_number();
    if(no < 1) return 0.;
    i = no*no;
    gscoef = malloc( i*sizeof(float));
    if( gscoef == NULL ) {
        aaerror("cannot allocate memory in Ritz_orth");
        exit(0); }
    overlaps = malloc( i*sizeof(float));
    if( overlaps == NULL ) {
        aaerror("cannot allocate memory in Ritz_orth");
        exit(0); }
    o_cache = malloc( i*sizeof(float));
    if( o_cache == NULL ) {
        aaerror("cannot allocate memory in Ritz_orth");
        exit(0); }
    H = 0.;
    total_ex = 0.;
    total_col = 0.;
    total_nuc = 0.;
    total_enuc = 0.;
    total_kinet = 0.;
    total_norm = 0.;
    /* force the expansion to be current and normal */
    o_update_normals();
    renormalize();
    ritzn = 0.;
    ritzd = 1.;
    total_norm = ritzd;
    total_nuc = nuclear();
    H = total_nuc;

    /* orthogonalize the bases */
    orth_invert( gscoef,overlaps,o_cache,no);
    /* calculate the trace */
    o4 = o_next(-1);
    for( l=0; l<no; l++)
    {
        o1 = o_next(-1);
        for( i=0; i< no; i++)
        {

            cv = orth_nuclear(o1,o4);
            ev = orth_kinet(o1,o4);
            total_enuc += cv*(*gscoef)[l*no+i];
            total_kinet += ev*(*gscoef)[l*no+i];

            /* this makes it INDO */
            if( (*o_cache)[l*no+i] > cutoff){
                o2 = o_next(-1);
                for( j=0; j< no; j++)
                {
                    o3 = o_next(-1);
                    for( k=0; k <no  ; k++)
                    {
                        /* this makes it INDO */
                        if( (*o_cache)[k*no+j] > cutoff){
                            cv = orth_coloumb(o1,o2,o3,o4);
                            total_col += cv*2.*(*gscoef)[l*no+i]*(*gscoef)[k*no+j];
                            total_ex  -= cv*(*gscoef)[l*no+k]*(*gscoef)[i*no+j];
                        }
                        o3 = o3->next;
                    }
                    o2 = o2->next;
                }/* j */
            } /* outer INDO */
            o1 = o1->next;
        }/* i */
        o4 = o4->next;
    }/* l */
    H = total_nuc+total_kinet+total_enuc+total_col+total_ex;

    free( o_cache);
    free( overlaps);
    free( gscoef);
    return H;
}/* end of Ritz_orth() */


void orth_invert( gs,ol,oc,no)
float(*gs)[];
float(*ol)[];
float(*oc)[];
int no;
{
    int i,j,k;
    float x,y;
    float phiphi();
    ORBITAL *o1,*o2,*o_next();

    for( i=0; i< no*no; i ++)
    { (*gs)[i] = 0.; (*ol)[i] = 0; (*oc)[i] = 0.;}

    for( i=0; i< no; i++)
    {
        (*gs)[i*no+i] = 1.;
        (*ol)[i*no+i] = 1.;
        (*oc)[i*no+i] = 1.;
    }

    o1 = o_next(-1);
    for( i=0;i< no;i++)
    {
        if( o1->ncouple > 0){
            o2 = o_next(-1);
            for( j=0; j< i+1; j++)
            {
                if( o2->ncouple > 0 ){
                    (*ol)[i*no+j] = phiphi(o1,o2);
                    (*oc)[i*no+j] = (*ol)[i*no+j];
                    (*ol)[j*no+i] = (*ol)[i*no+j];
                    (*oc)[j*no+i] = (*ol)[i*no+j];
                }
                o2 = o2->next;
            }
        }
        o1 = o1->next;
    }

    /* gauss jacobi inverse
    * rather ugly but
    * a) we're diagonally dominant
    * b) it works
    * c) matrix inverse is unique
    */

    /* to upper triagonal */
    for( i=0; i< no-1; i++)
    {
        x = 1./(*ol)[i*no+i];
        for( j=i+1; j<no; j++)
        {
            y = (*ol)[i*no+j]*x;
            for( k=0; k<no; k++)
            {
                (*ol)[k*no+j] -= y*(*ol)[k*no+i];
                (*gs)[k*no+j] -= y*(*gs)[k*no+i];
            }
        }
    }
    /* backsubs */

    for( i=no-1; i>=0; i--)
    {
        x = 1./(*ol)[i*no+i];
        for( k=0; k<no; k++)
        {
            (*gs)[k*no+i] *= x;
            (*ol)[k*no+i] *= x;
        }

        for( j=i-1; j>=0; j--)
        {
            y = (*ol)[i*no+j];
            for( k=0; k<no; k++)
            {
                (*ol)[k*no+j] -= y*(*ol)[k*no+i];
                (*gs)[k*no+j] -= y*(*gs)[k*no+i];
            }
        }
    }

    /* debugging
    	for( i=0; i< no; i++)
    	{
    	for( j=0; j< no; j++)
    	printf("%f ",(*ol)[i*no+j]);
    	printf("\n");
    	}
    	for( i=0; i< no; i++)
    	{
    	for( j=0; j< no; j++)
    	printf("%f ",(*gs)[i*no+j]);
    	printf("\n");
    	}
    */

}/* end of orth_invert() */

void orth_orb( gs,ol,ok,no)
float(*gs)[];
float(*ol)[];
float(*ok)[];
int no;
{
    float phiphi();
    ORBITAL *o1,*o2,*o_next();
    int i,j,k,l;

    for( i=0; i< no*no; i++)
    { (*ol)[i] = 0.;(*ok)[i] = 0.; (*gs)[i] = 0.;}
    for( i=0; i< no; i++)
        (*ol)[i*no+i] = 1.;

    o1 = o_next(-1);
    for( i=0;i< no;i++)
    {
        if( o1->ncouple > 0){
            o2 = o_next(-1);
            for( j=0; j< i+1; j++)
            {
                if( o2->ncouple > 0 ){
                    (*ol)[i*no+j] = phiphi(o1,o2);
                    (*ol)[j*no+i] = (*ol)[i*no+j];
                    (*ok)[i*no+j] = (*ol)[i*no+j];
                    (*ok)[j*no+i] = (*ol)[i*no+j];
                }
                o2 = o2->next;
            }
        }
        o1 = o1->next;
    }
    /* form the orthogonal basis using similarity
    * transformations  (Jacobi's method)
    * we need this basis for the overlap and coloumb integrals
    * but not for the others!
    */
    /* jacobi is in normal.c */

    jacobi( ol,gs,no,no*no*100,1.e-7);
    for(i=0; i< no; i++)
        printf("%f ",(*ol)[i*no+i]);
    printf("\n");
}




float orth_kinet(oo1,oo2,prenorm)
ORBITAL *oo1,*oo2;
float prenorm;
{

    ORBITAL *o1,*o2;
    int i,j,kk,kkk,l,m;
    int ll;
    float x,y,z,r;
    float x1,y1,z1,r1;
    float x2,y2,z2,r2;
    float ke,nuc;
    float rab,rcp,rpq,re,rs;
    float s4,s12,s34,ntemp;
    float Fzero();
    float phiphi();
    ATOM *ac,*a_next();
    float s1,s2,sign;
    int k,numatm,a_number();

    ke = 0.;
    nuc = 0.;
    if( oo1->ncouple < 1) return 0.;
    if( oo2->ncouple < 1) return 0.;

    for( ll=0;ll< oo2->ncouple; ll++)
    { o2 = oo2->couple[ll];
        s2 = 1;
        for( l = 0; l< o2->ncenter; l++)
        {
            if( o2->type == Or4p) s2 = -s2;
            for( kkk=0; kkk<oo1->ncouple; kkk++)
            {
                o1 = oo1->couple[kkk];
                s1 = 1;
                for( k = 0; k< o1->ncenter; k++)
                {
                    if( o1->type == Or4p) s1 = -s1;
                    sign = s1*s2;
                    x = o2->rx[l] - o1->rx[k];
                    y = o2->ry[l] - o1->ry[k];
                    z = o2->rz[l] - o1->rz[k];
                    rab = (x*x + y*y + z*z)*INVBOHR*INVBOHR;
                    numatm = a_number();
                    for( i=0; i< o1->n; i++)
                    {
                        for( j=0; j< o2->n; j++)
                        {

                            r = o1->r[i]*o2->r[j]/( o1->r[i] + o2->r[j]);
                            re = exp(-r*rab);
                            rs = (o1->r[i] + o2->r[j]);
                            /* x1,y1,z1 are in ANGSTROMS */
                            ke = ke +  sign*o1->a[i]*o2->a[j]*
                                 r*(three-two*r*rab)*pow( PI/rs,1.5)*re;
                        }
                    }
                }} /* k,kkk loops */
        }} /* ll,l */
    if( oo1->ipair == 2 && oo2->ipair == 2) {
        ke = 2*ke;
    } else {
        ke = ke;
    }

    return (ke);
}

float orth_nuclear(oo1,oo2,prenorm)
ORBITAL *oo1,*oo2;
float prenorm;
{

    ORBITAL *o1,*o2;
    int i,j,kk,kkk,l,m;
    int ll;
    float x,y,z,r;
    float x1,y1,z1,r1;
    float x2,y2,z2,r2;
    float ke,nuc;
    float rab,rcp,rpq,re,rs;
    float s4,s12,s34,ntemp;
    float Fzero();
    float phiphi();
    ATOM *ac,*a_next();
    float s1,s2,sign;
    int k,numatm,a_number();

    ke = 0.;
    nuc = 0.;
    if( oo1->ncouple < 1) return 0.;
    if( oo2->ncouple < 1) return 0.;

    for( ll=0;ll< oo2->ncouple; ll++)
    { o2 = oo2->couple[ll];
        s2 = 1;
        for( l = 0; l< o2->ncenter; l++)
        {
            if( o2->type == Or4p) s2 = -s2;
            for( kkk=0; kkk<oo1->ncouple; kkk++)
            {
                o1 = oo1->couple[kkk];
                s1 = 1;
                for( k = 0; k< o1->ncenter; k++)
                {
                    if( o1->type == Or4p) s1 = -s1;
                    sign = s1*s2;
                    x = o2->rx[l] - o1->rx[k];
                    y = o2->ry[l] - o1->ry[k];
                    z = o2->rz[l] - o1->rz[k];
                    rab = (x*x + y*y + z*z)*INVBOHR*INVBOHR;
                    numatm = a_number();
                    for( i=0; i< o1->n; i++)
                    {
                        for( j=0; j< o2->n; j++)
                        {

                            r = o1->r[i]*o2->r[j]/( o1->r[i] + o2->r[j]);
                            re = exp(-r*rab);
                            rs = (o1->r[i] + o2->r[j]);
                            /* nuclear integrals including partial charges for at a distance */
                            x1 = (o2->rx[l]*o2->r[j]+o1->rx[k]*o1->r[i])/rs;
                            y1 = (o2->ry[l]*o2->r[j]+o1->ry[k]*o1->r[i])/rs;
                            z1 = (o2->rz[l]*o2->r[j]+o1->rz[k]*o1->r[i])/rs;
                            ntemp =  sign*o1->a[i]*o2->a[j]*TWOPI/rs*re;
                            for(kk=0; kk< numatm; kk++)
                            { ac = a_next(kk);
                                /* ANGSTROM difference */
                                x2 = x1 - ac->x;
                                y2 = y1 - ac->y;
                                z2 = z1 - ac->z;
                                /* converted to BOHR here */
                                rcp = (x2*x2 +y2*y2 + z2*z2)*INVBOHR*INVBOHR;
                                if( ac->na > 0 )
                                    nuc = nuc + ntemp*ac->na *Fzero(rs*rcp);
                                else
                                {
                                    for( m=0; m< o1->myatom->dontuse; m++)
                                        if( ac == o1->myatom->excluded[m] )  goto SKIP;
                                    for( m=0; m< o2->myatom->dontuse; m++)
                                        if( ac == o2->myatom->excluded[m] )  goto SKIP;
                                    nuc = nuc - ntemp*ac->q *Fzero(rs*rcp);
SKIP:		;
                                }/*else */
                            } /* kk */
                        }
                    }
                }} /* k,kkk loops */
        }} /* ll,l */
    if( oo1->ipair == 2 && oo2->ipair == 2) {
        nuc = 2*nuc;
    } else {
        nuc = nuc;
    }

    return (-nuc );
}
float orth_justone_nuclear(ac,oo1,oo2,prenorm)
ATOM *ac;
ORBITAL *oo1,*oo2;
float prenorm;
{

    ORBITAL *o1,*o2;
    int i,j,kk,kkk,l,m;
    int ll;
    float x,y,z,r;
    float x1,y1,z1,r1;
    float x2,y2,z2,r2;
    float ke,nuc;
    float rab,rcp,rpq,re,rs;
    float s4,s12,s34,ntemp;
    float Fzero();
    float phiphi();
    float s1,s2,sign;
    int k,numatm,a_number();

    ke = 0.;
    nuc = 0.;
    if( oo1->ncouple < 1) return 0.;
    if( oo2->ncouple < 1) return 0.;

    for( ll=0;ll< oo2->ncouple; ll++)
    { o2 = oo2->couple[ll];
        s2 = 1;
        for( l = 0; l< o2->ncenter; l++)
        {
            if( o2->type == Or4p) s2 = -s2;
            for( kkk=0; kkk<oo1->ncouple; kkk++)
            {
                o1 = oo1->couple[kkk];
                s1 = 1;
                for( k = 0; k< o1->ncenter; k++)
                {
                    if( o1->type == Or4p) s1 = -s1;
                    sign = s1*s2;
                    x = o2->rx[l] - o1->rx[k];
                    y = o2->ry[l] - o1->ry[k];
                    z = o2->rz[l] - o1->rz[k];
                    rab = (x*x + y*y + z*z)*INVBOHR*INVBOHR;
                    numatm = a_number();
                    for( i=0; i< o1->n; i++)
                    {
                        for( j=0; j< o2->n; j++)
                        {

                            r = o1->r[i]*o2->r[j]/( o1->r[i] + o2->r[j]);
                            re = exp(-r*rab);
                            rs = (o1->r[i] + o2->r[j]);
                            /* nuclear integrals including partial charges for at a distance */
                            x1 = (o2->rx[l]*o2->r[j]+o1->rx[k]*o1->r[i])/rs;
                            y1 = (o2->ry[l]*o2->r[j]+o1->ry[k]*o1->r[i])/rs;
                            z1 = (o2->rz[l]*o2->r[j]+o1->rz[k]*o1->r[i])/rs;
                            ntemp =  sign*o1->a[i]*o2->a[j]*TWOPI/rs*re;

                            /* ANGSTROM difference */
                            x2 = x1 - ac->x;
                            y2 = y1 - ac->y;
                            z2 = z1 - ac->z;
                            /* converted to BOHR here */
                            rcp = (x2*x2 +y2*y2 + z2*z2)*INVBOHR*INVBOHR;
                            if( ac->na > 0 )
                                nuc = nuc + ntemp*ac->na *Fzero(rs*rcp);
                            else
                            {
                                for( m=0; m< o1->myatom->dontuse; m++)
                                    if( ac == o1->myatom->excluded[m] )  goto SKIP;
                                for( m=0; m< o2->myatom->dontuse; m++)
                                    if( ac == o2->myatom->excluded[m] )  goto SKIP;
                                nuc = nuc - ntemp*ac->q *Fzero(rs*rcp);
SKIP:		;
                            }/*else */
                        }
                    }
                }} /* k,kkk loops */
        }} /* ll,l */
    if( oo1->ipair == 2 && oo2->ipair == 2) {
        nuc = 2*nuc;
    } else {
        nuc = nuc;
    }

    return (-nuc);
}

float orth_coloumb(oo1,oo2,oo3,oo4,prenorm)
ORBITAL *oo1,*oo4;
ORBITAL *oo2,*oo3;
float prenorm;
{

    /* coloumb and exchange integrals are 4-center uugh!! */
    /* do all the integrals
    * o1 is assumed to be o4 
    o1   o2 o3 o4(4)
    and
    o1  o2 o3(4) o4(3)
    */
    /* just Hartree version !!! */
    int i,j,m,n,o,p;
    float x,y,z,r;
    float x1,y1,z1,r1;
    float x2,y2,z2,r2;
    float coloumb,exchange;
    float rab,rbc,rac,rcd,rpq,r3,r4,rs;
    float sa4,s12,s34;
    float Fzero();
    int k,l,no,o_number();
    int io3,io2;
    int ioo1,ioo4;
    int ioo3,ioo2;
    ORBITAL *o1,*o4,*o3,*o2,*o_next();
    float  phiphi();
    float s1,s2,s3,s4,sign;

    no = o_number();
    if( no == 1  && oo1->ipair == 1 ) return 0.;
    coloumb = 0.;
    exchange = 0.;

    if( oo1->ncouple < 1) return 0;
    if( oo4->ncouple < 1) return 0;
    if( oo2->ncouple < 1) return 0;
    if( oo3->ncouple < 1) return 0;
    for( ioo1=0; ioo1< oo1->ncouple; ioo1++)
    {
        o1 = oo1->couple[ioo1];
        for( ioo4=0; ioo4< oo4->ncouple; ioo4++)
        {
            o4 = oo4->couple[ioo4];
            for( io3=0; io3< oo3->ncouple; io3++)
            {
                o3 = oo3->couple[io3];
                for( io2=0; io2< oo2->ncouple; io2++)
                {
                    o2 = oo2->couple[io2];
                    /* o1 == o4 so only one overlap */


                    /*
                    	if there is one function then the coloumb and exchange
                    	integrals are the same and the 2*c -ex = 1*c
                    */
                    s1 = 1;
                    for( m = 0 ; m < o1->ncenter; m++)
                    {
                        if( o1->type == Or4p) s1 = -s1;
                        s2 = 1;
                        for( n = 0 ; n < o2->ncenter; n++)
                        {
                            s3 = 1;
                            if( o2->type == Or4p) s2 = -s2;
                            for( o = 0 ; o < o3->ncenter; o++)
                            {
                                if( o3->type == Or4p) s3 = -s3;
                                s4 = 1;
                                for( p = 0 ; p < o4->ncenter; p++)
                                {
                                    if( o4->type == Or4p) s4 = -s4;

                                    sign = s1*s2*s3*s4;
                                    x1 = o1->rx[m] -o4->rx[p];
                                    y1 = o1->ry[m] -o4->ry[p];
                                    z1 = o1->rz[m] -o4->rz[p];
                                    rab = (x1*x1 + y1*y1 + z1*z1)*INVBOHR*INVBOHR;
                                    x1 = o2->rx[n] -o3->rx[o];
                                    y1 = o2->ry[n] -o3->ry[o];
                                    z1 = o2->rz[n] -o3->rz[o];
                                    rac = (x1*x1 + y1*y1 + z1*z1)*INVBOHR*INVBOHR;
                                    x1 = o1->rx[m] -o3->rx[o];
                                    y1 = o1->ry[m] -o3->ry[o];
                                    z1 = o1->rz[m] -o3->rz[o];
                                    rbc = (x1*x1 + y1*y1 + z1*z1)*INVBOHR*INVBOHR;
                                    x1 = o4->rx[p] -o2->rx[n];
                                    y1 = o4->ry[p] -o2->ry[n];
                                    z1 = o4->rz[p] -o2->rz[n];
                                    rcd = (x1*x1 + y1*y1 + z1*z1)*INVBOHR*INVBOHR;
                                    for( i=0; i< o1->n; i++)
                                    {
                                        for( j=0; j< o2->n; j++)
                                        {
                                            for( k=0; k< o3->n; k++)
                                            {
                                                for( l=0; l< o4->n; l++)
                                                {
                                                    /*  1,2  3,4 coloumb */
                                                    /* this is the correct form at least for 2 e- */

                                                    s12 = o1->r[i] + o4->r[l];
                                                    s34 = o3->r[k] + o2->r[j];
                                                    sa4  = s12 + s34;
                                                    r1 = o1->r[i]*o4->r[l]/( o1->r[i]+o4->r[l]);
                                                    r2 = o3->r[k]*o2->r[j]/( o3->r[k]+o2->r[j]);
                                                    x1 = (o1->r[i]*o1->rx[m] + o4->r[l]*o4->rx[p])/s12;
                                                    y1 = (o1->r[i]*o1->ry[m] + o4->r[l]*o4->ry[p])/s12;
                                                    z1 = (o1->r[i]*o1->rz[m] + o4->r[l]*o4->rz[p])/s12;
                                                    x2 = (o3->r[k]*o3->rx[o] + o2->r[j]*o2->rx[n])/s34;
                                                    y2 = (o3->r[k]*o3->ry[o] + o2->r[j]*o2->ry[n])/s34;
                                                    z2 = (o3->r[k]*o3->rz[o] + o2->r[j]*o2->rz[n])/s34;
                                                    x2 -= x1;
                                                    y2 -= y1;
                                                    z2 -= z1;
                                                    rpq = (x2*x2 + y2*y2 + z2*z2)*INVBOHR*INVBOHR;
                                                    /* rab == 0 , rcd == 0 so just need rac or rbd or ... */
                                                    x2 = r1*rab+r2*rac;
                                                    if( x2 < 70.)
                                                        coloumb = coloumb + sign*
                                                                  o1->a[i]*o2->a[j]*o3->a[k]*o4->a[l] *
                                                                  34.986837 /(s12*s34*sqrt(sa4))*Fzero( s12*s34/sa4*rpq)
                                                                  *exp( -x2);
                                                } /* l */
                                            }/*k */
                                        } /* j */
                                    } /* i */
                                }}}} /* m,n,o,p loops */
                } /* io2 */
            } /* io3 */
        } /* ioo4 */
    } /* ioo1 */

    return (coloumb  );
}
float orth_exchange(oo1,oo2,oo3,oo4,prenorm)
ORBITAL *oo1,*oo4;
ORBITAL *oo2,*oo3;
float prenorm;
{

    /* coloumb and exchange integrals are 4-center uugh!! */
    /* do all the integrals
    * o1 is assumed to be o4 
    o1   o2 o3 o4(4)
    and
    o1  o2 o3(4) o4(3)
    */
    /* just Hartree version !!! */
    int i,j,m,n,o,p;
    float x,y,z,r;
    float x1,y1,z1,r1;
    float x2,y2,z2,r2;
    float coloumb,exchange;
    float rab,rbc,rac,rcd,rpq,r3,r4,rs;
    float sa4,s12,s34;
    float Fzero();
    int k,l,no,o_number();
    int io3,io2;
    int ioo1,ioo4;
    int ioo3,ioo2;
    ORBITAL *o1,*o4,*o3,*o2,*o_next();
    float phiphi();
    float s1,s2,s3,s4,sign;

    no = o_number();
    if( no == 1  && oo1->ipair == 1 ) return 0.;
    if( oo1 != oo4)
    {
        fprintf(stderr,"invalid call to phiH4phi\n");
        return 0.;
    }
    coloumb = 0.;
    exchange = 0.;

    if( oo1->ncouple < 1) return 0;
    if( oo4->ncouple < 1) return 0;
    if( oo2->ncouple < 1) return 0;
    if( oo3->ncouple < 1) return 0;
    for( ioo1=0; ioo1< oo1->ncouple; ioo1++)
    {
        o1 = oo1->couple[ioo1];
        for( ioo4=0; ioo4< oo4->ncouple; ioo4++)
        {
            o4 = oo4->couple[ioo4];
            for( io3=0; io3< oo3->ncouple; io3++)
            {
                o3 = oo3->couple[io3];
                for( io2=0; io2< oo2->ncouple; io2++)
                {
                    o2 = oo2->couple[io2];
                    /* o1 == o4 so only one overlap */


                    /*
                    	if there is one function then the coloumb and exchange
                    	integrals are the same and the 2*c -ex = 1*c
                    */
                    s1 = 1;
                    for( m = 0 ; m < o1->ncenter; m++)
                    {
                        if( o1->type == Or4p) s1 = -s1;
                        s2 = 1;
                        for( n = 0 ; n < o2->ncenter; n++)
                        {
                            s3 = 1;
                            if( o2->type == Or4p) s2 = -s2;
                            for( o = 0 ; o < o3->ncenter; o++)
                            {
                                if( o3->type == Or4p) s3 = -s3;
                                s4 = 1;
                                for( p = 0 ; p < o4->ncenter; p++)
                                {
                                    if( o4->type == Or4p) s4 = -s4;

                                    sign = s1*s2*s3*s4;
                                    x1 = o1->rx[m] -o4->rx[p];
                                    y1 = o1->ry[m] -o4->ry[p];
                                    z1 = o1->rz[m] -o4->rz[p];
                                    rab = (x1*x1 + y1*y1 + z1*z1)*INVBOHR*INVBOHR;
                                    x1 = o2->rx[n] -o3->rx[o];
                                    y1 = o2->ry[n] -o3->ry[o];
                                    z1 = o2->rz[n] -o3->rz[o];
                                    rac = (x1*x1 + y1*y1 + z1*z1)*INVBOHR*INVBOHR;
                                    x1 = o1->rx[m] -o3->rx[o];
                                    y1 = o1->ry[m] -o3->ry[o];
                                    z1 = o1->rz[m] -o3->rz[o];
                                    rbc = (x1*x1 + y1*y1 + z1*z1)*INVBOHR*INVBOHR;
                                    x1 = o4->rx[p] -o2->rx[n];
                                    y1 = o4->ry[p] -o2->ry[n];
                                    z1 = o4->rz[p] -o2->rz[n];
                                    rcd = (x1*x1 + y1*y1 + z1*z1)*INVBOHR*INVBOHR;
                                    for( i=0; i< o1->n; i++)
                                    {
                                        for( j=0; j< o2->n; j++)
                                        {
                                            for( k=0; k< o3->n; k++)
                                            {
                                                for( l=0; l< o4->n; l++)
                                                {
                                                    /*  1,2  3,4 coloumb */
                                                    /* this is the correct form at least for 2 e- */

                                                    s12 = o1->r[i] + o4->r[l];
                                                    s34 = o3->r[k] + o2->r[j];
                                                    sa4  = s12 + s34;
                                                    r1 = o1->r[i]*o4->r[l]/( o1->r[i]+o4->r[l]);
                                                    r2 = o3->r[k]*o2->r[j]/( o3->r[k]+o2->r[j]);
                                                    s12 = o1->r[i] + o3->r[k];
                                                    s34 = o4->r[l] + o2->r[j];
                                                    sa4  = s12 + s34;
                                                    x1 = (o1->r[i]*o1->rx[m] + o3->r[k]*o3->rx[o])/s12;
                                                    y1 = (o1->r[i]*o1->ry[m] + o3->r[k]*o3->ry[o])/s12;
                                                    z1 = (o1->r[i]*o1->rz[m] + o3->r[k]*o3->rz[o])/s12;
                                                    x2 = (o4->r[l]*o4->rx[p] + o2->r[j]*o2->rx[n])/s34;
                                                    y2 = (o4->r[l]*o4->ry[p] + o2->r[j]*o2->ry[n])/s34;
                                                    z2 = (o4->r[l]*o4->rz[p] + o2->r[j]*o2->rz[n])/s34;
                                                    x2 -= x1;
                                                    y2 -= y1;
                                                    z2 -= z1;
                                                    rpq = (x2*x2 + y2*y2 + z2*z2)*INVBOHR*INVBOHR;
                                                    r1 = o1->r[i]*o3->r[k]/( o1->r[i]+o3->r[k]);
                                                    r2 = o4->r[l]*o2->r[j]/( o4->r[l]+o2->r[j]);
                                                    x2 = r1*rbc + r2*rcd;
                                                    if( x2 < 70.)
                                                        coloumb = coloumb - sign*
                                                                  o1->a[i]*o2->a[j]*o3->a[k]*o4->a[l] *
                                                                  34.986837 /(s12*s34*sqrt(sa4))*Fzero( s12*s34/sa4*rpq)
                                                                  *exp( -x2);


                                                } /* l */
                                            }/*k */
                                        } /* j */
                                    } /* i */
                                }}}} /* m,n,o,p loops */
                } /* io2 */
            } /* io3 */
        } /* ioo4 */
    } /* ioo1 */

    /*
    float total_ex,total_col,total_nuc,total_enuc,total_kinet,total_norm;
    */
    return (coloumb  );
}

#define fast
void Ritz_orth_grad(grad,where,params,whatO,ngrad)
float grad[];
float *where[];
float params[];
ORBITAL *whatO[];
int ngrad;
{
    float (*gscoef)[];
    float (*o_cache)[];
    float (*overlaps)[];
    float (*twocent)[];
    float (*fourcent)[];
    float H,ritzd,ritzn;
    int no,o_number();
    int i,j,k,l,m,n;
    int igrad;
    float cv,ev,x,y;
    float orth_kinet();
    float orth_nuclear();
    float orth_exchange();
    float orth_coloumb();
    float nuclear(),phiphi();
    void o_update_normals(),orth_invert();
    ORBITAL *o1,*o2,*o3,*o4, *o_next();

    if( ngrad < 1) return ;

    no = o_number();
    if(no < 1) return ;
    i = no*no;
    gscoef = malloc( i*sizeof(float));
    if( gscoef == NULL ) {
        aaerror("cannot allocate memory in Ritz_orth");
        exit(0); }
    overlaps = malloc( i*sizeof(float));
    if( overlaps == NULL ) {
        aaerror("cannot allocate memory in Ritz_orth");
        exit(0); }
    o_cache = malloc( i*sizeof(float));
    if( o_cache == NULL ) {
        aaerror("cannot allocate memory in Ritz_orth");
        exit(0); }
    twocent = malloc( i*sizeof(float));
    if( twocent == NULL ) {
        aaerror("cannot allocate memory in Ritz_orth");
        exit(0); }
    fourcent = malloc( i*i*sizeof(float));
    if( fourcent == NULL ) {
        aaerror("cannot allocate memory in Ritz_orth");
        exit(0); }
    total_ex = 0.;
    total_col = 0.;
    total_nuc = 0.;
    total_enuc = 0.;
    total_kinet = 0.;
    total_norm = 0.;
    /* force the expansion to be current and normal */
    o_update_normals();
    renormalize();
    ritzn = 0.;
    ritzd = 1.;
    /* calculate the four center integrals */
    o4 = o_next(-1);
    for( l=0; l< no; l++)
    {
        o1 = o_next(-1);
        for( i=0; i<= l; i++)
        {
            cv = orth_nuclear(o1,o4)+ orth_kinet(o1,o4);
            (*twocent)[l*no+i] = cv;
            (*twocent)[i*no+l] = cv;
            o2 = o_next(-1);
            for( j=0; j< no; j++)
            {
                o3 = o_next(-1);
                for( k=0; k< no; k++)
                {
                    cv = orth_coloumb(o1,o2,o3,o4);
                    (*fourcent)[ ((l*no+i)*no+j)*no+k] = cv;
                    (*fourcent)[ ((i*no+l)*no+j)*no+k] = cv;
                    o3 = o3->next;
                } /*k */
                o2 = o2->next;
            } /*j */
            o1 = o1->next;
        } /*i */
        o4 = o4->next;
    } /*l */

    /* orthogonalize the bases */
    for( igrad = 0; igrad < ngrad; igrad++)
    {
        grad[igrad] = 0.;
        total_kinet = 0.;
        total_col = 0.;
        total_ex = 0.;

        for( i=0; i< ngrad; i++)
            *where[i] = params[i];

        *where[igrad] -= orth_orb_dstep;
        o_update_normals();
        renormalize();
        grad[igrad] -= nuclear();
        orth_invert( gscoef,overlaps,o_cache,no);
        /* calculate the trace */
        o4 = o_next(-1);
        for( l=0; l<no; l++)
        {
            o1 = o_next(-1);
            for( i=0; i<= l; i++)
            {
                x = 2.;
                if( i == l ) x = 1.;
#ifdef fast
                if( whatO[igrad] == o1 || whatO[igrad] == o4 ||
                        whatO[igrad] == o1->gang || whatO[igrad] == o4->gang ){
#endif
                    cv = orth_nuclear(o1,o4)+ orth_kinet(o1,o4);
#ifdef fast
                } else {
                    cv = (*twocent)[l*no+i];
                }
#endif
                total_kinet -= x*cv*(*gscoef)[l*no+i];

                o2 = o_next(-1);
                for( j=0; j< no; j++)
                {

                    o3 = o_next(-1);
                    for( k=0; k <no  ; k++)
                    {
#ifdef fast
                        if( whatO[igrad] == o1 || whatO[igrad] == o2 ||
                                whatO[igrad] == o3 || whatO[igrad] == o4
                                || whatO[igrad] == o1->gang || whatO[igrad] == o2->gang ||
                                whatO[igrad] == o3->gang || whatO[igrad] == o4->gang ){
#endif
                            cv = orth_coloumb(o1,o2,o3,o4);
#ifdef fast
                        }else{
                            cv = (*fourcent)[ ((l*no+i)*no+j)*no+k] ;
                        }
#endif
                        total_col -= x*cv*2.*(*gscoef)[l*no+i]*(*gscoef)[k*no+j];
                        total_ex += cv*(*gscoef)[l*no+k]*(*gscoef)[i*no+j];
                        if( i != l )
                            total_ex += cv*(*gscoef)[i*no+k]*(*gscoef)[l*no+j];
                        o3 = o3->next;
                    }
                    o2 = o2->next;
                }/* j */
                o1 = o1->next;
            }/* i */
            o4 = o4->next;
        }/* l */

        for( i=0; i< ngrad; i++)
            *where[i] = params[i];
        *where[igrad] += orth_orb_dstep;

        o_update_normals();
        renormalize();
        grad[igrad] += nuclear();
        orth_invert( gscoef,overlaps,o_cache,no);
        /* calculate the trace */
        o4 = o_next(-1);
        for( l=0; l<no; l++)
        {
            o1 = o_next(-1);
            for( i=0; i<= l; i++)
            {
                x = 2.;
                if( i == l ) x = 1.;

#ifdef fast
                if( o1 == whatO[igrad] || o4 == whatO[igrad]
                        || o1->gang == whatO[igrad] || o4->gang == whatO[igrad]){
#endif
                    cv = orth_nuclear(o1,o4)+ orth_kinet(o1,o4);
#ifdef fast
                } else {
                    cv = (*twocent)[l*no+i];
                }
#endif
                total_kinet += x*cv*(*gscoef)[l*no+i];

                o2 = o_next(-1);
                for( j=0; j< no; j++)
                {
                    o3 = o_next(-1);
                    for( k=0; k < no ; k++)
                    {
#ifdef fast
                        if( whatO[igrad] == o1 || whatO[igrad] == o2 ||
                                whatO[igrad] == o3 || whatO[igrad] == o4
                                || whatO[igrad] == o1->gang || whatO[igrad] == o2->gang ||
                                whatO[igrad] == o3->gang || whatO[igrad] == o4->gang ){
#endif
                            cv = orth_coloumb(o1,o2,o3,o4);
#ifdef fast
                        }else{
                            cv = (*fourcent)[ ((l*no+i)*no+j)*no+k] ;
                        }
#endif
                        total_col += x*cv*2.*(*gscoef)[l*no+i]*(*gscoef)[k*no+j];
                        total_ex -= cv*(*gscoef)[l*no+k]*(*gscoef)[i*no+j];
                        if( i != l )
                            total_ex -= cv*(*gscoef)[i*no+k]*(*gscoef)[l*no+j];
                        o3 = o3->next;
                    }
                    o2 = o2->next;
                }/* j */
                o1 = o1->next;
            }/* i */
            o4 = o4->next;
        }/* l */
        grad[igrad] = total_kinet + total_col + total_ex;
        grad[igrad] *= 	-0.5e2;
    }/* igrad */

    /*
    * also done in dscf.c
    	for( i=0; i< ngrad; i++)
    	*where[i] = params[i];
    */
    free( fourcent);
    free( twocent);
    free( o_cache);
    free( overlaps);
    free( gscoef);
    return;
}/* end of Ritz_orth_grad() */

void Ritz_orth_grad_virial(grad,where,params,whatO,ngrad)
float grad[];
float *where[];
float params[];
ORBITAL *whatO[];
int ngrad;
{
    float (*gscoef)[];
    float (*o_cache)[];
    float (*overlaps)[];
    float (*twocent)[];
    float (*fourcent)[];
    float H,ritzd,ritzn;
    int no,o_number();
    int i,j,k,l,m,n;
    int igrad;
    float cv,ev,x;
    float orth_kinet();
    float orth_nuclear();
    float orth_exchange();
    float orth_coloumb();
    float nuclear(),phiphi();
    void o_update_normals(),orth_invert();
    ORBITAL *o1,*o2,*o3,*o4, *o_next();
    float dV,kin;

    if( ngrad < 1) return ;

    no = o_number();
    if(no < 1) return ;
    i = no*no;
    gscoef = malloc( i*sizeof(float));
    if( gscoef == NULL ) {
        aaerror("cannot allocate memory in Ritz_orth");
        exit(0); }
    overlaps = malloc( i*sizeof(float));
    if( overlaps == NULL ) {
        aaerror("cannot allocate memory in Ritz_orth");
        exit(0); }
    o_cache = malloc( i*sizeof(float));
    if( o_cache == NULL ) {
        aaerror("cannot allocate memory in Ritz_orth");
        exit(0); }
    twocent = malloc( i*sizeof(float));
    if( twocent == NULL ) {
        aaerror("cannot allocate memory in Ritz_orth");
        exit(0); }
    fourcent = malloc( i*i*sizeof(float));
    if( fourcent == NULL ) {
        aaerror("cannot allocate memory in Ritz_orth");
        exit(0); }
    total_ex = 0.;
    total_col = 0.;
    total_nuc = 0.;
    total_enuc = 0.;
    total_kinet = 0.;
    total_norm = 0.;
    /* force the expansion to be current and normal */
    o_update_normals();
    renormalize();
    ritzn = 0.;
    ritzd = 1.;
    /* calculate the four center integrals */
    o4 = o_next(-1);
    for( l=0; l< no; l++)
    {
        o1 = o_next(-1);
        for( i=0; i< no; i++)
        {
            cv = orth_nuclear(o1,o4)+ orth_kinet(o1,o4);
            (*twocent)[l*no+i] = cv;
            o2 = o_next(-1);
            for( j=0; j< no; j++)
            {
                o3 = o_next(-1);
                for( k=0; k< no; k++)
                {
                    cv = orth_coloumb(o1,o2,o3,o4);
                    (*fourcent)[ ((l*no+i)*no+j)*no+k] = cv;
                    o3 = o3->next;
                } /*k */
                o2 = o2->next;
            } /*j */
            o1 = o1->next;
        } /*i */
        o4 = o4->next;
    } /*l */

    /* orthogonalize the bases */
    for( igrad = 0; igrad < ngrad; igrad++)
    {
        grad[igrad] = 0.;
        total_kinet = 0.;
        total_col = 0.;
        total_ex = 0.;

        for( i=0; i< ngrad; i++)
            *where[i] = params[i];
        H = 0.;
        dV = 0.;
        kin = 0.;
        printf("%d\n",(whatO[igrad])->osn);

        *where[igrad] -= orth_orb_dstep;
        o_update_normals();
        renormalize();
        grad[igrad] -= nuclear();
        H       = grad[igrad];
        orth_invert( gscoef,overlaps,o_cache,no);
        /* calculate the trace */
        o4 = o_next(-1);
        for( l=0; l<no; l++)
        {
            o1 = o_next(-1);
            for( i=0; i< no; i++)
            {
#ifdef fast
                if( o1 == whatO[igrad] || o4 == whatO[igrad]){
#endif
                    cv = orth_nuclear(o1,o4)+ orth_kinet(o1,o4);
#ifdef fast
                } else {
                    cv = (*twocent)[l*no+i];
                }
#endif
                total_kinet -= cv*(*gscoef)[l*no+i];
                kin += orth_kinet(o1,o4)*(*gscoef)[l*no+i];

                o2 = o_next(-1);
                for( j=0; j< no; j++)
                {

                    o3 = o_next(-1);
                    for( k=0; k <no  ; k++)
                    {
#ifdef fast
                        if( whatO[igrad] == o1 || whatO[igrad] == o2 ||
                                whatO[igrad] == o3 || whatO[igrad] == o4){
#endif
                            cv = orth_coloumb(o1,o2,o3,o4);
#ifdef fast
                        }else{
                            cv = (*fourcent)[ ((l*no+i)*no+j)*no+k] ;
                        }
#endif
                        total_col -= cv*2.*(*gscoef)[l*no+i]*(*gscoef)[k*no+j];
                        total_ex += cv*(*gscoef)[l*no+k]*(*gscoef)[i*no+j];
                        o3 = o3->next;
                    }
                    o2 = o2->next;
                }/* j */
                o1 = o1->next;
            }/* i */
            o4 = o4->next;
        }/* l */
        H += total_kinet + total_col + total_ex;
        /* H is negative H because we're makeing the gradient */
        dV = ( kin -H)*(kin-H);
        printf("dV %f kin %f H %f \n",dV,kin,H);

        H = 0.;
        kin = 0.;

        for( i=0; i< ngrad; i++)
            *where[i] = params[i];
        *where[igrad] += orth_orb_dstep;

        o_update_normals();
        renormalize();
        ev = nuclear();
        grad[igrad] += ev;
        H += ev ;
        orth_invert( gscoef,overlaps,o_cache,no);
        /* calculate the trace */
        o4 = o_next(-1);
        for( l=0; l<no; l++)
        {
            o1 = o_next(-1);
            for( i=0; i< no; i++)
            {

#ifdef fast
                if( o1 == whatO[igrad] || o4 == whatO[igrad]){
#endif
                    cv = orth_nuclear(o1,o4)+ orth_kinet(o1,o4);
#ifdef fast
                } else {
                    cv = (*twocent)[l*no+i];
                }
#endif
                total_kinet += cv*(*gscoef)[l*no+i];
                kin += orth_kinet(o1,o4)*(*gscoef)[l*no+i];

                o2 = o_next(-1);
                for( j=0; j< no; j++)
                {
                    o3 = o_next(-1);
                    for( k=0; k <no  ; k++)
                    {
#ifdef fast
                        if( whatO[igrad] == o1 || whatO[igrad] == o2 ||
                                whatO[igrad] == o3 || whatO[igrad] == o4){
#endif
                            cv = orth_coloumb(o1,o2,o3,o4);
#ifdef fast
                        }else{
                            cv = (*fourcent)[ ((l*no+i)*no+j)*no+k] ;
                        }
#endif
                        total_col += cv*2.*(*gscoef)[l*no+i]*(*gscoef)[k*no+j];
                        total_ex -= cv*(*gscoef)[l*no+k]*(*gscoef)[i*no+j];
                        o3 = o3->next;
                    }
                    o2 = o2->next;
                }/* j */
                o1 = o1->next;
            }/* i */
            o4 = o4->next;
        }/* l */
        grad[igrad] = total_kinet + total_col + total_ex;
        H += total_kinet + total_col + total_ex;
        dV = ( kin +H)*(kin+H);
        printf("dV %f kin %f H %f \n",dV,kin,H);
        /*
        	grad[igrad] -= dV;
        */
        grad[igrad] *= 	-0.5e2;
    }/* igrad */

    /*
    * also done in dscf.c
    	for( i=0; i< ngrad; i++)
    	*where[i] = params[i];
    */
    free( fourcent);
    free( twocent);
    free( o_cache);
    free( overlaps);
    free( gscoef);
    return;
}/* end of Ritz_orth() */
