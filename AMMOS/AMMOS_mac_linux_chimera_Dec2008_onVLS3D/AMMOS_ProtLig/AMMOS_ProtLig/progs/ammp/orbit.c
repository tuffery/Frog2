/* orbit.c
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

#include <assert.h>

#define ANSI 1
/* misc includes - ANSI and some are just to be safe */
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#ifdef ANSI
#include <stdlib.h>
#endif
#include "ammp.h" /* this includes numeric.h */

#include "orbit.h" 


ORBITAL *firstORBITAL = NULL , *lastORBITAL = NULL;
static int orbitNUMBER = 0;
static int orbitUPDATE = 0;
float total_ex,total_col,total_nuc,total_enuc,total_kinet,total_norm;

/* dummy function call  */
#define o_first() firstORBITAL;

/* dummy calls which call generic orbital routine */
/* this should be done in eval */
/*
int  orb1();
{
}
int  orb1o();
{
}
int  orb2();
{
}
int  orb3();
{
}
int  orb4s();
{
}
int  orb4u();
{
}

*/

int orbital( type, i1,i2,i3,i4,i5,i6, osn, along, x,y,z ,spin,pair)
int type,i1,i2,i3,i4,i5,i6,osn,spin,pair;
float along,x,y,z;
{
    static int highest = -1, lowest = -1;
    ORBITAL *new,*o_m_serial();
    ATOM *ap,*a_m_serial();

    if( pair == 0 ) pair = 2;
    ap = a_m_serial(i1);
    if( ap == NULL ) {
        aaerror("No atom for an orbital\n"); return 0;
    }
    new = NULL;
    if( highest >= osn && lowest <= osn) new = o_m_serial(osn);
    if( highest < osn ) highest = osn;
    if( lowest > osn ) lowest = osn;
    if( new == NULL )
    {
        orbitUPDATE = 1;
        new = malloc( sizeof(ORBITAL));
        if( new == NULL )
        {
            aaerror("cannot allocate memory in orbital\n");
            return 0;
        }
        new->next = NULL;
    }
    if( firstORBITAL == NULL ) {
        firstORBITAL = new;
        highest = osn; lowest = osn; }
    if( lastORBITAL == NULL ) lastORBITAL = new;
    if( new->next == NULL)
    {
        new->gang = NULL;
        new->next = new;
        new->n = 0;
        new->active = (1==1);
        lastORBITAL->next = new;
        lastORBITAL = new;
    }

    new->type = type;
    new->myatom = ap;
    new->ncenter = 1;
    new->couple[0] = new;
    new->ncouple = 1;
    new->a1 = NULL;
    new->a2 = NULL;
    new->a3 = NULL;
    new->a4 = NULL;
    new->a5 = NULL;
    if( i2 >= 0 ) new->a1 = a_m_serial(i2);
    if( i3 >= 0 ) new->a2 = a_m_serial(i3);
    if( i4 >= 0 ) new->a3 = a_m_serial(i4);
    if( i5 >= 0 ) new->a4 = a_m_serial(i5);
    if( i6 >= 0 ) new->a5 = a_m_serial(i6);
    if( type == Orm ){ /*  multi atom center */
        if( i2 >= 0 ) new->ncenter += 1;
        if( i3 >= 0 ) new->ncenter += 1;
        if( i4 >= 0 ) new->ncenter += 1;
        if( i5 >= 0 ) new->ncenter += 1;
        if( i6 >= 0 ) new->ncenter += 1;
    }
    new->osn = osn;
    new->along = along;
    new->along2 = x;
    new->x = x;
    new->y = y;
    new->z = z;
    new->spin = spin;
    new->ipair = pair;
    return 1;
}
int  expand(osn,n,a1,b1,a2,b2,a3,b3,a4,b4,a5,b5,a6,b6)
int osn,n;
float a1,a2,a3,a4,a5,a6;
float b1,b2,b3,b4,b5,b6;
{

    ORBITAL *op,*o_m_serial();

    op = o_m_serial(osn);
    if( op == NULL )
    { aaerror("no orbital for expansion\n"); return 0;}
    if( n > 6 ) n = 6;
    if( n < 0 ) n = 0;
    op->n = n;
    op->a[0] = a1;
    op->a[1] = a2;
    op->a[2] = a3;
    op->a[3] = a4;
    op->a[4] = a5;
    op->a[5] = a6;
    if( b1 > 0.)
        op->rl[0] = log(b1);
    else
        op->rl[0] = b1;
    if( b2 > 0.)
        op->rl[1] = log(b2);
    else
        op->rl[1] = b2;
    if( b3 > 0.)
        op->rl[2] = log(b3);
    else
        op->rl[2] = b3;
    if( b4 > 0.)
        op->rl[3] = log(b4);
    else
        op->rl[3] = b4;
    if( b4 > 0.)
        op->rl[4] = log(b4);
    else
        op->rl[4] = b4;


    return 1;
}
int dump_orbit(outp)
FILE *outp;
{
    ORBITAL *op,*op1,*o_next();
    int i,no,o_number(),natom,a_number();
    int j;
    int s1,s2,s3,s4,s5;
    int numatm,a_number();
    ATOM *ap,*a_next();;

    no = o_number();
    if( no <= 0) return; /*no error to dump nonextant*/
    for( i=0; i< no; i++)
    {
        op = o_next(i);
        switch( op->type )
        {
        case Or1:
            fprintf(outp,"orbit 1 %d %d %d %d ;\n",
                    (op->myatom)->serial, op->osn,op->spin,op->ipair);
            break;
        case Or1o:
            fprintf(outp,"orbit 1o %d %d %f %f %f %f %d %d;\n",
                    (op->myatom)->serial,op->osn,op->along,op->x,op->y,op->z,op->spin,op->ipair);
            break;
        case Or2:
            fprintf(outp,"orbit 2 %d %d %d %f %d %d;\n",
                    (op->myatom)->serial,
                    (op->a1)->serial,
                    op->osn,op->along,op->spin,op->ipair);
            break;
        case Or2p:
            fprintf(outp,"orbit 2p %d %d %d %f %d %d;\n",
                    (op->myatom)->serial,
                    (op->a1)->serial,
                    op->osn,op->along,op->spin,op->ipair);
            break;
        case Or3:
            fprintf(outp,"orbit 3 %d %d %d %d %f %f %d %d;\n",
                    (op->myatom)->serial,
                    (op->a1)->serial,
                    (op->a2)->serial,
                    op->osn,op->along,op->along2,op->spin,op->ipair);
            break;
        case Or4s:
            fprintf(outp,"orbit 4s %d %d %d %d %d %f %d %d;\n",
                    (op->myatom)->serial,
                    (op->a1)->serial,
                    (op->a2)->serial,
                    (op->a3)->serial,
                    op->osn,op->along,op->spin,op->ipair);
            break;
        case Or4p:
            fprintf(outp,"orbit 4p %d %d %d %d %d %f %d %d;\n",
                    (op->myatom)->serial,
                    (op->a1)->serial,
                    (op->a2)->serial,
                    (op->a3)->serial,
                    op->osn,op->along,op->spin,op->ipair);
            break;
        case Orm:
            s1 = -1;
            s2 = -1;
            s3 = -1;
            s4 = -1;
            s5 = -1;
            if( op->a1  != NULL ) s1 = (op->a1)->serial;
            if( op->a2  != NULL ) s2 = (op->a2)->serial;
            if( op->a3  != NULL ) s3 = (op->a3)->serial;
            if( op->a4  != NULL ) s4 = (op->a4)->serial;
            if( op->a5  != NULL ) s5 = (op->a5)->serial;
            fprintf(outp,"orbit m %d %d %d %d %d %d  %d %d %d;\n",
                    (op->myatom)->serial,
                    s1,s2,s3,s4,s5,
                    op->osn,op->spin,op->ipair);
            break;
        }/* end switch */
    }

    for( i=0; i< no; i++)
    {
        op = o_next(i);
        fprintf(outp,"expand %d %d ",op->osn,op->n);
        for( j=0; j< op->n; j++)
            if( op->rl[j] > 0.)
                fprintf(outp,"%g %g ",op->a[j],exp(op->rl[j]));
            else
                fprintf(outp,"%g %g ",op->a[j],(op->rl[j]));

        fprintf(outp,";\n");
    }

    natom = a_number();
    for( i=0; i< natom; i++)
    {
        ap =a_next(i) ;
        if( ap->na > 0 )
            fprintf(outp,"mov %d.na %f;\n",ap->serial,ap->na);
    }

    for( i=0; i< no; i++)
    {
        op = o_next(i);
        if( op->gang != NULL )
        {
            op1 = op->gang;
            fprintf(outp,"dscf gang %d %d ;\n", op->osn,op1->osn);
        }
        if( !op-> active )
            fprintf(outp,"dscf freeze %d ;\n",op->osn);

        if( op->ncouple > 1){
            for( j=1; j< op->ncouple; j++)
            {
                op1 = op->couple[j];
                fprintf(outp,"dscf couple %d %d;\n",op->osn,op1->osn);
            }
        }


    }

}
float Loose_Ritz_product(odummy)
ORBITAL *odummy;
{
    float Basic_Ritz_product();
    float Ritz_orth();
    float H,x;
    int no,o_number();
    /*	H = Basic_Ritz_product(odummy);
    */
    H = Ritz_orth();
    x = fabs( H +total_kinet/total_norm);
    no = o_number();
    return H;
    /*
    	return H + x*2/no;
    */
}

float Just_Ritz_product(odummy)
ORBITAL *odummy;
{
    float Ritz_orth();
    float Basic_Ritz_product();
    float H,x;
    /*
    	H = Basic_Ritz_product(odummy);
    */
    H = Ritz_orth();
    return H ;
}

float orbit_ktarget;
float set_Ritz_target(odummy)
ORBITAL *odummy;
{
    float Basic_Ritz_product();
    float Ritz_orth();
    float H,x;
    int no,o_number();
    return 0.;
    H = Ritz_orth();
    x = fabs( H +total_kinet/total_norm);
    orbit_ktarget = total_kinet/total_norm;
    return H + x;
}

float Ritz_product(odummy)
ORBITAL *odummy;
{
    float Basic_Ritz_product();
    float Ritz_orth();
    float H,x;
    int no,o_number();
    /*
    	H = Basic_Ritz_product(odummy);
    */
    H = Ritz_orth();
    x = fabs( H +total_kinet/total_norm);
    return H + x;
}
float INDO_Ritz_product(odummy)
ORBITAL *odummy;
{
    float Basic_INDO_Ritz_product();
    float Ritz_orth_INDO();
    float H,x;
    int no,o_number();
    /*
    	H = Basic_INDO_Ritz_product(odummy);
    */
    /*	x = fabs( H +total_kinet/total_norm);
    */
    H = Ritz_orth_INDO(0.0001);
    return H;
    /*	return H + x;
    */
}
float tight_INDO_Ritz_product(odummy)
ORBITAL *odummy;
{
    float Basic_INDO_Ritz_product();
    float Ritz_orth_INDO();
    float H,x;
    int no,o_number();
    /*
    	H = Basic_INDO_Ritz_product(odummy);
    */
    H = Ritz_orth_INDO(0.0001);
    x = fabs( H +total_kinet/total_norm);
    return H + x;
}


int renormalize()
{
    ORBITAL *op,*qp,*o_next();
    int i,j,k,no,o_number();
    float mag,phiphi();

    no = o_number();
    if( no < 1) return;
    op = o_first();
    for( i=0; i< no; i++)
    {
        if( op->ncouple > 0){
            mag = phiphi(op,op);
            mag = 1./sqrt(mag);
            for( k=0; k<op->ncouple; k++)
            {
                qp = op->couple[k];
                for(j=0; j< qp->n; j++)
                    qp->a[j] *= mag;
            } /* k */
        } /* op->ncouple > 0 */
        op = op->next;
    }
}

int report(out)
FILE *out;
{
    fprintf(out," sums \n");
    fprintf(out," nuclear    %f\n",total_nuc);
    fprintf(out," kinetic    %f\n",total_kinet/total_norm);
    fprintf(out," e-/nuclear %f\n",total_enuc/total_norm);
    fprintf(out," coloumb    %f\n",total_col/total_norm);
    fprintf(out," exchange   %f\n",total_ex/total_norm);
    fprintf(out," Normalizer %f\n",total_norm);
    fprintf(out," Energy H   %f\n", total_nuc+
            (total_kinet+total_enuc+total_col+total_ex)/total_norm);
    fprintf(out," Virial error  %f\n", fabs(total_nuc+
                                            (total_kinet+total_kinet+total_enuc+total_col+total_ex)/total_norm));
    return;
}



float phiphi(oo1,oo2)
ORBITAL *oo1,*oo2;
{
    int i,j,k,l;
    int kk,ll;
    float x,y,z,r;
    float rab,rs,re;
    float overlap;
    float s1,s2,sign;
    ORBITAL *o1,*o2;
    overlap = 0.;

    if( oo1->ncouple < 1 && oo1 == oo2) return 1.;
    if( oo1->ncouple < 1) return 0.;
    if( oo2->ncouple < 1) return 0.;
    for( kk=0; kk< oo1->ncouple; kk++)
    {
        o1 = oo1->couple[kk];
        s1 = 1.;
        for( k = 0; k < o1->ncenter; k++)
        {
            if( o1->type == Or4p) s1 = -s1;
            if( o1->type == Or2p) s1 = -s1;
            s2 = 1.;
            for( ll = 0; ll < oo2->ncouple; ll++)
            {
                o2 = oo2->couple[ll];
                for( l = 0; l < o2->ncenter; l++)
                {
                    if( o2->type == Or4p) s2 = -s2;
                    if( o2->type == Or2p) s2 = -s2;
                    sign = s1*s2;
                    x = o2->rx[l] - o1->rx[k];
                    y = o2->ry[l] - o1->ry[k];
                    z = o2->rz[l] - o1->rz[k];
                    rab = (x*x + y*y + z*z)*INVBOHR*INVBOHR;
                    for( i=0; i< o1->n; i++)
                    {
                        for( j=0; j< o2->n; j++)
                        {
                            rs = (o1->r[i] + o2->r[j]);
                            /* not the error
                            		if( rs < 1.e-6) rs = 1.e-6;
                            */
                            r = o1->r[i]*o2->r[j]/rs;
                            r = r*rab;
                            if( r > 70. ) r = 70.;
                            if( r < -70.) r = -70.;
                            re = exp(-r);
                            overlap += sign*o1->a[i]*o2->a[j]*pow(PI/rs,1.5)*re;
                        }
                    }
                }} /* ll,l */
        }}/* kk,k */
    return overlap;
}


float nuclear()
{
    ATOM *a1,*a2,*a_next();
    int i,j,natom,a_number();
    int l;
    float x,y,z,r,r2,r6,accum;
    float q1,q2;
    natom = a_number();
    if( natom < 2) return 0.;

    accum = 0.;
    a1 = a_next(-1);
    a1 = a1->next;
    for( i=1; i< natom ; i++)
    {
        q1 = a1->na;
        /*
        		if( q1 == 0 ) q1 = a1->q;
        */
        if( q1 != 0 ){
            for( j= 0; j< i; j++)
            {
                a2 = a_next(j);
                q2 = a2->na;
                if( q2 == 0.  ){
                    q2 = a2->q;
                    x = a1->x - a2->x;
                    y = a1->y - a2->y;
                    z = a1->z - a2->z;
                    r = 1./sqrt(x*x+y*y+z*z);
                    for(l=0; l< a1->dontuse; l++)
                    {
                        if(a2 == a1->excluded[l]) r = 0.;
                    }
                    r2 = r*r;
                    r6 = r2*r2*r2;
                    r2 = r6*r6;
                    r = q1*q2*r*332.17752 - a1->a*a2->a*r6 + a1->b*a2->b*r2 ;
                    /* accum = accum + r / 627.51; */
                    accum = accum + r * 0.0015936001;
                } else {
                    x = a1->x - a2->x;
                    y = a1->y - a2->y;
                    z = a1->z - a2->z;
                    r = sqrt(x*x+y*y+z*z)*INVBOHR;
                    accum = accum + q1*q2/r;
                }
            }
        }/* q1 != 0 */
        a1 = a1->next;
    }
    return accum;
}

float Fzero( x)
float x;
{
#ifdef notanyoldSGI
    float accum,etox;
    accum = 0.;
    etox = exp(-x);
    /*
            accum = 1./(23.)*(2.*x*accum + etox);
            accum = 1./(21.)*(2.*x*accum + etox);
    */
    accum = 1./(19.)*(2.*x*accum + etox);
    accum = 1./(17.)*(2.*x*accum + etox);
    accum = 1./(15.)*(2.*x*accum + etox);
    accum = 1./(13.)*(2.*x*accum + etox);
    accum = 1./(11.)*(2.*x*accum + etox);
    accum = 1./( 9.)*(2.*x*accum + etox);
    accum = 1./( 7.)*(2.*x*accum + etox);
    accum = 1./( 5.)*(2.*x*accum + etox);
    accum = 1./( 3.)*(2.*x*accum + etox);
    accum = 2.*x*accum + etox;
    return accum;
#else
    if( x < 1.e-7) return 1.;
    x = sqrt(x);
    return .5*ROOTPI/x*erf(x);

#endif

}/* end of Fzero */

/* Fone(x) is -d/dx(Fzero(x)) */
float Fone(x)
float x;
{
    float y ;
    if( x < 1.e-7) return 0.;
    y = sqrt(x);
    y = .5*ROOTPI/y*erf(y);   /* Fzero */
    /* Fone is defined by the recursion
    *  Fzero = (2xFone + exp(-x)) 
    * or 
    *  Fone = 1./(2x)*(Fzero-exp(-x))
    */
    return (y-exp(-x))/(x+x);
} /* end of Fone */


/* go through the list and force each normal to be the current one */
/* since i can't imagine just doing this also do the orbital origin */
int o_update_normals()
{
    ORBITAL *op,*op1,*o_next();
    ATOM *ap;
    int i,no,o_number();
    int j;
    float x,y,z,r;
    float x1,y1,z1,r1;
    float x2,y2,z2;

    no = o_number();
    if( no < 1) return;
    op = o_first();
    for( i=0; i< no; i++)
    {
        /*	op = o_next(i);
        */
        if( op->gang != NULL && op->gang != op )
        { /* if i'm ganged */
            op1 = op->gang;
            while(op1->gang != NULL ){ op1 = op1->gang; op->gang = op1;}
            op->x = op1->x;
            op->y = op1->y;
            op->z = op1->z;
            op->along = op1->along;
            op->along2 = op1->along2;
            op->n = op1->n;
            for( j=0; j< op->n; j++)
            {
                op->rl[j] = op1->rl[j];
                op->a[j] = op1->a[j];
            }
        }

        for( j=0; j< op->n; j++)
            if( op->rl[j] < -6.) op->rl[j] = -6.;
        /*
        		if( op->rl[j] < -10.) op->rl[j] = -10.;
        		if( op->rl[j] < -9.) op->rl[j] = -9.;
        		if( op->rl[j] < -4.) op->rl[j] = -4.;
        		if( op->rl[j] < -3.) op->rl[j] = -3.;
        */
        for( j=0 ; j < op->n; j++)
            op->r[j] = exp( op->rl[j]);

        switch( op->type)
        {
        case Or1:
            ap = op->myatom;
            op->rx[0] = ap->x;
            op->ry[0] = ap->y;
            op->rz[0] = ap->z;
            break;

        case Or1o:
            ap = op->myatom;
            op->rx[0] = ap->x + op->x;
            op->ry[0] = ap->y + op->y;
            op->rz[0] = ap->z + op->z;
            break;

        case Or2:  /* direction along a bond */
        case Or2p:  /* direction along a bond */
            ap = op->a1;
            x = ap->x;
            y = ap->y;
            z = ap->z;
            ap = op->myatom;
            x -= ap->x;
            y -= ap->y;
            z -= ap->z;
            op->x = x;
            op->y = y;
            op->z = z;
            if( op->type == Or2){
                if( op->along < -.5 ) op->along = -.5;
                if( op->along > 2.0 ) op->along = 2.0;
            }else{
                if( op->along < 0. && op->along > -0.05) op->along = -0.05;
                if( op->along > 0. && op->along < 0.05) op->along  = 0.05;
            }
            op->rx[0] =   ap->x + op->along*op->x;
            op->ry[0] =   ap->y + op->along*op->y;
            op->rz[0] =   ap->z + op->along*op->z;
            if( op->type == Or2p)
            {
                op->ncenter = 2;
                op->rx[1] =   ap->x - op->along*op->x;
                op->ry[1] =   ap->y - op->along*op->y;
                op->rz[1] =   ap->z - op->along*op->z;
            }
            break;

        case Or3: /* along a angle diagonal  and the cross product*/
            ap = op->a1;
            x = ap->x;
            y = ap->y;
            z = ap->z;
            ap = op->a2;
            x1 = ap->x;
            y1 = ap->y;
            z1 = ap->z;
            ap = op->myatom;
            x -= ap->x;
            y -= ap->y;
            z -= ap->z;
            x1 -= ap->x;
            y1 -= ap->y;
            z1 -= ap->z;
            r = sqrt(x*x + y*y + z*z);
            r1 = sqrt(x1*x1 + y1*y1 + z1*z1);
            if( r < 1.e-7) r = 1.;
            if( r1 < 1.e-7) r1 = 1.;
            op->x = (x/r +x1/r1)*ROOTHALF;
            op->y = (x/r +x1/r1)*ROOTHALF;
            op->z = (x/r +x1/r1)*ROOTHALF;
            /*
            	if( op->along > .6) op->along = .6;
            	if( op->along < -.6) op->along = -.6;
            */
            op->rx[0] =   ap->x + op->along*op->x;
            op->ry[0] =   ap->y + op->along*op->y;
            op->rz[0] =   ap->z + op->along*op->z;
            /* the cross product */
            x2 = y*z1 - z*y1;
            y2 = z*x1 - x*z1;
            z2 = x*y1 - y*x1;
            r =  sqrt(x2*x2 + y2*y2 + z2*z2);
            /*
            	if( op->along2 > .4) op->along2 = .4;
            	if( op->along2 < -.4) op->along2 = -.4;
            */
            op->rx[0] += op->along2*x2/r;
            op->ry[0] += op->along2*y2/r;
            op->rz[0] += op->along2*z2/r;
            break;

        case Or4s: /* perp to a plane */
        case Or4p: /* perp to a plane */
            if( op->type == Or4p && op->along < .05)  op->along = .05;
            ap = op->a2;
            x = ap->x;
            y = ap->y;
            z = ap->z;
            ap = op->a3;
            x1 = ap->x;
            y1 = ap->y;
            z1 = ap->z;
            ap = op->a1;
            x -= ap->x;
            y -= ap->y;
            z -= ap->z;
            x1 -= ap->x;
            y1 -= ap->y;
            z1 -= ap->z;
            x2 = y*z1 - z*y1;
            y2 = z*x1 - x*z1;
            z2 = x*y1 - y*x1;
            r =  sqrt(x2*x2 + y2*y2 + z2*z2);
            if( r > 0.){
                op->x = x2/r;
                op->y = y2/r;
                op->z = z2/r;
            } else {
                op->x = 0.;
                op->y = 0.;
                op->z = 0.;
            }
            ap = op->myatom;
            op->rx[0] =   ap->x + op->along*op->x;
            op->ry[0] =   ap->y + op->along*op->y;
            op->rz[0] =   ap->z + op->along*op->z;
            if( op->type == Or4p)
            {
                op->ncenter = 2;
                op->rx[1] =   ap->x - op->along*op->x;
                op->ry[1] =   ap->y - op->along*op->y;
                op->rz[1] =   ap->z - op->along*op->z;
            }

            break;
        case Orm:
            op->rx[0] = (op->myatom)->x;
            op->ry[0] = (op->myatom)->y;
            op->rz[0] = (op->myatom)->z;
            if( op->ncenter >1 ){
                op->rx[1] = (op->a1)->x;
                op->ry[1] = (op->a1)->y;
                op->rz[1] = (op->a1)->z;
            }
            if( op->ncenter >2 ){
                op->rx[2] = (op->a2)->x;
                op->ry[2] = (op->a2)->y;
                op->rz[2] = (op->a2)->z;
            }
            if( op->ncenter >3 ){
                op->rx[3] = (op->a3)->x;
                op->ry[3] = (op->a3)->y;
                op->rz[3] = (op->a3)->z;
            }
            if( op->ncenter >4 ){
                op->rx[4] = (op->a4)->x;
                op->ry[4] = (op->a4)->y;
                op->rz[4] = (op->a4)->z;
            }
            if( op->ncenter >5 ){
                op->rx[5] = (op->a5)->x;
                op->ry[5] = (op->a5)->y;
                op->rz[5] = (op->a5)->z;
            }
            break;
        }
        op = op->next;
    }
    /*
    	renormalize();
    */

}

/* useful functions from atom.c */
/* function o_number()
* returns number of orbits defined
*  this is just orbitNUMBER if orbitUPDATE == 0
*  other wise just figure it out
*/
int o_number()
{
    ORBITAL *op;
    if( orbitUPDATE )
    {
        orbitUPDATE = 0;
        orbitNUMBER = 0;
        if( firstORBITAL == NULL ) return 0 ;
        op = firstORBITAL;
        while(1)
        {
            if( op->next == NULL) break;
            orbitNUMBER++;
            if( op->next == op ) break;
            op = op->next;
        }
    }
    return orbitNUMBER;
}


/* function o_m_serial( serial )
* returns NULL on error or returns the address of the ORBITAL
* which matches serial
* cute?
*/
ORBITAL *o_m_serial( serial )
int serial;
{
    static ORBITAL *op = NULL;
    static ORBITAL *lastmatched = NULL;
    int i , n, o_number();

    if( orbitUPDATE) n= o_number();
    else n = orbitNUMBER;

    op = firstORBITAL;
    /* static pointer is hook for more efficient search */
    if( op == NULL) return NULL;
    if( lastmatched == NULL ) lastmatched = firstORBITAL;

    if( serial == lastmatched->osn) return lastmatched;
    if( serial > lastmatched->osn) op = lastmatched;
    for( i=0; i< n; i++ )
    {
        if( op-> osn == serial) {lastmatched = op;return op;}
        if( op == op->next)op = firstORBITAL ;
        else op = op->next;
    }
    return NULL;
}


/* function o_next( flag )
* returns NULL on error or last orbital
* then steps to the next
* cute?
* flag <= 0 starts it off
*/
ORBITAL *o_next( flag )
int flag;
{
    static ORBITAL *op = NULL;
    if( op == NULL) op = firstORBITAL ;
    if( op == NULL) return NULL;
if( flag <= 0){ op = firstORBITAL; return op;}
    if( op == op->next) return NULL;
    op = op->next;
    return op;
}
float total_orbit_norm(void)
{
    ORBITAL *o1,*o2,*o_next();
    int no, i,j,o_number();
    float phiphi();
    float accum_over,accum_elec;
    float x;

    no = o_number();
    if( no == 0 ) return 1.;
    accum_over = 0.;
    accum_elec = 0.;
    o1 = o_first();
    for( i = 0 ; i< no; i++)
    {
        if( o1->ipair == 2) accum_elec += 1.;
        accum_elec += 1.;
        for( j=0; j<no; j++)
        {
            o2 = o_next(j);
            x = phiphi(o1,o2);
            if( o1->ipair ==2) accum_over += x;
            accum_over += x;
        }
        o1 = o1->next;
    }

    return accum_elec/accum_over;
}
/*
* overlap
*  
* report on the overlap between orbitals
*/
void overlap( fp)
FILE *fp;
{
    int i,j,no,o_number();
    ORBITAL *o1,*o2,*o_next();
    float phiphi();

    no = o_number();
    o1 = o_next(-1);
    for( i=0; i< no; i++)
    {
        o2 = o_next(-1);
        for( j=0; j< i+1; j++)
        {
            fprintf(fp,"%f ",phiphi(o1,o2));
            o2  = o2->next;
        }
        fprintf(fp,"\n");
        o1 = o1->next;
    }
}/* end of overlap */
