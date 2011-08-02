/*  simplex.c
*
*  polytope simplex minimizer
*
* AMMP version 
*  given a range of atom ID's and uses them to define hull
*  also uses the total potential, 
*
*
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


void simplex( toler,niter,var, potent, inpotent, imin, imax)
float var,toler;
int niter, imin,imax, inpotent;
int (*potent[])();
{
    float (*polytope)[];
    float (*fvals)[];
    float (*reflect)[];
    ATOM *(*ap)[],  *a_m_serial();
    int i,j,k,l,natom,nargs,mycase;
    int best, simplex_get_best();
    int worst, simplex_get_worst();
    void simplex_reflect();
    float simplex_get_var();
    float randf(); /* rng assumed to be seeded */
    float x;
    int (*ptem[2])(), (*ftem[2])();
    int f_bond(),v_bond();
    int f_angle(), v_angle();
    int v_nonbon(),u_v_nonbon();

    ptem[0] = v_bond; ftem[0] = f_bond;
    ptem[1] = v_angle; ftem[1] = f_angle;

    /* figure out how many atoms there are in the range given */
    if( imin > imax ){ i= imin; imin = imax ; imax = i;}
    natom = 0;
    for( i=imin; i <= imax; i ++)
    {
        if( a_m_serial(i) != NULL) natom ++;
    }
    if( natom == 0) return ;
    nargs = natom *3  ;


    polytope = malloc( nargs*(nargs + 1)*sizeof( float));
    if( polytope == NULL)
    { fprintf(stderr," cannot allocate memory in simplex.c \n");
        return ;}
    fvals = malloc( (nargs + 1)*sizeof( float));
    if( fvals == NULL)
    { fprintf(stderr," cannot allocate memory in simplex.c \n");
        return ;}
    reflect = malloc( (nargs )*sizeof( float));
    if( reflect == NULL)
    { fprintf(stderr,"reflect  cannot allocate memory in simplex.c \n");
        return ;}
    ap = malloc( (natom )*sizeof( ATOM *));
    if( ap == NULL)
    { fprintf(stderr," ap cannot allocate memory in simplex.c \n");
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

    for( i= 0 ; i< nargs+1 ; i++)
    {
        for( j = 0; j < nargs; j+=3)
        {
            k = j/3;
            (*polytope)[i*nargs + j] = (*ap)[k]->gx + var*(2*randf()-1);
            (*polytope)[i*nargs + j+1] = (*ap)[k]->gy + var*(2*randf()-1);
            (*polytope)[i*nargs + j+2] = (*ap)[k]->gz + var*(2*randf()-1);
            (*ap)[k] ->x = (*polytope)[ i*nargs + j ];
            (*ap)[k] ->y = (*polytope)[ i*nargs + j +1];
            (*ap)[k] ->z = (*polytope)[ i*nargs + j +2];
        }
        steep( ptem,ftem,2,5,0.);
        for( j=0; j< nargs ; j+=3)
        {
            k = j/3;
            (*polytope)[ i*nargs + j ] = (*ap)[k]->x;
            (*polytope)[ i*nargs + j+ 1 ] = (*ap)[k]->y;
            (*polytope)[ i*nargs + j+ 2 ] = (*ap)[k]->z;
        }
        (*fvals)[i] = 0.;
        for( k =0 ; k< inpotent; k++)
        {
            if( (*potent[k]) != v_nonbon && (*potent[k]) != u_v_nonbon)
                (*potent[k])(&(*fvals)[i],0.);
            else
                zone_nonbon( &(*fvals)[i],0., ap,natom);
        }

    }
    for( k=0; k< niter; k++ )
    {

        best = simplex_get_best( fvals,nargs);
        if( (var =simplex_get_var( fvals,nargs+1 )) < toler ) goto DONE;
        worst = simplex_get_worst( fvals,nargs  );

        printf(" %d best %d energy %f\n", k,best,(*fvals)[best]);
        printf(" %d worst %d energy %f\n", k,worst,(*fvals)[worst]);
        mycase = 0;

NEW_HULL:
        simplex_reflect( polytope, reflect , nargs, worst , 2.);
        for( j = 0; j < nargs; j+=3)
        {
            i = j/3;
            (*ap)[i] ->x = (*reflect)[  j ];
            (*ap)[i] ->y = (*reflect)[  j +1];
            (*ap)[i] ->z = (*reflect)[ j +2];
        }
EVALUATE:
        steep( ptem,ftem,2,10,0.);
        for( j=0; j< nargs ; j+=3)
        {
            i = j/3;
            (*reflect)[  j ] = (*ap)[i]->x;
            (*reflect)[  j+ 1 ] = (*ap)[i]->y;
            (*reflect)[ j+ 2 ] = (*ap)[i]->z;
        }
        x = 0.;
        for( j =0 ; j< inpotent; j++)
        {
            if( (*potent[j]) != v_nonbon && (*potent[j]) != u_v_nonbon)
                (*potent[j])(&x,0.);
            else
                zone_nonbon( &x,0., ap,natom);
        }
        if( x >= (*fvals)[worst] )
        {
            if( mycase == 0 ) { mycase = 1;
                simplex_reflect( polytope, reflect , nargs, worst , 1.);
            }
            if( mycase == 1 ) { mycase = 2;
                simplex_reflect( polytope, reflect , nargs, worst , .25);
            }
            if( mycase == 2 ) { mycase = 3;
                simplex_reflect( polytope, reflect , nargs, worst , .0);
            }
            if( mycase == 3 ) { mycase = 4;
                simplex_reflect( polytope, reflect , nargs, worst , -.25);
            }
            if( mycase == 4 ) { mycase = 5;
                simplex_reflect( polytope, reflect , nargs, worst , -.5);
            }
            if( mycase == 5 )
            { /* desparation !!! */
                simplex_contract( polytope, reflect, nargs,best,worst, .5);
                mycase = 0;
                for( i= 0 ; i< nargs+1 ; i++)
                {
                    for( j = 0; j < nargs; j+=3)
                    {
                        l = j/3;
                        (*ap)[l] ->x = (*polytope)[ i*nargs + j ];
                        (*ap)[l] ->y = (*polytope)[ i*nargs + j +1];
                        (*ap)[l] ->z = (*polytope)[ i*nargs + j +2];
                    }
                    steep( ptem,ftem,2,5,0.);
                    for( j=0; j< nargs ; j+=3)
                    {
                        l = j/3;
                        (*polytope)[ i*nargs + j ] = (*ap)[l]->x;
                        (*polytope)[ i*nargs + j+ 1 ] = (*ap)[l]->y;
                        (*polytope)[ i*nargs + j+ 2 ] = (*ap)[l]->z;
                    }
                    (*fvals)[i] = 0.;
                    for( l =0 ; l< inpotent; l++)
                    {
                        if( (*potent[l]) != v_nonbon && (*potent[l]) != u_v_nonbon)
                            (*potent[l])(&(*fvals)[i],0.);
                        else
                            zone_nonbon( &(*fvals)[i],0., ap,natom);
                    }
                }
                best = simplex_get_best( fvals,nargs );
                if( (var =simplex_get_var( fvals,nargs+1 )) < toler ) goto DONE;
                worst = simplex_get_worst( fvals,nargs );
                printf(" reflect %d best %d energy %f\n", k,best,(*fvals)[best]);
                printf(" reflect %d worst %d energy %f\n", k,worst,(*fvals)[worst]);
                printf(" refelct %f var %f toler\n", var,toler);
                goto NEW_HULL;
            }

            for( j = 0; j < nargs; j+=3)
            {
                i = j/3;
                (*ap)[i] ->x = (*reflect)[ j ];
                (*ap)[i] ->y = (*reflect)[ j + 1 ];
                (*ap)[i] ->z = (*reflect)[ j + 2 ];
            }

            goto EVALUATE;
        }
        (*fvals)[worst] = x;

        for( j = 0; j < nargs ; j++)
            (*polytope)[worst*nargs +j ] = (*reflect)[j];


    }

DONE:
    printf(" putting best %d into coordinates \n",best);
    fflush(stdout);
    for( i=0; i< nargs; i+=3)
    {
        j = i /3;
        (*ap)[j] ->x = (*polytope)[ best*nargs + i ];
        (*ap)[j] ->y = (*polytope)[ best*nargs + i +1];
        (*ap)[j] ->z = (*polytope)[ best*nargs + i +2];
    }
    free( ap );
    free( reflect );
    free( fvals );
    free( polytope );
}/* end of the routine */
/* simplex_get_best( fvals, number );
*  return the array index in fvals of the best (most negative)
* value of f
*/
int simplex_get_best( fvals, n )
int n;
float (*fvals)[];
{
    int i ;
    int who;
    float x;
    who = -1;
    x = 10e10;
    for( i=0; i< n+1; i++)
    {
        if( x > (*fvals)[i]) {
            who = i; x = (*fvals)[i];
        }
    }
    if( who < 0 ) who = 0;
    return who;
}
/* simplex_get_worst( fvals, number );
*  return the array index in fvals of the worst (least negative)
* value of f
*/
int simplex_get_worst( fvals, n )
int n;
float (*fvals)[];
{
    int i ;
    int who;
    float x;
    who = -1;
    x = -10e10;
    for( i=0; i< n+1; i++)
    {
        if( x < (*fvals)[i]) {
            who = i; x = (*fvals)[i];
        }
    }
    if( who < 0 ) who = 0;
    return who;
}

float simplex_get_var( fvals,n )
int n;
float (*fvals)[];
{
    int i;
    float sx, sx2;

    sx = 0.; sx2 =0.;
    for( i=0; i< n; i++)
    {
        sx += (*fvals)[i];
        sx2 += (*fvals)[i]*(*fvals)[i];
    }
    sx /= n;
    sx2 /= n;
    sx2 = sx2 - sx*sx;

    return sx2;

}
/*	simplex_reflect( polytope, reflect , nargs, worst , -1.);
*  reflect the worst around the mean 
*/
void simplex_reflect( pt, rf, n, worst, wait)
float (*pt)[],(*rf)[],wait;
int n, worst;
{

    int i,j ;

    for( j=0 ; j < n; j++)
        (*rf)[j] = 0.;
    for( i=0; i< worst; i++)
    {
        for( j=0 ; j< n; j++)
            (*rf)[j] += (*pt)[ i*n + j];
    }
    for( i=worst + 1; i< n+1; i++)
    {
        for( j=0 ; j< n; j++)
            (*rf)[j] += (*pt)[ i*n + j];
    }

    for( j = 0; j < n; j++)
    {
        (*rf)[j] /= n;
        (*rf)[j] =
            (*rf)[j] +wait*((*rf)[j] - (*pt)[ worst*n  + j]);
    }
}
simplex_contract( pt,rf, n,best,worst , howmuch )
int n,best,worst;
float howmuch;
float (*pt)[], (*rf)[];
{
    int i ,j ;
    float r1,r2,r3,randf();

    for( i=0; i< n+1; i++)
    {
        /*
        r1 = randf(); r2 = randf() +1;
        r3 = 1./(r1 + r2);
        */
        r1 = 2.; r2 = 1.; r3 = 1./(r1+r2);
        for( j = 0; j < n; j ++ )
            (*pt)[i*n + j] =r3*( r1*(*pt)[i*n +j] + r2*(*pt)[n*best +j])  ;
    }
}
/* simplex/quarternion rigid body solver */
void rigid( toler,niter,var, potent, inpotent, imin, imax)
float var,toler;
int niter, imin,imax, inpotent;
int (*potent[])();
{
    float polytope[56];
    float fvals[8];
    float reflect[8];
    ATOM *(*ap)[],  *a_m_serial();
    int i,j,k,l,natom,nargs,mycase;
    int best, simplex_get_best();
    int worst, simplex_get_worst();
    void simplex_reflect();
    void quarternion_rot_tran();
    float simplex_get_var();
    float randf(); /* rng assumed to be seeded */
    float x;
    int v_nonbon(),u_v_nonbon();

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
    { fprintf(stderr," ap cannot allocate memory in simplex.c \n");
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

    /* initialize the hull */
    for( i=0; i< 8; i++)
    {
        for( j=0; j< 7; j++)
        {
            polytope[ i*nargs + j ] = var*(2*randf()-1.);
        }
        /*
        		polytope[i*nargs + 3 ] *= .1;
        		polytope[i*nargs + 4 ] *= .1;
        		polytope[i*nargs + 5 ] *= .1;
        		polytope[i*nargs + 6 ] *= .1;
        */
        polytope[i*nargs + 3 ] += 1./var; /* I = 1 0 0 0*/
        quarternion_rot_tran( &polytope[i*nargs],ap,natom);
        fvals[i] = 0.;
        for( j=0; j < inpotent; j++)
        {
            if( (*potent[j]) != v_nonbon && (*potent[j]) != u_v_nonbon)
                (*potent[j])(&fvals[i],0.);
            else
                zone_nonbon(&fvals[i],0., ap, natom);
        }

    }

    for( k=0; k< niter; k++ )
    {

        best = simplex_get_best( fvals,nargs );
        if( (var =simplex_get_var( fvals,nargs+1 )) < toler ) goto DONE;
        worst = simplex_get_worst( fvals,nargs );

        printf(" %d best %d energy %f\n", k,best,fvals[best]);
        printf(" %d worst %d energy %f\n", k,worst,fvals[worst]);
        mycase = 0;

NEW_HULL:
        simplex_reflect( polytope, reflect , nargs, worst , 2.);
        quarternion_rot_tran( &reflect[0],ap,natom);
EVALUATE:
        x = 0.;
        for( j =0 ; j< inpotent; j++)
        {
            if( (*potent[j]) != v_nonbon && (*potent[j]) != u_v_nonbon)
                (*potent[j])(&x,0.);
            else
                zone_nonbon(&x,0., ap, natom);
        }
        if( x >= fvals[worst] )
        {
            if( mycase == 0 ) { mycase = 1;
                simplex_reflect( polytope, reflect , nargs, worst , 1.);
            }
            if( mycase == 1 ) { mycase = 2;
                simplex_reflect( polytope, reflect , nargs, worst , .25);
            }
            if( mycase == 2 ) { mycase = 3;
                simplex_reflect( polytope, reflect , nargs, worst , .0);
            }
            if( mycase == 3 ) { mycase = 4;
                simplex_reflect( polytope, reflect , nargs, worst , -.25);
            }
            if( mycase == 4 ) { mycase = 5;
                simplex_reflect( polytope, reflect , nargs, worst , -.5);
            }
            if( mycase == 5 )
            { /* desparation !!! */
                simplex_contract( polytope, reflect, nargs,best,worst, .5);
                mycase = 0;
                for( i= 0 ; i< nargs+1 ; i++)
                {

                    quarternion_rot_tran( &polytope[i*nargs],ap,natom);
                    fvals[i] = 0.;
                    for( l =0 ; l< inpotent; l++)
                    {
                        if( (*potent[l]) != v_nonbon && (*potent[l]) != u_v_nonbon)
                            (*potent[l])(&fvals[i],0.);
                        else
                            zone_nonbon(&fvals[i],0., ap, natom);
                    }
                }
                best = simplex_get_best( fvals,nargs );
                if( (var =simplex_get_var( fvals,nargs+1 )) < toler ) goto DONE;
                worst = simplex_get_worst( fvals,nargs );
                printf(" reflect %d best %d energy %f\n", k,best,fvals[best]);
                printf(" reflect %d worst %d energy %f\n", k,worst,fvals[worst]);
                printf(" reflect %f var %f toler\n", var,toler);
                goto NEW_HULL;
            }

            quarternion_rot_tran( &reflect[0],ap,natom);

            goto EVALUATE;
        }
        fvals[worst] = x;

        for( j = 0; j < nargs ; j++)
            polytope[worst*nargs +j ] = reflect[j];


    }

DONE:
    printf(" putting best %d into coordinates \n",best);
    fflush(stdout);
    quarternion_rot_tran( &polytope[best*nargs],ap,natom);
    free( ap );
}/* end of the routine */

/*		quarternion_rot_tran( &reflect[0],ap,natom);
*/
void quarternion_rot_tran( what,who,howmany)
float what[];
ATOM *(*who)[];
int howmany;
{

    int i;
    float norm;
    float rot[3][3];
    float x,y,z;
    float x1,x2,x3,x4,x5,x6;

    norm = what[3]*what[3] +   what [4]*  what [4]
           +what[5]*what[5]+  what [6]*  what [6];
    norm = sqrt(norm);
    if( norm == 0.) return;
    norm = 1./norm;
    what[3] *= norm ;
    what[4] *= norm ;
    what[5] *= norm ;
    what[6] *= norm ;
    x1 = what[3];
    x2 = what[4];
    x3 = what[5];
    x4 = what[6];
    rot[0][0] = x1*x1+x2*x2 - x3*x3-x4*x4;
    rot[1][0] = 2*(x2*x3+x1*x4);
    rot[2][0] = 2*(x2*x4-x1*x3);
    rot[0][1] = 2*(x2*x3-x1*x4);
    rot[1][1] = x1*x1-x2*x2+x3*x3-x4*x4;
    rot[2][1] = 2*(x3*x4+x1*x2);
    rot[0][2] = 2*(x2*x4+x1*x3);
    rot[1][2] = 2*(x3*x4-x1*x2);
    rot[2][2] = x1*x1-x2*x2-x3*x3+x4*x4;
    x4 = 0; x5 = 0; x6 = 0;
    for( i=0; i< howmany; i++)
    {
        x4 += (*who)[i]->gx;
        x5 += (*who)[i]->gy;
        x6 += (*who)[i]->gz;
    }
    x4 = x4/howmany;
    x5 = x5/howmany;
    x6 = x6/howmany;
    for (i = 0; i < howmany; i++)
    {
        x = what[0]+x4;
        y = what[1]+x5;
        z = what[2]+x6;
        x1 = (*who)[i]->gx-x4;
        x2 = (*who)[i]->gy-x5;
        x3 = (*who)[i]->gz-x6;
        x += rot[0][0]*x1 + rot[0][1]*x2 + rot[0][2]*x3;
        y += rot[1][0]*x1 + rot[1][1]*x2 + rot[1][2]*x3;
        z += rot[2][0]*x1 + rot[2][1]*x2 + rot[2][2]*x3;
        (*who)[i]->x = x;
        (*who)[i]->y = y;
        (*who)[i]->z = z;
    }

}/*end of routine */
