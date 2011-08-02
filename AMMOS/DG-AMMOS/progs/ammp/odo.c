
#include <math.h>
#include <string.h>
#include <stdlib.h>

/* this is a quicky to test and debug solve.c */
#include <stdio.h>
main ( )
{
    float matrix[4][4];
    float eig[4][4];
    float vector[4];
    int i,j;
    for( i=0; i< 4; i++)
        for( j=0; j< 4; j++)
            matrix[i][j] = 0.;
    for (i=0;i<4;i++)
    {
        matrix[i][i] = 2.;
        if( i < 3 )
        {matrix[i][i+1]=1.; }
        if( i > 1 )
        { matrix[i-1][i] = 1.;}
        vector[i] = 1.;  }

    for ( i=0;i<4;i++)
    { printf(" %f %f \n",matrix[i][i],vector[i]);}
    i = solve( &matrix,&vector,4,4);
    printf(" %d \n ",i);

    for ( i=0;i<4;i++)
    { printf(" %f %f \n",matrix[i][i],vector[i]);}
    for( i=0; i< 4; i++)
        for( j=0; j< 4; j++)
            matrix[i][j] = 0.;
    for (i=0;i<4;i++)
    {
        matrix[i][i] = 2.;
        if( i < 3 )
        {matrix[i][i+1]=1.; }
        if( i > 1 )
        { matrix[i-1][i] = 1.;}
        vector[i] = 1.;  }

    for ( i=0;i<4;i++)
    { printf(" %f %f \n",matrix[i][i],vector[i]);}

    i = odo_solve( &matrix,&eig,&vector,4,4);

    printf(" %d \n ",i);

    for ( i=0;i<4;i++)
    { printf(" %f %f \n",matrix[i][i],vector[i]);}
}


/* solve a poorly conditioned matrix with
*  orthogonal factorization 
* this is a sledgehammer approach to the solution of
* an otherwise simple problem.  
*/
int odo_solve( matrix,eig,vector,i1,i2)
int i1,i2;
float (*matrix)[],(*eig)[],(*vector)[];
{
    int jacobi();
    int i,j;
    float temp[10];

    /* form the Ot D O transform the eigenvectors (O) are returned in the
    *  row space of eig (i.e. (*eig)[ i*i1 + j ] */

    jacobi( matrix,eig,i1, 100*i1, -1.e-7);

    /* form D-1 O b */
    for( i=0; i< i1; i++)
        temp[i] = 0.;
    for( i = 0; i< i1; i++)
    {
        for( j=0; j< i1; j++)
            temp[j] += (*eig)[i*i1+j] * (*vector)[j];
    }
    for( i=0; i< i1; i++)
    {
        if( fabs( (*matrix)[i*i1 +i] ) > 1.e-7 )
        {
            temp[i] = temp[i]/(*matrix)[i*i1 + i];
        } else {
            /* svd-like treatment just carry the zero */
            printf("Warning odo_solve() has detected a singular value \n");
        }
    }
    for( i=0; i< i1; i++)
    { (*vector)[i] = 0.; }/* i*/
    /* now back convert Ot (D-1) O b */
    for( i=0; i< i1; i++)
    {
        for( j=0; j< i1; j++)
            (*vector)[i] += (*eig)[i*i1 + j]* temp[j];
    }/* i (2nd time ) */
}/* end of odo_solve() */



/*
* jacobi method for eigenvalue/eigenvector calculations
*/

int jacobi( am,em, n, maxit, toler)
float toler;
int n,maxit;
float (*am)[], (*em)[];
{

    float (*s1)[],(*s2)[];
    int iindex,jindex,i,j;
    int imax,jmax;
    int iter;
    float emax,r,sa,ca;

    s1 = malloc( n* sizeof(float));
    s2 = malloc( n* sizeof(float));
    if( s1 == NULL || s2 == NULL )
    { return 1; }

    /* set em (the eigenvector matrix)  to I */
    for( i=0; i< n; i++)
    {
        for( j=0; j< n; j++)
        {
            (*em)[i*n+j] = 0.;
        }
        (*em)[i*n+i] = 1.;
    }

    for( iter=0; iter< maxit; iter++)
    {

        emax = -1;
        imax = 0; jmax = 0;
        for( i=0; i< n; i++)
            for( j=i+1; j< n; j++)
            {
                if( fabs((*am)[i*n+j]) > emax)
                {emax = fabs((*am)[i*n+j]); imax = i; jmax = j;}
            }

        if( emax < toler)
        {  free(s1); free(s2); return 0; }

        r = (*am)[imax*n + imax] - (*am)[jmax*n+jmax];
        r = r*r + 4*(*am)[imax*n+jmax]*(*am)[imax*n+jmax];
        if( r <= 0. ) /* error return */
        { free(s1); free(s2); return 1; }
        r = sqrt(r);
        iindex = imax*n + imax;
        jindex = jmax*n + jmax;
        if( (*am)[iindex] > (*am)[jindex])
        {
            ca = .5*(1.+((*am)[iindex]-(*am)[jindex])/r);
            ca = sqrt(ca);
            if( (*am)[imax*n + jmax] < 0.) ca = -ca;
            sa = (*am)[imax*n+jmax]/r/ca;
        }else{
            sa = .5*(1.-((*am)[iindex]-(*am)[jindex])/r);
            sa = sqrt(sa);
            ca = (*am)[imax*n+jmax]/r/sa;
        }
        /* use the transformation */
        /* the rows */
        for( i=0; i< n; i++)
        {
            iindex = i*n;
            (*s1)[i] =  ca*(*am)[iindex +imax] + sa*(*am)[iindex + jmax];
            (*s2)[i] = -sa*(*am)[iindex +imax] + ca*(*am)[iindex + jmax];
        }
        for( i=0; i< n; i++)
        {
            iindex = i*n;
            (*am)[iindex  + imax] = (*s1)[i];
            (*am)[iindex  + jmax] = (*s2)[i];
        }
        /* the columns */
        for( i=0; i< n; i++)
        {
            iindex = imax*n;
            jindex = jmax*n;
            (*s1)[i] =  ca*(*am)[iindex +i] + sa*(*am)[jindex + i];
            (*s2)[i] = -sa*(*am)[iindex +i] + ca*(*am)[jindex + i];
        }
        for( i=0; i< n; i++)
        {
            iindex = imax*n;
            jindex = jmax*n;
            (*am)[iindex  + i] = (*s1)[i];
            (*am)[jindex  + i] = (*s2)[i];
        }
        /* and finally update v */
        for( i=0; i< n; i++)
        {
            iindex = i*n;
            (*s1)[i] =  ca*(*em)[iindex +imax] +sa*(*em)[iindex + jmax];
            (*s2)[i] = -sa*(*em)[iindex +imax] +ca*(*em)[iindex + jmax];
        }
        for( i=0;i< n; i++)
        {
            iindex = i*n;
            (*em)[iindex+imax] = (*s1)[i];
            (*em)[iindex+jmax] = (*s2)[i];
        }

    } /* end of iter loop */
    free(s1); free(s2);
    return 0;
}/*end of jacobi */


/*  this is a routine to solve a linear equation by
    guassian elimination.  (basically solve.for translated) */
/* in order to have the  array matrix be of any length it must be passed as
   a linear array.  Since C has the opposite convention for array packing from 
   FORTRAN ( row fastest rather than column fastest) the leading dimension
   ilead is the row size of the array to which matrix points */
solve( matrix,vector,irow,ilead )
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
