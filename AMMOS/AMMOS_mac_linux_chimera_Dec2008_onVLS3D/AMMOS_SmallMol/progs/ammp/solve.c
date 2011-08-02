
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
