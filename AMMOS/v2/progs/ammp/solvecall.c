/* this is a quicky to test and debug solve.c */
#include stdio
main ( )
{
    float matrix[4][4];
    float vector[4];
    int i,j;
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

    i = gsolve( &matrix,&vector,4,4);

    printf(" %d \n ",i);

    for ( i=0;i<4;i++)
    { printf(" %f %f \n",matrix[i][i],vector[i]);}
}
#include "duck:[rob.c]solve.c"
