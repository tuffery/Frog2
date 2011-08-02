
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
