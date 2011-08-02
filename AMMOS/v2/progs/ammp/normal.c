/*  normal.c
*
*   calculate the normal modes of a molecule
*   this is best for small  molecules
*   
*  use finite differences for d2V/dxdy
*  use jacobi method to solve the equations
*
*  report the spectrum and modes, but what to do with them?
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

#define ANSI 1
/* misc includes - ANSI and some are just to be safe */
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#ifdef ANSI
#include <stdlib.h>
#endif
#include "ammp.h"


#define STEP  1.e-6
#define DSTEP  one/(STEP+STEP)
/*
#define DSTEP  one/(STEP)
*/
int FDnormal(forces,nfs, echo,op, normal )
int (*forces[])();
int nfs,echo;
float normal;
FILE *op;
{
    int numatm,a_number();
    float (*norm)[];
    float (*vect)[];
    float x;
    int i,j,ii;
    int iindex,jindex;
    int iforce,iatom;
    ATOM *bp;
    ATOM *ap,*a_next();
    char resid[5],atid[5],*np;
    float x1,x2,x3; /* for checking orthogonal */
    /* don't if there are no potentials */
    if( nfs < 1) return ;
    /* get the space */
    numatm = a_number();
    if( numatm < 2) return ;
    i = 3*numatm;
    j = i*i;
    norm = malloc( j*sizeof(float));
    vect = malloc( j*sizeof(float));
    if( norm == NULL || vect == NULL )
    { aaerror("cannot allocate memory in FDnormal"); return ;}

    /*  now do the finite differences */
    bp = a_next(-1); /* bp loops over the atoms with iatom */
    /* bp cannot be done with a_next because ap will */
    a_d_zero();
    for( iatom = 0; iatom < numatm; iatom++)
    {
        a_f_zero();
        bp->dx = 1.;
        for( iforce = 0; iforce < nfs; iforce++)
            (*forces[iforce])(STEP);
        a_ftogx( 1.,0.);
        a_f_zero();
        for( iforce = 0; iforce < nfs; iforce++)
            (*forces[iforce])(-STEP);
        a_ftogx( -1.,1.);
        ap = a_next(-1);
        for( i=0; i< numatm; i++)
        {
            j =  (iatom*3  )*3*numatm + i*3;
            (*norm)[j  ] = -0.5* ap->gx *DSTEP;
            (*norm)[j+1] = -0.5* ap->gy *DSTEP;
            (*norm)[j+2] = -0.5* ap->gz *DSTEP;
            ap = ap->next;
        }
        bp->dx = 0.;
        a_f_zero();
        bp->dy = 1.;
        for( iforce = 0; iforce < nfs; iforce++)
            (*forces[iforce])(STEP);
        a_ftogx( 1.,0.);
        a_f_zero();
        for( iforce = 0; iforce < nfs; iforce++)
            (*forces[iforce])(-STEP);
        a_ftogx( -1.,1.);
        ap = a_next(-1);
        for( i=0; i< numatm; i++)
        {
            j =  (iatom*3 + 1  )*3*numatm + i*3;
            (*norm)[j  ] = -0.5* ap->gx *DSTEP;
            (*norm)[j+1] = -0.5* ap->gy *DSTEP;
            (*norm)[j+2] = -0.5* ap->gz *DSTEP;
            ap = ap->next;
        }
        bp->dy = 0.;
        a_f_zero();
        bp->dz = 1.;
        for( iforce = 0; iforce < nfs; iforce++)
            (*forces[iforce])(STEP);
        a_ftogx( 1.,0.);
        a_f_zero();
        for( iforce = 0; iforce < nfs; iforce++)
            (*forces[iforce])(-STEP);
        a_ftogx( -1.,1.);
        ap = a_next(-1);
        for( i=0; i< numatm; i++)
        {
            j =  (iatom*3 + 2 )*3*numatm  + i*3;
            (*norm)[j  ] = -0.5* ap->gx *DSTEP;
            (*norm)[j+1] = -0.5* ap->gy *DSTEP;
            (*norm)[j+2] = -0.5* ap->gz *DSTEP;
            ap = ap->next;
        }
        bp->dz = 0.;
        bp = bp->next;
    }
    /* now symmetrize the matrix */
    for( i=0; i< numatm*3; i++)
        for( j=i; j < numatm*3; j++)
        {
            iindex = i*3*numatm + j;
            jindex = j*3*numatm + i;
            x = (*norm)[iindex] +(*norm)[jindex];
            x *= .5;
            (*norm)[iindex] = x;
            (*norm)[jindex] = x;
        }
    /* force the diagonal to be the sum of the others */
    for( i=0; i< numatm*3; i++)
    {
        x = 0.;
        for( j=0; j < i; j++)
        {
            jindex = i*3*numatm + j;
            x += (*norm)[jindex] ;
        }
        for( j=i+1; j < numatm*3; j++)
        {
            jindex = i*3*numatm + j;
            x += (*norm)[jindex] ;
        }
        (*norm)[i*3*numatm +i] = -x;
    }
    /* now mass-weight it */
    bp = a_next(-1);
    for( i=0; i< numatm*3; i+=3)
    {
        for( j=0; j< numatm; j++)
        {
            ap = a_next(j);
            x = one/sqrt(ap->mass*bp->mass);
            iindex = i*3*numatm + j*3;
            (*norm)[iindex    ] *= x;
            (*norm)[iindex + 1] *= x;
            (*norm)[iindex + 2] *= x;
            iindex = i*3*numatm +  3*numatm  + j*3;
            (*norm)[iindex    ] *= x;
            (*norm)[iindex + 1] *= x;
            (*norm)[iindex + 2] *= x;
            iindex = i*3*numatm + 6*numatm + j*3;
            (*norm)[iindex    ] *= x;
            (*norm)[iindex + 1] *= x;
            (*norm)[iindex + 2] *= x;
        }
        bp = bp->next;
    }

    /*
    	if( echo )
    	{
    		fprintf(op,"The mass-weighted Force Matrix \n");
    		for( i=0; i< numatm*3; i++)
    		{
    		for( j=0; j< numatm*3; j++)
    		{
    			fprintf(op,"%f ",(*norm)[i*numatm*3 + j ]);
    		}
    		fprintf(op,"\n");
    		}
    	}
    */

    if( jacobi( norm,vect,numatm*3, 100*numatm*numatm, 1.e-10) != 0)
    {/* error condition */
        aaerror(" Jacobi in FDnormal returns an error ");
    }
    /* check orthogonal */
    /*
    	for( i=0; i< numatm; i++)
    	{
    	x1 =0; x2 = 0; x3 = 0;
    	for (j=0; j< numatm; j++)
    	{
    	x1 +=(*vect)[ j*3]*(*vect)[j*3];
    	x1 +=(*vect)[ j*3+1]*(*vect)[j*3+1];
    	x1 +=(*vect)[ j*3+2]*(*vect)[j*3+2];
    	x2 +=(*vect)[i*numatm*3 +j*3]*(*vect)[i*numatm*3+j*3];
    	x2 +=(*vect)[i*numatm*3 +j*3+1]*(*vect)[i*numatm*3+j*3+1];
    	x2 +=(*vect)[i*numatm*3 +j*3+2]*(*vect)[i*numatm*3+j*3+2];
    	x3 +=(*vect)[j*3]*(*vect)[i*numatm*3+j*3];
    	x3 +=(*vect)[j*3+1]*(*vect)[i*numatm*3+j*3+1];
    	x3 +=(*vect)[j*3+2]*(*vect)[i*numatm*3+j*3+2];
    	}
    	fprintf(op,"normality check %f %d %f >%f< should be zero\n",x1,i,x2,x3);
    	}
    */
    /* end check */
    if( echo)
    {
        fprintf(op,"The %d  Eigenvalues\n",3*numatm);
        for( i=0; i< numatm*3; i++)
        {
            fprintf(op,"%f kcal/A^2g", (*norm)[i*3*numatm + i]);
            if( (*norm)[i*3*numatm+i] > 0 )
            {
                fprintf(op," %f cm-1\n",
                        sqrt(4.184E26*(*norm)[i*3*numatm + i]*.5)/2.997924e10/3.14159265);
            }else{
                fprintf(op," ***\n");
            }
        }
    }

    if( echo && normal > 0. )
    {
        fprintf(op,"The Eigenvectors \n");
        for( i=0; i< numatm*3; i++)
        {
            if( (*norm)[i*numatm*3 + i] > 0 )
            {
                fprintf(op,"REMARK %f cm-1\n",
                        sqrt(4.184E26*(*norm)[i*3*numatm + i]*.5)/2.997924e10/3.14159265);

                for( j=0; j< numatm; j++)
                {
                    ap = a_next(j);
                    np = ap->name;
                    for( ii= 0; ii< 5; ii++)
                        { if( *np != '.')
                        { if( islower(*np)) {resid[ii] = toupper(*np);}
                            else {resid[ii] = *np;}
                        } else { resid[ii] = '\0'; break;}
                        if(*np == '\0') break;
                        np++;
                    }
                    np++;
                    for( ii= 0; ii< 5; ii++)
                        { if( *np != '.')
                        { if( islower(*np)) {atid[ii] = toupper(*np);}
                            else {atid[ii] = *np;}
                        } else { atid[ii] = '\0'; break;}
                        if(*np == '\0') break;
                        np++;
                    }

                    fprintf(op,
                            "ATOM  %5d %-4s%c%-3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                            j,atid,' ',resid,ap->serial/100+i+1,
                            ap->x+(*vect)[i*numatm*3 +j*3]*normal,
                            ap->y+(*vect)[i*numatm*3 +j*3+1]*normal,
                            ap->z+(*vect)[i*numatm*3 +j*3+2]*normal,1.,10.);

                }
            }
            else {
                /*
                		fprintf(op," ***** cm-1\n" );
                		for( j=0; j< numatm*3; j++)
                		{
                			fprintf(op,"%f ",(*vect)[i*numatm*3 + j ]);
                		}
                		fprintf(op,"\n");
                */
            }
        }
    }


    free( vect );
    free( norm );

}/* end of FDnormal */

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
