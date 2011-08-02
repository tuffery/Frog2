/* dgeom()
*
* standard distance geometry with the full eigen expansion
* of the whole big matrix.
* 
*  UUGH.
*  maybe if we can get this to work we can go to the power method
*
*/
/* experiments with deflation */
/*
*  copyright 1993,1994,1995,2001 Robert W. Harrison
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

/* this is based on an older iterative version that never
* quite worked.
* but to save time i've kept the call structure
* -this means it can be added to the other versions of AMMP without
* -a big headache
*/
void dgeom(echo,op, vfs,nfs,niter,origin,eigshift )
FILE *op;
int (*vfs[])(),nfs;
int niter,echo;
int origin;
float eigshift;
{
ATOM *ap,*bp,*a_next(),*aorig,*a_m_serial();
int numatm , a_number();
int  iter,i,j,k,l;
float x,y,z,xo,yo,zo,x1,y1,z1,r;
int v_bond(),v_angle(),u_v_nonbon(),v_noel();
int v_mmbond(),v_mmangle(),v_c_angle(),v_nonbon();
int v_hard(),v_abc();
int v_torsion(), v_step(), v_ttarget();

int jacobi(); /* in normal.c */
int a_v_zero();

int use_bond, use_angle, use_noel;
int use_torsion, use_ttarget, use_step;


float (*matrix)[];
float (*eig)[];
float (*raw)[];
float (*floyd)[]; /* keep the current floyded matrix */
float (*scratch)[]; /* keep the current floyded matrix */
float (*bat)[];/* deflation matrix (die fleidermaus) */
ATOM *(*allatom)[];

numatm = a_number();
if( numatm < 2 ) return; /* the solution for 1 atom is trivial (0,0,0) */

matrix = (float (*)[])malloc( numatm*numatm*sizeof(float));
raw = (float (*)[])malloc( numatm*numatm*sizeof(float));
eig = (float (*)[])malloc( numatm*numatm*sizeof(float));
floyd = (float (*)[])malloc( numatm*numatm*sizeof(float));
scratch = (float (*)[])malloc( numatm*numatm*sizeof(float));
bat = (float (*)[])malloc( numatm*numatm*sizeof(float));

allatom = (ATOM *(*)[])malloc( numatm *sizeof(ATOM *));


for( i=0; i< numatm; i++)
{ (*allatom)[i] = a_next(i);}


use_bond = (1==0);
use_angle = (1==0);
use_noel = (1==0);
use_step = (1==0);
use_torsion = (1==0);
use_ttarget = (1==0);
for( i=0; i< nfs; i++)
{
	if( vfs[i] == v_bond) use_bond = (1==1);
	if( vfs[i] == v_angle) use_angle = (1==1);
	if( vfs[i] == v_noel) use_noel = (1==1);
	if( vfs[i] == v_step) use_step = (1==1);
	if( vfs[i] == v_ttarget) use_ttarget = (1==1);
	if( vfs[i] == v_torsion) use_torsion = (1==1);
}

	for( i=0; i< numatm*numatm; i++)
		(*bat)[i] = 0.;

	iter = 0;
	for( i=0; i< numatm; i++)
	{
	a_v_zero();
	if( use_bond) gsdg_bond((*allatom)[i]);
	if( use_angle) gsdg_angle((*allatom)[i]);
//	if( use_torsion) gsdg_torsion((*allatom)[i]);
//	if( use_ttarget) gsdg_ttarget((*allatom)[i]);
	if( use_noel) gsdg_noel((*allatom)[i]);
	if( use_step) gsdg_step((*allatom)[i]);

	(*eig)[i*numatm+i]  = 0.;
	(*raw)[i*numatm+i] = 0.;
	ap = (*allatom)[i];
	for( j=i+1; j< numatm; j++)
	{
		bp= (*allatom)[j];
/* if i and j are both inactive calculate the distance and use it */
/* otherwise use the vx value for the j'th atom (which was set by gsdg_) */
	
	if( !ap->active && !bp->active)
	{
		xo = ap->x -bp->x;
		yo = ap->y -bp->y;
		zo = ap->z -bp->z;
		r = xo*xo + yo*yo +zo*zo;
		if( bp->vx > 0.001)
			fprintf(op,"%d %d %f here %f theoretical\n",ap->serial, bp->serial, r, bp->vx);
		iter +=1;
	}else{
		r = bp->vx;
		if( r < 0.001 ) r = 10.e10;
		else iter += 1;
	}
	(*raw)[i*numatm + j] = r;
	(*eig)[i*numatm + j] = r;
	(*raw)[j*numatm + i] = r;
	(*eig)[j*numatm + i] = r;
	}/* j */
	}/* i */
	fprintf(op, "%d distances in %d unique %f fraction\n", iter,
		numatm*(numatm-1)/2,(float)iter/(numatm*(numatm-1)/2));

#define FLOYD
#ifdef FLOYD
/* now use Floyd's algorithm to fill in the matrix */

	l = 0;
REFLOYD:;
		for( k=0; k< numatm; k++)
		{
	for( i=0; i<numatm; i++)
	for( j=0; j<numatm; j++)
	{

/* if we don't know the distance */	
	if( (*raw)[i*numatm + j] > 10.e6)
	{
		r = (*eig)[i*numatm + j];
		x = (*eig)[i*numatm +k] + (*eig)[k*numatm + j];
		if( x < r ) r = x;
		(*eig)[i*numatm + j] = r;
		(*eig)[j*numatm + i] = r;

	}


	}
		}
	if( l< numatm)
	for( i=0; i<numatm*numatm; i++)
	{
		r = (*eig)[i];
		(*floyd)[i] = r;
		if( r > 10.e6)
		{
		goto REFLOYD;
		} 
	}
	
#else
	for( i=0; i< numatm*numatm; i++)
		if( (*raw)[i] > 10.e6 ){ (*eig)[i] = 1. ; (*floyd)[i] = 1.;}
#endif

/* this is now where we do the deflation steps */
/* initially bat is a zero matrix and we correct it by deflation */

	if( niter <1) niter = 1;
	for( iter=0; iter< niter; iter++)
	{
	
	for( i=0; i< numatm*numatm; i++)
		{
		if( (*raw)[i] > 10.e6)
		(*eig)[i] = (*floyd)[i] + (*bat)[i];
		else
		(*eig)[i] = (*floyd)[i];
		}

/* now this is not a well-set matrix
*  we need to apply a constraint that the center of mass is zero */
	r = 0.;
	x = 0.;
	for( i=0; i< numatm; i++)
	for( j=0; j< numatm; j++)
	{
		   x1 = (*eig)[i*numatm +j];
		 {r += x1; x += 1.;}
	}
	r /= (x);

	for( i=0; i< numatm; i++)
	for( j=0; j< numatm; j++)
	{
		xo = 0.;
		yo = 0.;
		x = 0.;
		y = 0.;
		for( k=0; k< numatm; k++)
		{ 
		   x1 = (*eig)[i*numatm +k];
		   y1 = (*eig)[k*numatm +j];
		 { xo += x1; x += 1.;}
		 { yo += y1; y += 1.;}
		}
		if( x > 0.) xo /= x;
		if( y > 0.) yo /= y;
		(*matrix)[i*numatm + j] = -0.5*((*eig)[i*numatm+j] -xo - yo+r);
	}

/*
	for( i=0; i< numatm; i++)
	for( j=0; j< numatm; j++)
	{  if( j == 0)fprintf(op,"\n");
	 fprintf(op,"%f ", (*matrix)[i*numatm + j]); 
	}
	*/
/* finally we call jacobi and get the coordinates */
/*	jacobi(  matrix,eig, numatm, 10000, 0.0001);
 *	*/
	simultaneous_iteration( matrix, eig, scratch, NULL,(1==1), numatm,3,1000000,0.001);

if( echo)
{
	fprintf(op," the eigenvalues of the system \n");
	for( i=0; i< 3; i++)
	{ fprintf(op,"%f ", (*scratch)[i]); 
	  if( i%10 == 9) fprintf(op,"\n");}
	fprintf(op,"\n");
}
/* trace of square == sum of squares of eigenvalues */
	z = 0.;
	for( l=0; l< numatm; l++)
	{ 
		x = 0.;
		for( i=0; i< numatm; i++)
		x += (*matrix)[l*numatm+i]*(*matrix)[l*numatm+i];
			 z+= x;}
	z -= (*scratch)[0]*(*scratch)[0];
	z -= (*scratch)[1]*(*scratch)[1];
	z -= (*scratch)[2]*(*scratch)[2];

	fprintf(op," the current dimensional error %f\n\n", z);
	fflush(op);

	x = sqrt( (*scratch)[0]);
	y = sqrt( (*scratch)[ 1]);
	z = sqrt( (*scratch)[ 2]);
	for( i=0; i< numatm; i++)
	{
		ap = (*allatom)[i];
		ap->x = (*eig)[i]*x;
		ap->y = (*eig)[numatm+i]*y;
		ap->z = (*eig)[numatm+numatm+i]*z;
	}


/* 
 * this deflation attempt did not work 
 * need to search the space of (*eig) after scaling
	for( l=numatm-1; l < numatm; l++)
	{
	z = sqrt( fabs( (*matrix)[numatm*l+l]))*10;

	for( i=0; i< numatm; i++)
	for( j=0; j< numatm; j++)
	{
		(*bat)[i*numatm+j] -= z*(*eig)[numatm*l+i]*(*eig)[numatm*l+j];
	}
	}

*/
/*
int dgeom_deflate( floyd, bat,raw, eig,matrix, scratch, n, first,last)

*/
	if( iter < niter-1)
	dgeom_deflate(floyd, bat,raw,eig,matrix,scratch,numatm, 3,numatm);


	} /* iter */
/* clean up at the end of the routine */
free(allatom);
free(scratch);
free(bat);
free(floyd);
free(raw);
free(eig);
free(matrix);
}/* end of dgeom */


int dgeom_deflate( floyd, bat,raw, eig,matrix, scratch, n, first,last)
	float (*floyd)[];
	float (*bat)[];
	float (*raw)[];
	float (*matrix)[];
	float (*eig)[];
	float (*scratch)[];
	int n,first,last;
{
	float (*ls)[];
	float (*values)[];
	float (*guess)[];
	float r,x,y,z;	
	float xo,yo,zo;	
	float x1,y1,z1;	
	int iter,i,j,k,l,istep,changed;
	int iout,jout;
	int initial;
	float dstep, upper, middle,lower;
	float randf();
	float current_error;


	guess = (float (*)[])malloc(n*n *sizeof(float));
	values = (float (*)[])malloc(n *sizeof(float));
	ls = (float (*)[])malloc(n *sizeof(float));


	initial = (1==1);
	for( i=0; i< n*n; i++)
		(*guess)[i] = 0.;


	/* build the distance matrix in scratch */
		for( i=0; i< n*n; i++)
		{ if( (*raw)[i] < 10.e7) (*scratch)[i] = (*floyd)[i];
			else (*scratch)[i] = (*floyd)[i] + (*bat)[i];}

	r = 0.;
	x = 0.;
	for( i=0; i< n; i++)
	for( j=0; j< n; j++)
	{
		   x1 = (*scratch)[i*n +j];
		 {r += x1; x += 1.;}
	}
	r /= (x);

	for( i=0; i< n; i++)
	for( j=0; j< n; j++)
	{
		xo = 0.;
		yo = 0.;
		x = 0.;
		y = 0.;
		for( k=0; k< n; k++)
		{ 
		   x1 = (*scratch)[i*n +k];
		   y1 = (*scratch)[k*n +j];
		 { xo += x1; x += 1.;}
		 { yo += y1; y += 1.;}
		}
		if( x > 0.) xo /= x;
		if( y > 0.) yo /= y;
		(*matrix)[i*n + j] = -0.5*((*scratch)[i*n+j] -xo - yo+r);
	}
/* finally we call jacobi and get the coordinates */
	simultaneous_iteration( matrix, scratch, values, ls,initial, n,3,10000,0.001);
	initial = (1==0);
/* now calculate the cost of the first,...,last eigenvector products */
	z = 0.;

/* this is the trace of the square == sum of the squares */
	for( l=0; l<n; l++)
	{
		x = 0.;
		for( i=0; i< n; i++)
			x+= (*matrix)[l*n+i]*(*matrix)[l*n+i]; /* remember its a symmetric matrix*/

		z +=  x;
	}
	z -= (*values)[0]*(*values)[0];
	z -= (*values)[1]*(*values)[1];
	z -= (*values)[2]*(*values)[2];

	current_error = z;

	for( iter =0; iter< 100; ) /* iter will be incremented only when a hit*/
	{


		/*
		for( i=0; i< n-1; i++)
		for( j=i+1; j< n; j++)
		{
			if( (*raw)[i] < 10.e7)
			{ (*guess)[i*n+j] = 2.*randf()-1.; } 
			else 
			{ (*guess)[i*n+j] = 0.;}

			(*guess)[j*n+i] = (*guess)[i*n+j];
		}
		*/
		for( i=0; i< n*n; i++)
			(*guess)[i] = 0.;
		x = randf()*n;
		iout = x;
		x = randf()*n;
		jout = x;
		(*guess)[iout*n+jout] = 2.*randf()-1.;
		(*guess)[jout*n+iout] =  (*guess)[iout*n+jout];
	/* build the distance matrix in scratch */
		for( i=0; i< n*n; i++)
		{ if( (*raw)[i] < 10.e7) (*scratch)[i] = (*floyd)[i];
			else (*scratch)[i] = (*floyd)[i] + (*bat)[i] +(*guess)[i];}

	r = 0.;
	x = 0.;
	for( i=0; i< n; i++)
	for( j=0; j< n; j++)
	{
		   x1 = (*scratch)[i*n +j];
		 {r += x1; x += 1.;}
	}
	r /= (x);

	for( i=0; i< n; i++)
	for( j=0; j< n; j++)
	{
		xo = 0.;
		yo = 0.;
		x = 0.;
		y = 0.;
		for( k=0; k< n; k++)
		{ 
		   x1 = (*scratch)[i*n +k];
		   y1 = (*scratch)[k*n +j];
		 { xo += x1; x += 1.;}
		 { yo += y1; y += 1.;}
		}
		if( x > 0.) xo /= x;
		if( y > 0.) yo /= y;
		(*matrix)[i*n + j] = -0.5*((*scratch)[i*n+j] -xo - yo+r);

	}
/* finally we call jacobi and get the coordinates */
	simultaneous_iteration( matrix, scratch, values, ls,initial, n,3,10000,0.001);
	initial = (1==0);
/* now calculate the cost of the first,...,last eigenvector products */
	z = 0.;

/* this is the trace of the square == sum of the squares */
	for( l=0; l<n; l++)
	{
		x = 0.;
		for( i=0; i< n; i++)
			x+= (*matrix)[l*n+i]*(*matrix)[l*n+i]; /* remember its a symmetric matrix*/

		z +=  x;
	}
	z -= (*values)[0]*(*values)[0];
	z -= (*values)[1]*(*values)[1];
	z -= (*values)[2]*(*values)[2];
/* for sparse systems lets try a geometry measure */
	/*
	r = 0.;
	xo = (*values)[0];
	yo = (*values)[1];
	zo = (*values)[2];
	for( l=0; l < n-1; l++)
	{
		for( i=l; i< n; i++)
		{
			if( (*raw)[l*n+i] > 10.e6) continue;
		x1   = (*scratch)[l] -(*scratch)[i];
		y1   = (*scratch)[n+l] -(*scratch)[n+i];
		z1   = (*scratch)[n+n+l] -(*scratch)[n+n+i];
		x1 = xo*x1*x1 + yo*y1*y1 + zo*z1*z1;

		 x1 = ( x1 - (*floyd)[l*n+i]);
		 r += x1*x1;
		}	

	}
	z += r/n;
	*/
	/*
	printf(" %f %f \n", z,current_error);
	*/
	if( z < current_error)
	{
			iter += 1;
			current_error = z;
			for( i=0; i< n*n; i++)
				(*bat)[i] += (*guess)[i];
	}
		
	}/* iter */





	free(ls);
	free(values);
	free(guess);
}/* end of dgeom_deflate*/

int simultaneous_iteration( matrix, eig, values,scratch,initialize, n,nfit, niter, toler)
	float (*matrix)[], (*eig)[], (*values)[], (*scratch)[], toler;
	int n,nfit,niter;
{
	int iter,i,j,k,free_scratch;
	float deig,x,y;
	float randf();

		free_scratch = (1==0);
		if( scratch  == NULL)
		{
			free_scratch = (1==1);
			scratch = ( float (*)[])malloc( n* sizeof(float));
		}		
		if( initialize)
		{
			for( i=0; i< nfit*n; i++)
				(*eig)[i] = randf();
		}
			for( i=0; i< nfit; i++)
				(*values)[i] = 0.; /* this ensures that it iterate at least once */
		
		/* normalize */
		for( i=0; i< nfit; i++)
		{
			x = 0.;
			for( j=0; j< n; j++)
				x += (*eig)[i*n+j] * (*eig)[i*n+j];
			if( x < 1.e-7) aaerror("suspect zero eigenvector in simultaneous_iteration input");
			x = 1./sqrt(x);
			for( j=0; j< n; j++)
				(*eig)[i*n+j] *= x;

		}
		for( iter=0; iter < niter; iter++)
		{

			/*first find the products */
			for( k=0; k< nfit; k++)
			{
			for( i=0; i< n; i++)
				{
				(*scratch)[i] = 0.;
				for( j=0; j< n; j++)
				 (*scratch)[i] += (*matrix)[i*n+j]*(*eig)[k*n+j];
				} /*  i */

			for( i=0; i< n; i++)
				(*eig)[k*n+i] = (*scratch)[i];
			}/* k */
			/* then ortho-normalize them with Graham-schmidt orthogonalization */
			deig = 0.;
			for( i=0; i< nfit; i++)
			{
				x = 0.;
				for( j=0; j< n; j++)
					x += (*eig)[i*n+j] * (*eig)[i*n+j];
				x = sqrt(x);
				deig += fabs( x - (*values)[i]);	
				(*values)[i] = x;
				x = 1./x;
				for( j=0; j< n; j++)
					(*eig)[i*n+j] *= x;
			}/* i */
			/* now normal, but not orthogonal */
			for( k=0; k< nfit-1; k++)
			for( i=k+1; i < nfit; i++)
			{
				/* dot product */
				x = 0.;
				for( j=0; j< n; j++)
					x += (*eig)[k*n+j] * (*eig)[i*n+j];
				/* subtract */
				for( j=0; j< n; j++)
					(*eig)[i*n+j] -= x* (*eig)[k*n+j];
				/* re-normalize */
				x = 0.;
				for( j=0; j< n; j++)
					x += (*eig)[i*n+j] * (*eig)[i*n+j];
				x = 1./sqrt(x);
				for( j=0; j< n; j++)
					(*eig)[i*n+j] *= x;
			}/* k i */
			if( deig < toler ) break;

		}/* iter */	


		if( free_scratch) free(scratch);
}/* end of simultaneous_iteration */
