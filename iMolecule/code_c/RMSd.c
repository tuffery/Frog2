/* ============================================================================
 *  This software is part of Frog, a chemo informatics class able to build 
 *  3D coordinates for small compounds
 *   Copyright (C) 2006-2007 P. Tuffery, B.O. Villoutreix, Th. Bohme Leite, D. Gomes, M. Miteva, J. Chomilier
 *
 *   Frog2 (C) 2009-2010 by P. Tuffery, M. Miteva, F. Guyon
 *
 *   Using this software, please cite:
 *       Frog2: Efficient 3D conformation ensemble generator for small compounds.
 *       Miteva MA, Guyon F, Tuffery P.
 *       Nucleic Acids Res. 2010 Jul;38(Web Server issue):W622-7. Epub 2010 May 5.
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *   ==========================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define EPS 1.e-4
#define min(a,b) (a<b?a:b)
#define max(a,b) (a<b?b:a)

#define DtFloat		double	        /* Le type flottant.	        */
typedef DtFloat         DtPoint3[3];

static double *random_vect(int n) {
  int i;
  double *v=(double *) calloc(n, sizeof(double));
  for (i=0; i<n; i++) 
    v[i]= (double)rand()/(RAND_MAX+1.0);
  return(v);
}

static double inner(double *x, double *y, int n) {
  int i;
  double sum;

  for (sum=0, i=0; i<n; sum+=x[i]*y[i],i++);
  return sum;
    
}

static void product(double *result, double *A, double *x, int n) {
  int i, j;
  double sum;

  for (i=0; i<n; i++) {
    sum=0;
    for (j=0; j<n; j++)
      sum+=A[i+n*j]*x[j];
    result[i]=sum;
  }
}

static int power(double *a, int n, int maxiter, double eps, double *v, double *w) {
  int niter,i;
  double *y;
  double sum, l, normy, d;
  y=random_vect(n);
  niter=0;
  do {
    normy=sqrt(inner(y,y,n));
    for (i=0; i<n; i++) w[i]=y[i]/normy;
    product(y, a, w, n);
    l=inner(w,y,n);
    niter++;
    for (sum=0,i=0; i<n; i++) {
      d=y[i]-l*w[i];
      sum+=d*d;
    }
    d=sqrt(sum);
  } while (d>eps*fabs(l) && niter<maxiter);
  free(y);
  *v=l;
  return niter;
}

static double best_shift(double *a, int n) {
  double m, M, s;
  double t, sum;
  int i, j;
  t=a[0];
  for (i=1; i<n; i++) t=max(t, a[i+n*i]);
  M=t;
  t=a[0];
  for (i=0; i<n; i++) {
    for (sum=0,j=0; j<n; j++)
      if (j!=i) sum+=fabs(a[i+n*j]);
    t=min(t, a[i+n*i]-sum);
  }
  m=t;
  s=-0.5*(M+m);
  for (i=0; i<n; i++) 
    a[i+n*i]=a[i+n*i]+s;
  return s;
}


static int shift_power(double *a, int n, int maxiter, double eps, double *v, double *w) {  
  double sh;
  int niter;
  sh=best_shift(a, n);   

  niter=power(a, n, maxiter, eps, v, w);
  // fprintf(stderr,"niter=%d\n", niter);
  *v=*v-sh;
  return niter;
}

/*
 * calculs des coefficients d'Euler
 *
 */
static void euler(DtPoint3 *X,  DtPoint3 *Y, int n, double *v, double *w) {
  double Xm[3], Ym[3], K[3][3], A[16];
  double x, y, sum;
  int i,j,k;


  for (i=0; i<3;i++) {
    sum=0;
    for (j=0; j<n;j++)
      sum+=X[j][i];
    Xm[i]=sum/n;
  }
  for (i=0; i<3;i++) {
    sum=0;
    for (j=0; j<n;j++)
      sum+=Y[j][i]; 
    Ym[i]=sum/n;
  }

  for (i=0; i<3;i++)
    for (j=0; j<3; j++) {
      sum=0;
      for (k=0; k<n;k++)
	sum+=(Y[k][i]-Ym[i])*(X[k][j]-Xm[j]);
      K[i][j]=sum;
    }

  A[0+4*0]=K[0][0]+K[1][1]+K[2][2];
  A[1+4*1]=K[0][0]-K[1][1]-K[2][2];
  A[2+4*2]=-K[0][0]+K[1][1]-K[2][2];
  A[3+4*3]=-K[0][0]-K[1][1]+K[2][2];
  A[0+4*1]=A[1+4*0]=K[1][2]-K[2][1];
  A[0+4*2]=A[2+4*0]=K[2][0]-K[0][2];
  A[0+4*3]=A[3+4*0]=K[0][1]-K[1][0];
  A[1+4*2]=A[2+4*1]=K[0][1]+K[1][0];
  A[1+4*3]=A[3+4*1]=K[0][2]+K[2][0];
  A[2+4*3]=A[3+4*2]=K[1][2]+K[2][1];
  shift_power(A, 4, 10000, EPS, v, w);
}

/*
 * rmsd  entre X:3xN et Y:3xN stockés dans des vecteurs
 * par colonnes mises bout à bout.
 * effet de bord: centrage de X et Y
 */

void frmsd(DtPoint3 *X,  DtPoint3 *Y, int n, double *r) {
  double Xm[3], Ym[3];
  double x, y, sum, v, w[4];
  int i,j,k;

  euler(X, Y, n, &v, w);

  for (i=0; i<3;i++) {
    sum=0;
    for (j=0; j<n;j++)
      sum+=X[j][i]; 
    Xm[i]=sum/n;
  }
  for (i=0; i<3;i++) {
    sum=0;
    for (j=0; j<n;j++)
      sum+=Y[j][i]; 
    Ym[i]=sum/n;
  }

  sum=0;
  for (i=0; i<3; i++) 
    for (j=0; j<n; j++) {
      x=X[j][i]-Xm[i];
      y=Y[j][i]-Ym[i];
      sum+=x*x+y*y;
    }
  *r=sqrt(fabs(sum-2*v)/n);
}

