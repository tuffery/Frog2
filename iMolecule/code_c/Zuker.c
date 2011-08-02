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

/*
 * This program is for 3D superposition of sets of coordinates.
 * It is based on the quaternions, as proposed by:
 * Zuker & Somorjai, Bulletin of Mathematical Biology,
 * vol. 51, No 1, p 55-78, 1989.
 * 
 * The search for the largest eigen value is based on a QR, and not
 * on the Jacobian based solvation approach used by most programmers.
 * 
 * In our experience, the routine is very robust.
 *
 * If you use this program standalone, the call should be on the form:
 * QBestFit <-bfxyz, -bfpdb> < crdfile.xyz > fitcrd.xyz
 * 
 * Input:
 *   first line: N, the number of coordinates
 *   following 2 x N lines: N coordinates (x y z, e.g. "0.134 2.345 12.3245", i.e. no commas)
 *   the first N xyz lines will be the template, the last N xyz lines will 
 *   undergo the best fit transformation.
 * 
 * Output:
 *   first line: original rmsd, before fit.
 *   second line : bestfit rmsd.
 *   nextt 4 lines: the 4x4 transformation matrix
 *   N following lines: the N coordinates transformed.
 *
 * If you want it directly, the main function is:
 * 
 * zuker_superpose(DtPoint3 *c1, DtPoint3 *c2, int len, DtMatrix4x4 M)
 *
 * This program can be freely used, modified and distributed, given that this 
 * header is maintained.
 *
 * If using this routine for publication, please cite F. Guyon and P. Tuffery.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define DtFloat		double	        /* Le type flottant.	        */
typedef DtFloat         DtPoint3[3];
typedef DtFloat         DtPoint4[4];
typedef DtFloat         DtMatrix3x3[3][3];
typedef DtFloat         DtMatrix4x4[4][4];


#define boucle(i,min,max) for((i)=(min);(i)<(max);(i)++)
#define max(a,b) (a<b?b:a)
#define min(a,b) (a<b?a:b)
#define sign(a,b) ((b)>=0.0 ? fabs(a) : -fabs(a))
#define EPS 1.e-10


double **alloc_mat(int n, int m) {
  int i;
  double **mat = (double **)calloc(m,sizeof(double *));
  for (i=0; i<n; i++)
    mat[i]=(double *)calloc(n, sizeof(double));
  return mat;
}

void free_mat(double **mat, int n){
  int i;
  for (i=0; i<n; i++)
    free(mat[i]);
  free(mat);
}

void print_matrice(double **a, int n, int m, char *name){
  int i,j;

  for (i=0; i<n; i++) {
    for (j=0; j<m; j++) 
      printf ("%s[%d,%d]=%lf ", name, i, j, a[i][j]);
    printf ("\n");
  }
}

void print_vector(double *a, int n, char *name){
  int i;

  for (i=0; i<n; i++) 
      printf ("%s[%d]=%lf\n ", name, i, a[i]);
}


double *random_vect(int n) {
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

static double *product(double **A, double *x, int n) {
  int i, j;
  double sum;
  double *y=(double *)calloc(n, sizeof(double));

  for (i=0; i<n; i++) {
    sum=0;
    for (j=0; j<n; j++)
      sum+=A[i][j]*x[j];
    y[i]=sum;
  }
  return y;
}

static int lu_c (double a[4][4],  int n)
{
 int i,j,k,err;
 double pivot,coef,s;

 err=1;
 k=0;
 while (err==1 && k<n) {
  pivot=a[k][k];
  if(fabs(pivot)>=EPS) {
    for(i=k+1;i<n;i++) {
      coef=a[i][k]/pivot;
      for(j=k;j<n;j++)
	a[i][j] -= coef*a[k][j];
      a[i][k]=coef;
    }
  }
  else err=0;
  k++;
 }
 if(a[n-1][n-1]==0) err=0;
 return err;
}


static void resol_lu(double a[4][4], double *b, int n)
{
 int i,j;
 double sum;
 double y[n];
 y[0]=b[0];
 for(i=1;i<n;i++) {
  sum=b[i];
  for(j=0;j<i;j++)
    sum-=a[i][j]*y[j];
  y[i]=sum;
 }
 b[n-1]=y[n-1]/a[n-1][n-1];
 for(i=n-1;i>=0;i--) {
  sum=y[i];
  for(j=i+1;j<n;j++)
     sum-=a[i][j]*b[j];
  b[i]=sum/a[i][i];
 }
}

int power(double *a[], int n, int maxiter, double eps, double *v, double *w) {
  int niter,i;
  double *y;
  double sum, l, normy, d;
  y=random_vect(n);
  niter=0;
  do {
    normy=sqrt(inner(y,y,n));
    for (i=0; i<n; i++) w[i]=y[i]/normy;
    y=product(a, w, n);
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

double best_shift(double *a[], int n) {
  double m, M, s;
  double t, sum;
  int i, j;
  t=a[0][0];
  for (i=1; i<n; i++) t=max(t, a[i][i]);
  M=t;
  t=a[0][0];
  for (i=0; i<n; i++) {
    for (sum=0,j=0; j<n; j++)
      if (j!=i) sum+=fabs(a[i][j]);
    t=min(t, a[i][i]-sum);
  }
  m=t;
  s=-0.5*(M+m);
  for (i=0; i<n; i++) 
    a[i][i]=a[i][i]+s;
  return s;
}

double lmax_estim (double a[4][4], int n) {
  double t, sum;
  int i, j;
  t=a[0][0];
  for (i=0; i<n; i++) {
    for (sum=0,j=0; j<n; j++)
      if (j!=i) sum+=fabs(a[i][j]);
    t=max(t, a[i][i]+sum);
  }
  return t;
}


int shift_power(double *a, int n, int maxiter, double eps, double *v, double *w) {  
  double **tmp;
  double sh;
  int niter;
  int i,j;

  /* fprintf(stderr,"ESSAI avec shift_power.\n"); */

  tmp=alloc_mat(n, n);    
  /* copyMat(a, tmp, n, n); */
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      tmp[i][j] = a[i*n+j];
    }
  }
  sh=best_shift(tmp, n);   

  niter=power(tmp, n, maxiter, eps, v, w);
  *v=*v-sh;
  free_mat(tmp, n);
  return niter;
}

int inverse_power(double a[4][4], int n, int maxiter, double eps, double *v, double *w) {
  int niter,i;
  double *y;
  double r, sum, l, normy, d;
  y=random_vect(n);
  niter=0;
  
  r=lmax_estim(a, n);
  for (i=0; i<n; i++) a[i][i]=a[i][i]-r;
  if (lu_c(a, n)==0) {
    /* fprintf(stderr,"ATTENTION ! cas singulier de inverse_power.\n"); */
    free(y);
    //exit(0);
    return 0;
  }

  do {
    normy=sqrt(inner(y,y,n));
    for (i=0; i<n; i++) {
      w[i]=y[i]/normy;
      y[i]=w[i];
    }
    resol_lu(a, y, n);
    l=inner(w,y,n);
    niter++;
    for (sum=0,i=0; i<n; i++) {
      d=y[i]-l*w[i];
      sum+=d*d;
    }
    d=sqrt(sum);
  } while (d>eps*fabs(l) && niter<maxiter);
  free(y);
  *v=r+1.0/l;
  return niter;
}

int largestEV4(double R[4][4], double v[4], double *vp)
{
  double iw[4];
  double dx,dy,dz, dt;
  int aCycle;
  double M2[4][4];
  int rs;

  memcpy(M2,R,sizeof(double)*16);

  rs = inverse_power(R, 4, 10000, 1.e-8, vp, v);

  if (!rs) 
    return shift_power(&M2[0][0], 4, 10000, 1.e-8, vp, v);

  return rs;

}


/* ==================================================================
 * Calcule la matrice d'inertie des points
 *
 * id      : indices des atomes dans At
 * from,tto: indices extremes
 * P       : le barycentre
 * ==================================================================
 */
DtMatrix3x3 *XYCov(DtMatrix3x3 *pM, 
		   DtPoint3 *X, 
		   DtPoint3 *Y, 
		   DtPoint3 Xmean, 
		   DtPoint3 Ymean, 
		   int aSze)
{
  int i,j;
  double Xx,Xy,Xz;
  double Yx,Yy,Yz;
  double daSze;

/*   fprintf(stdout,"inertia, len %d\n",aSze); */

  /* X average*/
  Xmean[0] = Xmean[1] = Xmean[2] = 0.;
  for (i=0;i<aSze;i++) {
    Xmean[0] += X[i][0];
    Xmean[1] += X[i][1];
    Xmean[2] += X[i][2];
  }
  if (aSze) {
    daSze = (double) aSze;
    Xmean[0] /= daSze;
    Xmean[1] /= daSze;
    Xmean[2] /= daSze;
  }

  /* Y average*/
  Ymean[0] = Ymean[1] = Ymean[2] = 0.;
  for (i=0;i<aSze;i++) {
    Ymean[0] += Y[i][0];
    Ymean[1] += Y[i][1];
    Ymean[2] += Y[i][2];
  }
  if (aSze) {
    daSze = (double) aSze;
    Ymean[0] /= daSze;
    Ymean[1] /= daSze;
    Ymean[2] /= daSze;
  }

  /* Covariance matrix */
  if (pM == NULL) {
    pM = (DtMatrix3x3 *) calloc(1,sizeof(DtMatrix3x3));
  } else {
    memset((void *) pM, 0, sizeof(DtMatrix3x3));
  }

  for (i=0;i<aSze;i++) {
    Xx = (double) X[i][0] - Xmean[0];
    Xy = (double) X[i][1] - Xmean[1];
    Xz = (double) X[i][2] - Xmean[2];
    Yx = (double) Y[i][0] - Ymean[0];
    Yy = (double) Y[i][1] - Ymean[1];
    Yz = (double) Y[i][2] - Ymean[2];

    (*pM)[0][0] += Xx*Yx;
    (*pM)[0][1] += Xx*Yy;
    (*pM)[0][2] += Xx*Yz;

    (*pM)[1][0] += Xy*Yx;
    (*pM)[1][1] += Xy*Yy;
    (*pM)[1][2] += Xy*Yz;

    (*pM)[2][0] += Xz*Yx;
    (*pM)[2][1] += Xz*Yy;
    (*pM)[2][2] += Xz*Yz;
  }

  return pM;
}

/* ------------------------------------------------------------------ 
    Matrice 4x4 Translation
   ---------- PREMULTIPLIE -> Y = XM (X vecteur ligne) ----------- 
   ------------------------------------------------------------------ */
void MkTrnsIIMat4x4(DtMatrix4x4 m, DtPoint3 tr)
{
 m[0][0] = 1.;    m[0][1] = 0.;    m[0][2] = 0.;    m[0][3] = 0.;
 m[1][0] = 0.;    m[1][1] = 1.;    m[1][2] = 0.;    m[1][3] = 0.;
 m[2][0] = 0.;    m[2][1] = 0.;    m[2][2] = 1.;    m[2][3] = 0.;
 m[3][0] = tr[0]; m[3][1] = tr[1]; m[3][2] = tr[2]; m[3][3] = 1.;
}

/* Matrices 4 x 4 (a x b dans c) ------------------------------------ */
void mulMat4x4(DtMatrix4x4 a,DtMatrix4x4 b,DtMatrix4x4 c)
{
  int i,j,k;
  boucle(i,0,4) {
    boucle(j,0,4) {
      c[i][j] = 0.;
      boucle(k,0,4) c[i][j] += a[i][k]*b[k][j];
    }
  }
}


/* ===============================================================
 * ---------- Transformation d'une coordonnee par rmat. ----------
 * =============================================================== */

void crdRotate(DtPoint3 p,DtMatrix4x4 rmat)
{
  double x,y,z;

  x = p[0] * rmat[0][0] + p[1] * rmat[1][0] + p[2] * rmat[2][0] + rmat[3][0]; 
  y = p[0] * rmat[0][1] + p[1] * rmat[1][1] + p[2] * rmat[2][1] + rmat[3][1]; 
  z = p[0] * rmat[0][2] + p[1] * rmat[1][2] + p[2] * rmat[2][2] + rmat[3][2]; 
  p[0] = x;
  p[1] = y;
  p[2] = z;
}

#define SQUARED_DISTANCE(K,R) ((K)[0] - (R)[0]) * ((K)[0] - (R)[0]) + ((K)[1] - (R)[1]) * ((K)[1] - (R)[1]) + ((K)[2] - (R)[2]) * ((K)[2] - (R)[2])

DtFloat squared_distance(DtPoint3 R,DtPoint3 K)
{
  return (DtFloat) ((K[0] - R[0]) * (K[0] - R[0]) +
		    (K[1] - R[1]) * (K[1] - R[1]) +
		    (K[2] - R[2]) * (K[2] - R[2]));
}


/* 
 * Best fit matrix as proposed by:
 * Zuker & Somorjai, Bulletin of Mathematical Biology,
 * vol. 51, No 1, p 55-78, 1989.
 * 
 * c1, c2 two sets of coordinates to superpose.
 * len : the number of coordiantes of c1, c2.
 *
 * M, the 4x4 transforamtion matrix to pass from c2 to c2 superposed onto c1.
 *
 * Input: 
 *   c1, c2, and len must be given. 
 *   M should be passed, but is meaningless.
 * 
 * Output:
 *   c2 is superposed on c1.
 *   M, a transformation matrix ready for use.
 *
 */
double zuker_superpose(DtPoint3 *c1, DtPoint3 *c2, int len, DtMatrix4x4 M)
{
  DtMatrix3x3  C;
  DtMatrix3x3 *pC; 
  DtMatrix4x4  RM;
  DtMatrix4x4  TM;
  DtMatrix4x4  TMP;
  DtMatrix4x4  TX;
  DtMatrix4x4  TY;
  DtMatrix4x4  P;
  DtMatrix4x4  Pd;
  DtPoint4     lambdas;
  DtPoint4     V;
  DtPoint3     bc1, bc2;
  DtPoint3     trx, try;
  
  double eval;
  double squared_rms = 0.;

  int nCycles;
  int aDot;
  int i,j;

//   fprintf(stderr, "zuker_superpose: \n");
//   fprintf(stderr, "%.3lf %.3lf %.3lf\n", c1[len-1][0], c1[len-1][1], c1[len-1][2]); 
//   fprintf(stderr, "%.3lf %.3lf %.3lf\n", c2[len-1][0], c2[len-1][1], c2[len-1][2]); 

  /* Compute transformation matrix as proposed by zuker */
  pC = &C;
  pC = XYCov(pC, (DtPoint3 *) c1, (DtPoint3 *) c2, bc1, bc2, len);

#if 0
  fprintf(stderr, "zuker_superpose: XYCov done ...\n");
  for (i=0;i<3;i++) {
    for (j=0;j<3;j++) {
      fprintf(stderr,"%.2lf ", C[i][j]);
    }
    fprintf(stderr,"\n");
  }
#endif

  P[0][0] = -C[0][0]+C[1][1]-C[2][2];
  P[0][1] = P[1][0] = -C[0][1]-C[1][0];
  P[0][2] = P[2][0] = -C[1][2]-C[2][1];
  P[0][3] = P[3][0] =  C[0][2]-C[2][0];

  P[1][1] = C[0][0]-C[1][1]-C[2][2];
  P[1][2] = P[2][1] = C[0][2]+C[2][0];
  P[1][3] = P[3][1] = C[1][2]-C[2][1];

  P[2][2] = -C[0][0]-C[1][1]+C[2][2];
  P[2][3] = P[3][2] = C[0][1]-C[1][0];

  P[3][3] = C[0][0]+C[1][1]+C[2][2];

#if 0
  printMat4x4("zuker P", P);
#endif

  // fprintf(stderr, "zuker_superpose: Will largestEV4 ...\n");
  nCycles = largestEV4(P, V, &eval);
  // fprintf(stderr, "zuker_superpose: largestEV4 done ...\n");

  RM[0][0] = -V[0]*V[0]+V[1]*V[1]-V[2]*V[2]+V[3]*V[3];
  RM[1][0] =  2*(V[2]*V[3]-V[0]*V[1]);
  RM[2][0] =  2*(V[1]*V[2]+V[0]*V[3]);
  RM[3][0] =  0.;

  RM[0][1] = -2*(V[0]*V[1]+V[2]*V[3]);
  RM[1][1] = V[0]*V[0]-V[1]*V[1]-V[2]*V[2]+V[3]*V[3];
  RM[2][1] =  2*(V[1]*V[3]-V[0]*V[2]);
  RM[3][1] =  0.;

  RM[0][2] =  2*(V[1]*V[2]-V[0]*V[3]);
  RM[1][2] = -2*(V[0]*V[2]+V[1]*V[3]);
  RM[2][2] = -V[0]*V[0]-V[1]*V[1]+V[2]*V[2]+V[3]*V[3];
  RM[3][2] =  0.;

  RM[0][3] =  0.;
  RM[1][3] =  0.;
  RM[2][3] =  0.;
  RM[3][3] =  1.;

  /* printMat4x4("zuker RM", RM); */

  /* Solution eprouvee ! */
  try[0] = - bc2[0]; try[1] = - bc2[1]; try[2] = - bc2[2];
  MkTrnsIIMat4x4(TY, try);
  MkTrnsIIMat4x4(TX, bc1);
  mulMat4x4(TY,RM,TMP);
  mulMat4x4(TMP, TX, M);

  // printf("ZUCKER OK 1\n");

  /* Now superpose the coordinates */
  for (aDot=0;aDot<len;aDot++) {
    crdRotate(c2[aDot],M);
  }

  /* Compute squared RMSd */
  for (aDot=0;aDot<len;aDot++) {
    squared_rms += squared_distance(c1[aDot],c2[aDot]);
  }


  // printf("ZUCKER: Found %.2lf %d\n", squared_rms, len);

  return squared_rms / (double) len;
}

